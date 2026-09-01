import re
import argparse
import itk
from itk import PCT as pct
import difflib
import inspect
from typing import Optional

__all__ = ["PCTArgumentParser"]


class PCTHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _format_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = pct.__version__ + "\n\nusage: "
        return super()._format_usage(usage, actions, groups, prefix)


class PCTArgumentParser(argparse.ArgumentParser):
    def __init__(self, description=None, **kwargs):
        super().__init__(description=description, **kwargs)
        self.formatter_class = PCTHelpFormatter
        # allow negative numeric tokens to be treated as values, not options. This mirrors CPython behavior in python 3.14
        self._negative_number_matcher = re.compile(r"-\.?\d")
        self.add_argument("-V", "--version", action="version", version=pct.__version__)

    def required_dests(self):
        return sorted(
            a.dest
            for a in self._actions
            if getattr(a, "required", False)
            and a.dest
            and a.dest not in ("help", "version")
        )

    def build_signature(self) -> inspect.Signature:
        """Build a compact Python signature: only required kwargs + **kwargs."""
        params = [
            inspect.Parameter(name=d, kind=inspect.Parameter.KEYWORD_ONLY)
            for d in self.required_dests()
        ]
        # Add catch-all for the many optional CLI options to keep help inline
        params.append(inspect.Parameter("kwargs", kind=inspect.Parameter.VAR_KEYWORD))
        return inspect.Signature(params)

    def build_usage_examples(self, app_name: Optional[str] = None) -> str:
        """Return a Usage examples block for Python help()."""
        name = app_name or self.prog
        req = self.required_dests()
        shell = name + "(" + " ".join(f"--{d} {d.upper()}" for d in req) + ")"
        py = name + "(" + ", ".join(f"{d}={d.upper()}" for d in req) + ")"
        return f"Usage:\n    • Shell-style: {shell}\n    • Python API:  {py}\n\n"

    def apply_signature(self, func):
        """Apply the built signature to a callable and return it."""
        func.__signature__ = self.build_signature()
        return func

    def parse_args(self, args=None, namespace=None):
        """Parse args with optional single-token comma list support for multi-value options.
        Supported forms:
          --opt A B C      (space separated)
          --opt A,B,C      (single token, comma separated)
        Applies to all nargs="+" options, including string (path) lists.
        """
        multi_valued = {}
        for action in self._actions:
            dest = getattr(action, "dest", None)
            if not dest or dest in ("help", "version"):
                continue
            if getattr(action, "nargs", None) != "+":
                continue
            # Neutralize all types (including str) so comma tokens like "a,b" or
            # "1,2,3" can be split after parsing. The original type is restored below.
            cast = action.type or str
            multi_valued[dest] = cast
            action.type = str
        try:
            namespace = super().parse_args(args, namespace)
        finally:
            for action in self._actions:
                dest = getattr(action, "dest", None)
                if dest in multi_valued:
                    action.type = multi_valued[dest]
        for dest, cast in multi_valued.items():
            val = getattr(namespace, dest, None)
            # Case 1: user supplied a single token containing commas (e.g. "1,2,3" or
            # "a.nrrd,b.nrrd"). Split on commas, strip whitespace, drop empty pieces.
            if isinstance(val, list) and len(val) == 1 and isinstance(val[0], str) and "," in val[0]:
                pieces = [s for s in (p.strip() for p in val[0].split(",")) if s]
                setattr(namespace, dest, [cast(piece) for piece in pieces])
            else:
                # Case 2: normal space-separated form (e.g. "1 2 3"). Just cast every token.
                setattr(namespace, dest, [cast(piece) for piece in val])
        return namespace

    def parse_kwargs(self, func_name: Optional[str] = None, **kwargs):
        """Convert Python kwargs to argv and parse them.
        Lists/tuples for multi-value options are serialized as a single comma token."""
        actions = {
            a.dest: a
            for a in self._actions
            if a.dest and a.dest not in ("help", "version")
        }
        for key in kwargs:
            if key not in actions:
                matches = difflib.get_close_matches(
                    key, actions.keys(), n=3, cutoff=0.5
                )
                name = func_name or self.prog or "function"
                msg = f"{name}() got an unexpected keyword argument '{key}'"
                if matches:
                    msg += f"\nDid you mean: {', '.join(matches)}?"
                else:
                    msg += f"\nValid arguments are: {', '.join(sorted(actions.keys()))}"
                raise TypeError(msg)
        argv = []
        for key, val in kwargs.items():
            action = actions[key]
            flag = next(
                (o for o in action.option_strings if o.startswith("--")),
                action.option_strings[0],
            )
            if isinstance(val, bool):
                if val:
                    argv.append(flag)
            elif isinstance(val, (list, tuple)):
                argv += [flag, ",".join(map(str, val))]
            else:
                argv += [flag, str(val)]
        return self.parse_args(argv)
