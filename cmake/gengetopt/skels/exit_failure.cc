/*
 * File automatically generated by
 * gengen 1.4.1 by Lorenzo Bettini 
 * http://www.gnu.org/software/gengen
 */

#include "exit_failure.h"

void
exit_failure_gen_class::generate_exit_failure(ostream &stream, unsigned int indent)
{
  string indent_str (indent, ' ');
  indent = 0;

  if (handle_error)
    {
      stream << "if (result == EXIT_FAILURE)";
      stream << "\n";
      stream << indent_str;
      stream << "  {";
      stream << "\n";
      stream << indent_str;
      indent = 4;
      stream << "    ";
      generate_string (parser_name, stream, indent + indent_str.length ());
      stream << "_free (args_info);";
      indent = 0;
      stream << "\n";
      stream << indent_str;
      stream << "    exit (EXIT_FAILURE);";
      stream << "\n";
      stream << indent_str;
      stream << "  }";
      stream << "\n";
      stream << indent_str;
      stream << "\n";
      stream << indent_str;
    }
  stream << "return result;";
}
