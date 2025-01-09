# PCT

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/RTKConsortium/PCT/blob/master/LICENSE)

PCT (Proton Computed Tomography) is a toolkit used to process proton CT data and reconstruct proton stopping power maps. PCT is written both in C++ and Python, and is designed to be used either as a code library, or with command-line applications.

## Installation

PCT can be installed from source from its [git repository](https://github.com/SimonRit/PCT) hosted at GitHub. The command to clone the repository is:
```bash
git clone https://github.com/RTKConsortium/PCT.git
```

Configuration and compilation is handled via [CMake](https://cmake.org/). Configuration can be achieved using
```bash
cmake path/to/source/folder  # where PCT was cloned
```
or using a graphical tool such as `ccmake`. Once configuration completes, the compilation can be achieved using
```bash
make
```
Keep in mind that `make` can be sped-up using the `-j` option, allowing to parallelize the compilation, for instance using `make -j 4` if your computer has 4 cores.

Optionally, PCT can be added to the user's `$PATH` variable using for instance
```bash
export PATH=path/to/build/folder:${PATH}  # where PCT was built
```
in the user's `.bashrc` file.

Some functionalities of PCT require [GATE](http://www.opengatecollaboration.org/) to run. [Instructions on how to install GATE can be found in the GATE documentation.](https://opengate-python.readthedocs.io/en/master/user_guide/user_guide_installation.html) An example of proton CT GATE simulation can be found [here](../gate/protonct.py).

## Usage

Usage of each PCT application can be described using the `--help`/`-h` option. For instance, running
```bash
pctfdk --help
```
displays the help for `pctfdk`.

Reconstruction typically involves the following steps:
- `pctpairprotons` in order to arrange ROOT data in a format described [here](pct_format.md).
- `pctpaircuts` in order to remove nuclear collisions.
- `pctbinning` in order to compute the distance-driven binning as described [here]( https://doi.org/10.1118/1.4789589).
- `pctfdk` in order to reconstruct the data generated in the previous step using distance-driven FDK (as described [here](https://doi.org/10.1118/1.4789589)).

## Copyright Centre National de la Recherche Scientifique

Licensed under the Apache License, Version 2.0 (the "[License](https://www.apache.org/licenses/LICENSE-2.0.txt)"); you may not use this file except in compliance with the [License](https://www.apache.org/licenses/LICENSE-2.0.txt).

Unless required by applicable law or agreed to in writing, software distributed under the [License](https://www.apache.org/licenses/LICENSE-2.0.txt) is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the [License](https://www.apache.org/licenses/LICENSE-2.0.txt) for the specific language governing permissions and limitations under the [License](https://www.apache.org/licenses/LICENSE-2.0.txt).