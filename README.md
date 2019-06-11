# Lattice Tester

_A software package for testing the uniformity of integral lattices in the real space_

## What this software is about

_Lattice Tester_ is a library offering facilities to manipulate integral lattices 
and apply certain algorithms on them, for example to compute a reduced basis,
find a shortest nonzero vector, dualize the lattice, etc.
The lattices considered are in t dimensions and must be specified by giving a 
basis of t vectors with rational coordinates. By multiplying these basis vectors 
by some appropriate integer m, one obtains a lattice whose points have all integer coordinates.  
This is how the lattices are represented in Lattice Tester: all basis vectors 
have integer coordinates.

Lattices that satisfy this property are encountered for example when studying 
lattice rules for quasi-Monte Carlo integration and when studying certain types 
of linear congruential generators whose vectors of successive output values have
a lattice structure. The purpose of Lattice Tester is to compute various 
**figures of merit** that serve a measures of quality for integral lattices. 
These figures of merit can compute measures of unformity for several projections 
of the lattice over lower-dimensional subspaces, and combine them in some way.
Examples of such measures of uniformity include the length of the shortest 
nonzero vector in the dual lattice (the spectral test), the number of hyperplanes 
that contain all the points in the unit hypercube, the beyer quotient, etc.

_Lattice Tester_ was built primarily as a base library for the software packages
[LatNet Builder](https://github.com/umontreal-simul/latbuilder),
and [LatMRG](https://github.com/umontreal-simul/latmrg), designed to analyze
the lattice structure of lattice rules (for quasi-Monte-Carlo) and linear 
congruential random number generators, respectively. 
These packages require shared functionnalities which are provided by _Lattice Tester_.

More details on  *Lattice Tester* can be found in:
- [The tutorial](http://umontreal-simul.github.io/latticetester/df/d1d/examples_page.html) which 
  provides examples of how to use the library.
- [The API documentation](http://umontreal-simul.github.io/latticetester/namespaces.html),
  which specifies the interface.
- [The theoretical background](http://umontreal-simul.github.io/latticetester/da/d18/a_intro.html),
  which gives an overview of the underlying theory.
- [The complete user's guide](http://umontreal-simul.github.io/latticetester/), that contains all of the above.

_Lattice Tester_ is free open source software, distributed under the Apache License.

## Compiling

### Software Dependencies

Compiling *Lattice Tester* requires the following software to be installed:

* [NTL](http://www.shoup.net/ntl/index.html) 10.4.0 or later
* [GMP](https://gmplib.org/) compatible version with your NTL installation
* [Python](https://www.python.org/) *(Needed by waf to build the program)*
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the documentation)*

You will also need a recent compiler compliant with the C++14 standard.

### Configuring the Build

*Lattice Tester* relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source. Waf is included in the *Lattice Tester* source 
tree, but it depends on [Python](http://python.org/download), which must be 
available on the system on which *Lattice Tester* is to be compiled.

The commands below should work verbatim under Linux and MacOS systems.
**Microsoft Windows** users should replace every instance of `./waf` 
with the path under which the Python executable
(`python.exe`) or simply with `python waf`
if the Python installation path is accessible from the system `%PATH%`
environment variable.

Change the current directory to the root directory of the package, for example:

    cd latticetester

if you obtained the source code with the `git` command.
If you obtained the source code from the ZIP archive, the directory should be
named `latticetester-master` instead of `latticetester`.
At the root of the source tree lies the `waf` script, which manages the build
process.

Try:

	./waf --help

to see the various commands and options.

There are 6 options that you might want or need to use:
- `--out /path/to/build/location` allows you to specify in which directory the
  build process will operate. The default is `./build`. You will need permission
  to write in that directory.
- `--prefix /path/to/installation/location` allows you to specify in which 
  directory you would like to install *Lattice Tester* after it's compilation.
  The default is `/usr/local` on Linux (waf's default). You will need permission
  to write in that directory.
- `--ntl /path/to/NTL` allows you to specify the location of your NTL 
  installation. You will only need this flag if waf does not find your NTL
  installation automatically.
- `--gmp /path/to/gmp` allows you to specify the location of your gmp
  installation. You will only need this flag if waf doesn't find your gmp
  installation automatically.
- `--build-docs` waf will build the documentation if this flag is specified and 
  will not build it if it is omitted.
- `--link-static` if this flag is specified, the compiler will link all the 
  libraries statically to the executable program. This might be practical if
  you installed NTL in non standard paths.

First, the project must be configured with:

	./waf configure --needed --flags

For example, if NTL, and GMP are not part of the standard system installation and were
manually installed under, say, the `/opt/ntl`, and `/opt/gmp` directories —
which means that `/opt/ntl` and `/opt/gmp` all contain subdirectories named
`include` and `lib` — the following command indicates `waf` where to find these
two libraries:

        ./waf configure --ntl /opt/ntl --gmp /opt/gmp

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build Lattice Tester, before running the `waf
configure` command.

A simple 
    ./waf configure
command should be enough to configure `waf` for a minimal build,
without documentation. The documentation can be built by
appending the `--build-docs` option to `waf configure`, if
  [Doxygen](http://www.stack.nl/~dimitri/doxygen/) is available on the system.

Errors will be reported if required software components cannot be found.  In
that case, you should check the dependencies installation paths.

If a UNIX shell is available, it is also possible to run the simple `configure.sh`
script with `./configure.sh` to avoid typing the configure command by hand 
(that is especially usefull if you have a few flags to include).

### Building and Installing

Once everything is configured correctly, the following command will build the
*Lattice Tester* library and command-line tool:

    ./waf build

If the build process completed without errors, *Lattice Tester* can be installed to the
directory specified with the `--prefix` option during the configuration step,
with:

    ./waf install


## Running Lattice Tester

A *Lattice Tester* executable can be found in the `bin` subdirectory, under 
the installation prefix.  It includes:

- `lattest`: study lattice properties;

The [user guide](http://umontreal-simul.github.io/latticetester/) provides
further details.

**NOTE:** Under Windows, the programs have an additional `.exe` extension.

Before executing lattest, it may be necessary to inform the dynamic
linker where to find the NTL and GMP shared libraries.  Under Linux
this is done by appending the paths to the `LD_LIBRARY_PATH` environment
variable, e.g.,

    export LD_LIBRARY_PATH=/opt/ntl/lib:/opt/gmp/lib:$LD_LIBRARY_PATH

with a Bash-compatible shell.

**Microsoft Windows** users might need to copy the NTL and GMP DLLs into the
same directory (`$HOME/latticetestersoft/bin`, for example) as the executable programs.
