# LatticeTester

_A software package for testing the uniformity of integral lattices in the real space_

## Documentation

_LatticeTester_ is a library that is intended to facilitate the manipulation
of lattices and the execution of certain algorithms on such objects. Lattices
are the discrete subgroup of a vector space obtained by taking all the integer
linear combinations of a set of vectors. They have a lot of structure and can be
found, for example, in certain algorithmicaly generated point sets such as
some types of random number generators and quasi-Monte-Carlo integration point
sets.

Studying those point sets is mostly about trying to mesure their uniformity in
the unit hypercube with what is called **figures of merit**. Figures of merit
are mesures built to capture the uniformity of point sets in a number that can
easily be compared for different point sets. Namely, *LatticeTester* implements
the computation of the **Spectral Test** and of the **Beyer quotient**.

To perform such computations, _LatticeTester_ implements numerous algorithms
to manipulate lattices. Notable ones include a Branch-and-Bound procedure
that computes the **shortest non-zero vector** in a lattice and various **basis
reduction algorithms** using criterions like LLL and BKZ.
In *LatticeTester*, and in general, lattices are represented by a **basis**, a
set of vectors of minimal cardinality generating the entire lattice. Basis
reduction algorithms are algorithms that take the basis of a lattice as input
and try to find another basis of the same lattice, but with vector of shorter
length.

Finally _LatticeTester_ also implements lower level operations on lattices.
Those include the extraction of a basis for a lattice from a set of vectors
generating this lattice, norm computation in both L1 and L2 norms and vector 
swapping in the basis.
Concurrently to the basis computations, _LatticeTester_ can also compute what is
called the **dual basis**, a set of vectors spanning all vectors normal to
hyperplans in the lattice.

_LatticeTester_ is primarily built as a library for the software packages
[LatNet Builder](https://github.com/umontreal-simul/latbuilder),
and [LatMRG](https://github.com/umontreal-simul/latmrg) that
respectively study quasi-Monte-Carlo point sets and congruential random number
generators. The applications targeted by these programs require specific shared
functionnalities and _LatticeTester_ aims at provinding them.

For more details about *LatticeTester* can be found in:
- [The tutorial] containing a few examples showing of the main classes of the
  library
- [The API documentation], for specific implementation details of each class in
  the library
- [The theoretical background], should you need to remind yourself details of the above
- [The full documentation](http://umontreal-simul.github.io/latticetester/) Containing all the above and more


This software package allows it's user to calculate various figures of merit 
such as the normalized spectral test and the Beyer quotient on integral lattices.
It also contains classes that can represent such a lattice and that can be used
to manipulate and reduce it's basis. On the other hand, this program cannot
compute the basis of a lattice by itself. To do this, depending on your use case,
you will either need to use external software or to use some of the other
software available at [https://github.com/umontreal-simul](https://github.com/umontreal-simul),
namely [LatNet Builder](https://github.com/umontreal-simul/latbuilder) for QMC
point sets and [LatMRG](https://github.com/umontreal-simul/latmrg) for PRNG
point sets. A detailed documentation of *LatticeTester* is available 
[here].

## Compiling

### Software Dependencies

Compiling *LatticeTester* requires the following softwares to be installed on
the system:

* [NTL](http://www.shoup.net/ntl/index.html) 10.4.0 or later
* [GMP](https://gmplib.org/) compatible version with your NTL installation
* [Python](https://www.python.org/) *(Needed by waf to build the program)*
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the documentation)*

You will also need a recent compiler compliant with the C++14 standard.

### Configuring the Build

*LatticeTester* relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source. Waf is included in the *LatticeTester* source 
tree, but it depends on [Python](http://python.org/download), which must be 
available on the system on which *LatticeTester* is to be compiled.

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
At the root of the source tree lies the `waf` script, manages the build
process.

Try:

	./waf --help

to see the various commands and options.

There are 6 options that you might want/need to use:
- `--out /path/to/build/location` allows you to specify in which directory the
  build process will operate. The default is `./build`. You will need permission
  to write in that directory.
- `--prefix /path/to/installation/location` allows you to specify in which 
  directory you would like to install *LatticeTester* after it's compilation.
  The default is `/usr/local` on Linux (waf's default). You will need permission
  to write in that directory.
- `--ntl /path/to/NTL` allows you to specify the location of your NTL 
  installation. You will only need this flag if waf doesn't find your NTL
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
C++ compiler to be used to build LatticeTester, before running the `waf
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
*LatticeTester* library and command-line tool:

    ./waf build

If the build process completed without errors, *LatticeTester* can be installed to the
directory specified with the `--prefix` option during the configuration step,
with:

    ./waf install


## Running LatticeTester

The *LatticeTester* executable can be found in the `bin` subdirectory, under 
the installation prefix. These include:

- `lattest`: study lattice properties;

Refer to the [user guide](http://umontreal-simul.github.io/latticetester/) for 
further detail.

**NOTE:** Under Windows, the programs have an additional `.exe` extension.

Before executing the lattest program, it may be necessary to inform the dynamic
linker where to find the NTL and GMP shared libraries.  Under Linux
this is done by appending the paths to the `LD_LIBRARY_PATH` environment
variable, e.g.,

    export LD_LIBRARY_PATH=/opt/ntl/lib:/opt/gmp/lib:$LD_LIBRARY_PATH

with a Bash-compatible shell.

**Microsoft Windows** users might need to copy the NTL and GMP DLLs into the
same directory (`$HOME/latticetestersoft/bin`, for example) as the executable programs.
