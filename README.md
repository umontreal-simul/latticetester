# LatticeTester
_A software package for testing the uniformity of integral lattices in the real space_

## Compiling

### Software Dependencies

Compiling *LatticeTester* requires the following softwares to be installed on
the system:

* [Python](https://conda.io/docs/user-guide/install/download.html)
* [Boost C++ Libraries](http://www.boost.org/) 1.57.0 or later
  <small>
  Installation instructions (**section 5 is important**):
  [for Linux / MacOS](http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html),
  [for Microsoft Windows](http://www.boost.org/doc/libs/release/more/getting_started/windows.html)
  </small>
* [NTL](http://www.shoup.net/ntl/index.html) 10.4.0 or later
* [GMP](http://www.shoup.net/ntl/index.html) compatible version with your NTL installation
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the documentation)*

You will also need a recent compiler compliant with the C++14 standard.

### Configuring the Build

*LatticeTester* relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source.
Waf is included in the *LatticeTester* source tree, but it depends on
[Python](http://python.org/download), which must be available on the system
on which *LatticeTester* is to be compiled.

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
The most relevant options include `--out` to specify the directory in which the
files created during the build process will be placed, `--prefix` to specify
the directory under which you wish to install the LatNet Builder software
**once it is compiled**, and `--boost`, `--ntl` and `--gmp` to specify the directories
under which Boost, NTL, GMP and FFTW were installed, if not under standard system
directories.  First, the project must be configured with:

	./waf configure --prefix $HOME/latticetestersoft

with `$HOME/latticetestersoft` replaced with the directory into which you wish to install
*LatticeTester*.
Here, `$HOME` will expand to your own home directory; you can specify any other
directory to which you have permissions for write access, e.g., with `--prefix
/tmp/latticetestersoft`.

If Boost, NTL, and GMP are not part of the standard system installation and were
manually installed under, say, the `/opt/boost`, `/opt/ntl`, and `/opt/gmp` directories —
which means that `/opt/boost`, `/opt/ntl` and `/opt/gmp` and both contain subdirectories named
`include` and `lib` — the following command indicates `waf` where to find these
two libraries:

  ./waf configure --prefix $HOME/latnetsoft --boost /opt/boost configure --link-static

The following options can also be added to `./waf configure`:

- `--out`: directory in which the files created during the build process will
  be placed;
- `--prefix`: directory under which the *LatticeTester* software will be installed after
  compilation;
- `--boost`: directory under which the Boost libraries are installed;
- `--ntl`: directory under which the NTL library is installed;


The `--link-static` option suggested above will cause the 
libraries to be linked statically to the executable program, which may be
desirable especially if these are not installed in standard locations.

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build LatNet Builder, before running the `waf
configure` command.

The above `waf configure` commands configures `waf` for a minimal build,
without documentation.  The documentation can be built by
appending the `--build-docs` option to `waf configure`, if
  [Doxygen](http://www.stack.nl/~dimitri/doxygen/) is available on the system.

Errors will be reported if required software components cannot be found.  In
that case, you should check the dependencies installation paths.


### Building and Installing

Once everything is configured correctly, the following command will build the
*LatticeTester* library and command-line tool:

    ./waf build

If the build process completed without errors, *LatticeTester* can be installed to the
directory specified with the `--prefix` option during the configuration step,
with:

    ./waf install


## Running LatticeTester

The *LatticeTester* executable can be found in the `bin` subdirectory, under the installation prefix.
These include:

- `lattest`: study lattice properties;

Refer to the user guide (`http://umontreal-simul.github.io/latticetester/`) for further detail.

**NOTE:** Under Windows, the programs have an additional `.exe` extension.

Before executing the lattest program, it may be necessary to inform the dynamic
linker where to find the Boost, NTL and GMP shared libraries.  Under Linux
this is done by appending the paths to the `LD_LIBRARY_PATH` environment
variable, e.g.,

    export LD_LIBRARY_PATH=/opt/boost/lib:/opt/ntl/lib:/opt/gmp/lib:$LD_LIBRARY_PATH

with a Bash-compatible shell.

**Microsoft Windows** users might need to copy the Boost, NTL and GMP DLLs into the
same directory (`$HOME/latticetestersoft/bin`, for example) as the executable programs.
