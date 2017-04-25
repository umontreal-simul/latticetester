# LatCommon Library
*Tools for studying lattice properties*

LatCommon regroups the tools that are used by the software packages
[Lattice Builder](https://github.com/umontreal-simul/latbuilder)
and
[LatMRG](https://github.com/umontreal-simul/latmrg).
LatCommon is already embedded in these and does not need to be installed
separately.


## Compiling the Source Code


### Software Dependencies

Compiling LatCommon requires the following software to be installed on
the system:

* [Python](http://python.org/) 2.7 or later
* [Boost C++ Libraries](http://www.boost.org/) 1.57.0 or later
  <small>
  Installation instructions (**section 5 is important**):
  [for Linux / MacOS](http://www.boost.org/doc/libs/release/more/getting_started/unix-variants.html),
  [for Microsoft Windows](http://www.boost.org/doc/libs/release/more/getting_started/windows.html)
  </small>
* [NTL](http://shoup.net/ntl/) 6.2.1 or later *(optional)*
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the documentation)*

**MacOS** users can install these dependencies through
[MacPorts](http://www.macports.org/) by [installing the MacPorts
software](http://www.macports.org/install.php), then by installing the
following packages:
[python27](https://trac.macports.org/browser/trunk/dports/lang/python27/Portfile),
[boost](https://trac.macports.org/browser/trunk/dports/devel/boost/Portfile),
[git](https://trac.macports.org/browser/trunk/dports/devel/git/Portfile) and
[doxygen](https://trac.macports.org/browser/trunk/dports/textproc/doxygen/Portfile).

You will also need a recent compiler compliant with the C++11 standard.


### Obtaining the Source Code

Get the latest source tree from GitHub, either by typing:

	git clone https://github.com/umontreal-simul/latcommon.git

If [Git](http://git-scm.com/) is not available on your system, you can click on
the [Download ZIP](https://github.com/umontreal-simul/latcommon/archive/master.zip)
link from the [LatCommon GitHub page](https://github.com/umontreal-simul/latcommon),
then by unzipping the downloaded archive.

### Configuring the Build

LatCommon relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source.
Waf is included in the LatCommon source tree, but it depends on
[Python](http://python.org/download), which must be available on the system
on which LatCommon is to be compiled.

The commands below should work verbatim under Linux and MacOS systems.
**Microsoft Windows** users should replace every instance of `./waf` 
with `C:\Python27\python waf`, assuming that the Python executable
(`python.exe`) was installed under `C:\Python27`, or simply with `python waf`
if the Python installation path is accessible from the system `%PATH%`
environment variable.

Change the current directory to the root directory of the package, for example:

	cd latcommon

if you obtained the source code with the `git` command.
If you obtained the source code from the ZIP archive, the directory should be
named `latcommon-master` instead of `latcommon`.
At the root of the source tree lies the `waf` script, manages the build
process.
Try:

	./waf --help

to see the various commands and options.
The most relevant options include `--out` to specify the directory in which the
files created during the build process will be placed, `--prefix` to specify
the directory under which you wish to install the LatCommon software
**once it is compiled**, and `--boost` and `--ntl` to specify the directories
under which Boost and NTL were installed, if not under standard system
directories.  First, the project must be configured with:

	./waf configure --prefix $HOME/latsoft

with `$HOME/latsoft` replaced with the directory into which you wish to install
LatCommon.
Here, `$HOME` will expand to your own home directory; you can specify any other
directory to which you have permissions for write access, e.g., with `--prefix
/tmp/latsoft`.

If Boost and NTL are not part of the standard system installation and were
manually installed under, say, the `/opt/boost` and `/opt/ntl` directories —
which means that `/opt/boost` and `/opt/ntl` both contain subdirectories named
`include` and `lib` — the following command indicates `waf` where to find these
two libraries:

	./waf configure --prefix $HOME/latsoft --boost /opt/boost --ntl /opt/ntl configure --link-static

The `--link-static` option suggested above will cause the Boost and NTL
libraries to be linked statically to the executable program, which may be
desirable especially if these are not installed in standard locations.

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build LatCommon, before running the `waf
configure` command.

The above `waf configure` commands configures `waf` for a minimal build,
without documentation nor code examples.  These can be built by
appending the following options to `waf configure`:

* `--build-docs` to generate the documentation, if
  [Doxygen](http://www.stack.nl/~dimitri/doxygen/) is available on the system.
* `--build-examples` to compile and install example code, including
  code from the tutorial, which will also be verified to yield correct output.
  The expected outputs are stored in text files with names matching those of
  programs, under the `examples/tutorial/output` subdirectory.

Errors will be reported if required software components cannot be found.  In
that case, you should check the Boost and NTL installation paths.

### Building and Installing

Once everything is configured correctly, the following command will build the
LatCommon library and command-line tool:

	./waf build

If the build process completed without errors, LatCommon can be installed
to `$HOME/latsoft`, or any directory specified with the `--prefix` options
during the configuration step, with:

	./waf install

