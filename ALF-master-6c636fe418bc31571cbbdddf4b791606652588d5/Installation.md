## Installation.
ALF is pretty self-contained, you only need a LAPACK and a BLAS implementation
as external libraries. For compiling the source code you need make and a
Fortran 2003 compatible Compiler.

### Linux
In the following we give hints on how to install relevant packages.

#### Debian/Ubuntu/Linux Mint
- sudo apt-get install gfortran liblapack-dev make

#### Red Hat/Fedora/CentOS
- sudo yum install gcc-gfortran make liblapack-devel

#### OpenSuSE / SLES
- sudo zypper install gcc-gfortran make lapack-devel

#### Arch Linux
- pacman -S make gcc-fortran lapack


### Other Unixes
gfortran and the lapack implementation from netlib.org should be available for
your system. Consult the documentation of your system on how to install the
relevant packages. The package names from linux should give good starting points
for your search.

### MacOS
gfortran for MacOS can be found at https://gcc.gnu.org/wiki/GFortranBinaries#MacOS.
Detailed information on how to install the package  can be found at: https://gcc.gnu.org/wiki/GFortranBinariesMacOS.
You will need to have Xcode as well as the  Apple developer tools installed. 

### Windows
The easiest way to compile Fortran code in Windows is trough Cygwin, which 
provides a Unix-like environment for Windows. The installer also works as a 
package manager, providing an extensive collection of software from the Unix 
ecosystem. For convenience, we provide a zip archive containing the Cywin
installer and a local repository with all the additional software needed for 
ALF. 

Steps for installing Cygwin:
- Unizip the archive
- Execute "setup-x86_64.exe". If administrator rights are missing, execute it 
from the command line as "setup-x86.exe --no-admin".
- In the setup choose "Install from local directory".
- Choose root directory, where cygwin will be installed. In this directory you 
will also find the home directory of the Unix environment.
- Choose the diretory "cygwin_ALF" as local package Directory.
- At the "Select Packages" screen, in "Categories" view, at the line marked 
"All", click on the word "default" so that it changes to "install".
- Finish installation
- To add, remove or update installed packages, rerun the installer setup and 
chose "Install from Internet".

## Building
After you have obtained the source code and have set up the build environment
you need to build the source files. Executing make in the root directory of ALF
does the Job.
