# Tahoe 

Tahoe is an open-source research-oriented software platform for the development of numerical methods and material models. The goal of the work surrounding Tahoe is the simulation of materials physics involving measures such as stress, deformation, velocity, temperature, or other state variables of interest, for situations that cannot be treated by standard continuum simulation techniques. These situations include material fracture or failure, interfacial adhesion and debonding, shear banding, length-scale dependent elasticity, and plasticity, deformation in small-scale structures, solid-liquid-gas interactions, and other multi-physical phenomena observed at multiple length and time scales for a wide range of materials. Aside from a collection of standard finite elements, Tahoe includes meshfree simulation capability (Reproducing Kernel Particle Method (RKPM)) and other particle methods, such as ellipsoidal and poly-ellipsoidal Discrete Element Method (DEM), Smoothed Particle Hydrodynamics (SPH), PeriDynamics (PD), coupled DEM-SPH, coupled DEM-PD, poly-ellipsoidal DEM with particle fracture, and coupled ellipsoidal Discrete Element Method - Computational Fluid Dynamics (DEM-CFD). Tahoe also includes a number of "cohesive" approaches for modeling fracture. These include both surface and bulk constitutive models that incorporate cohesive behavior. Tahoe is capable of performing static and transient dynamic coupled-physics analysis in two and three dimensions, along with contact. Many capabilities support parallel execution.

This is a custom fork of Tahoe which was originally hosted on [SourceForge](https://tahoe.sourceforge.net). The development branch of Tahoe called `tahoe-dev` on SourceForge has been merged into the main `Tahoe` directory of this fork. However, the input files in the branch tahoe-dev-benchmark are not included in the repository because of the size restriction of GitHub. Interested users can obtain the files from SourceForge. For any general questions about using Tahoe, please refer to the [user's guide](guide_user/user_guide.pdf).


## Notes for Developers

If you are a developer who is interested in contributing to the Tahoe project, then fork this repository and create a development branch for yourself. Perform all the developments in that branch and once you are sure about making the changes available, commit your changes, and send a pull request. We will merge your request to the main branch after review. Learn more about the [GitHub flow process](https://docs.github.com/en/get-started/using-github/github-flow).


## Installing Tahoe

Tahoe is available for installation on macOS and Linux but not natively on Windows. SourceForge-based installation was tested on macOS Sonoma 14.3 and Ubuntu 22.04 LTS. The procedure to obtain the prerequisite on macOS and Linux is slightly different, however, the steps to build Tahoe are the same on both platforms.

> [!NOTE]
> Tahoe can be installed on Windows using [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). We recommend using Ubuntu as the Linux distribution on WSL. Once you have WSL installed, you can follow the procedure for Linux-based installation of Tahoe.

> [!WARNING]
>  This installation guide is created based on the old documentation. GitHub-based installation and building of Tahoe have not been tested yet. Soon it will be tested and updated. Experienced users can still use the guide.

### To-do list
- [ ] Updating the makefile.template to install the parallel version of Tahoe. The following line needs to be changed from
     `"AX_LFLAGS  = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf"` to `"AX_LFLAGS  = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf -lhdf5"`.


### Obtaining prerequisites on macOS

- In case you don't have a zsh profile set up on macOS, open a profile from the terminal. We will use this file later to save environment variables and paths. 
```bash
open ~/.zshrc
```
Save the file and close it. 
- Apple's developer platform is known as Xcode which includes a lot of packages for program development. While there is a GUI version of it, we will need the command line interface only for Tahoe.
```bash
 xcode-select --install
```
Xcode will install the tools in `/usr/bin` directory. Navigate to the directory using `cd` command and use `ls` command there to check the tools installed by Xcode. Among others, you can see g++ and gcc are available but gfortran is not there.
- Download and install homebrew on macOS
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
- Give read and write permission to homebrew
```bash
sudo chown -R $(whoami) $(brew --prefix)/*
```
-  Apple Clang compiler in the Xcode does not include gfortran. Use homebrew to install the following packages including gcc for gfortran.
```bash
brew install --cask xquartz
brew install gcc
brew install wget
brew install cmake
brew install automake
brew install open-mpi
```
homebrew will install the compilers in `/opt/homebrew/bin/` directory. Check the installed version of the compilers by navigating to the directory and using `ls` command before you do the next step 
- Create symlinks of the installed compilers (gcc, g++, and gfortran). At this time, the stable versions of the installed GNU compilers are 13. Following commands will the gcc and g++ from clang and gfortran from homebrew installation.
```bash
ln -s gcc-13 gcc
ln -s g++-13 g++
ln -s gfortran-13 gfortran
```
- Log out of the Zsh terminal (close) to activate these and log back in (reopen).




###  Obtaining prerequisites on Linux (Ubuntu)


Download and run the bash script `install_pre_ubuntu.sh` included in the Tahoe repository to install all the prerequisites. Please make sure to update the **installation directory** in the script. You may have to do a `dos2unix` conversion for the script for which the instructions are included in the script. To run the script, use 
```bash
./install_pre_ubuntu.sh
```

Alternatively, you can copy and paste those commands inside the script one by one to install the necessary packages on Linux.



### Downloading and building Tahoe (macOS and Ubuntu) with necessary modules

- Clone Tahoe from this repository using
```bash
git clone https://github.com/NguyenLabJHU/Tahoe-FEM.git
```
- Now navigate to `Tahoe-FEM` directory using `cd` command
```bash
cd Tahoe-FEM
```
- Now clone and install the SEACAS package with the third-party libraries (TPL).
> [!WARNING]
> Make sure to check the [SEACAS GitHub page](https://github.com/sandialabs/seacas) for updated instructions for building SEACAS library
-  At the time of writing this instruction for Tahoe, the following commands are needed. Use `cmake-exodus` instead of `cmake-config` for the exodus file format to be used in `Tahoe` when building SEACAS. 
```bash
git clone https://github.com/sandialabs/seacas.git
cd seacas
./install-tpl.sh
mkdir build && cd build
../cmake-exodus
make && make install
```
- Go back to the main `Tahoe-FEM` directory and add the following directories as system paths to your zsh profile.
```bash
TAHOE_MAIN=$PWD
echo export TAHOE_MAIN=$TAHOE_MAIN >> ~/.zshrc
echo export ACCESS=$TAHOE_MAIN/seacas >> ~/.zshrc
echo export LD_LIBRARY_PATH=$ACCESS/lib >> ~/.zshrc
source ~/.zshrc
```
- Create symlinks for the following libraries (only for the first-time installation of SEACAS)
```bash
cd $ACCESS
ln -s include/ inc
cd lib
ln -s libexodus.dylib libexoIIv2c.dylib
ln -s libexodus.a libexoIIv2c.a
```
- Log out of your zsh terminal to activate the changes and log back in (close and reopen). Navigate to `Tahoe-FEM` directory.
- Convert the `tahoe-manager` to an executable.
```bash
sudo chmod 755 ./tahoe-manager
```
- Now run `tahoe-manager` to build Tahoe
```bash
./tahoe-manager init build
```
Select GNU-GCC-13 (for serial installation) as architecture when asked. Confirm the choice of all the main modules shown on the screen. Additionally, for optional modules, select 6 (ACCESS), 8 (benchmark_XML), 10 (contrib), 11 (development), 15 (spooles). For the first time, it will take 10-15 minutes to build all the packages including the Tahoe executable. Subsequent building will be faster. 
> [!NOTE]
> For parallel installation, select  GNU-GCC-MPI-13 as the architecture and add metis and spoolesMPI as the optional modules.
- Add the following directories as your system path.
```bash
echo export TAHOE_MOD=$TAHOE_MAIN/tahoe >> ~/.zshrc
echo export TAHOE_DIR=$TAHOE_MOD/tahoe >> ~/.zshrc
echo export PATH=$PATH:$TAHOE_MOD/bin >> ~/.zshrc
source ~/.zshrc 
```
 > [!NOTE]
 > On macOS with Apple Silicon processors, you may have to add the following to your .zshrc file. We have not found a workaround for it yet. Make sure to use the right path including `user-name` appearing in the following command
```bash
install_name_tool -add_rpath “$ACCESS/lib” /Users/user-name/Tahoe-FEM/tahoe/tahoe
```

