# Tahoe 

Tahoe is an open-source research-oriented software platform for the development of numerical methods and material models. The goal of the work surrounding Tahoe is the simulation of materials physics involving measures such as stress, deformation, velocity, temperature, or other state variables of interest, for situations that cannot be treated by standard continuum simulation techniques. These situations include material fracture or failure, interfacial adhesion and debonding, shear banding, length-scale dependent elasticity, and plasticity, deformation in small-scale structures, solid-liquid-gas interactions, and other multi-physical phenomena observed at multiple length and time scales for a wide range of materials. Aside from a collection of standard finite elements, Tahoe includes meshfree simulation capability (Reproducing Kernel Particle Method (RKPM)) and other particle methods, such as ellipsoidal and poly-ellipsoidal Discrete Element Method (DEM), Smoothed Particle Hydrodynamics (SPH), PeriDynamics (PD), coupled DEM-SPH, coupled DEM-PD, poly-ellipsoidal DEM with particle fracture, and coupled ellipsoidal Discrete Element Method - Computational Fluid Dynamics (DEM-CFD). Tahoe also includes a number of "cohesive" approaches for modeling fracture. These include both surface and bulk constitutive models that incorporate cohesive behavior. Tahoe is capable of performing static and transient dynamic coupled-physics analysis in two and three dimensions, along with contact. Many capabilities support parallel execution.

This is a custom fork of Tahoe which was originally hosted on [SourceForge](https://tahoe.sourceforge.net). The development branch of Tahoe called tahoe-dev om SourceForge has been merged to the main directory of this fork, however, the input files in the branch tahoe-dev-benchmark are not included in the repository because of its size restriction of Github. Interested users should obtain the file from SourceForge. For any general questions about using Tahoe, please refer to the [user's guide](guide_user/user_manual.pdf).

## Installing Tahoe

Tahoe is available for installation on macOS and Linux but not natively on Windows. SourceForge-based installation was tested on macOS Sonoma 14.3 and Ubuntu 22.04 LTS. The procedure to obtain the prerequisite on macOS and Linux is slightly different, however, the steps to build Tahoe are the same on both platforms.

> [!NOTE]
> Tahoe can be installed on Windows using [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). We recommend using Ubuntu as the Linux distribution on WSL. Once you have WSL installed, you can follow the procedure for Linux-based installation of Tahoe.

> [!WARNING]
>  This installation guide is created based on the old documentation. GitHub-based installation and building of Tahoe have not been tested yet. Soon it will be tested and updated. Experienced users can still use the guide.

### To-do list
- [ ] Need more clarification on using compilers available via Clang vs. Homebrew.
- [ ] Updating the makefile.template to install the parallel version of Tahoe. The following line needs to be changed from
     `"AX_LFLAGS  = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf"` to `"AX_LFLAGS  = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf -lhdf5"`.
- [ ] Addressing the issue with `libf2c_.a` file


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
- Create symlinks of the installed compilers (gcc, g++, and gfortran). At this time, the stable versions of the installed compilers are 13.
```bash
ln -s gcc-13 gcc
ln -s g++-13 g++
ln -s gfortran-13 gfortran
```
- Log out of the Zsh terminal (close) to activate these and log back in (reopen).
- Create a dedicated directory for Tahoe in your home directory.
```bash
cd && mkdir tahoe-install && cd tahoe-install
```
- Clone and install SEACAS with the third-party libraries (TPL). Check the [SEACAS GitHub page](https://github.com/sandialabs/seacas) for updated instructions. However, please use `cmake-exodus` instead of `cmake-config` for the exodus file format to be used in `Tahoe`.
```bash
git clone https://github.com/sandialabs/seacas.git
cd seacas
./install-tpl.sh
mkdir build && cd build
../cmake-exodus
make && make install
```
- Add the following directories as system paths to your zsh profile.
```bash
TAHOE_MAIN=$PWD
echo export TAHOE_MAIN=$TAHOE_MAIN >> ~/.zshrc
echo export ACCESS=$TAHOE_MAIN/seacas >> ~/.zshrc
echo export LD_LIBRARY_PATH=$ACCESS/lib >> ~/.zshrc
echo export CVS_RSH=ssh >> ~/.zshrc
source ~/.zshrc
```
- Create symlinks for the following libraries (only for the first-time installation of SEACAS)
```bash
cd $ACCESS
ln -s include/ inc
cd lib
ln -s libexodus.dylib libexoIIv2c.dylib
$ ln -s libexodus.a libexoIIv2c.a
```
- Log out of your zsh terminal to activate the changes and log back in (close and reopen).



###  Obtaining prerequisites on Linux (Ubuntu)


Download and run the bash script `install_pre_ubuntu.sh` included in the Tahoe repository to install all the prerequisites. Please make sure to update the **installation directory** in the script. You may have to do a `dos2unix` conversion for the script for which the instructions are included in the script. To run the script, use 
```bash
./install_pre_ubuntu.sh
```

Alternatively, you can copy and paste those commands inside the script one by one to do the installation of the necessary packages on Linux.



### Downloading and building Tahoe (macOS and Ubuntu)

- Clone Tahoe from this repository using
```bash
git clone https://github.com/NguyenLabJHU/Tahoe.git
```
- Navigate to the cloned directory and convert `tahoe-manager` to an executable.
```bash
cd tahoe && sudo chmod 755 ./tahoe-manager
```
- Now run `tahoe-manager` to build Tahoe
```bash
./tahoe-manager build
```
Select GNU-GCC-MPI-9.3 (for parallel installation) or GNU-GCC9.3 (for serial installation) as architecture when asked. Additionally, for optional modules, select 0 (CBLAS), 6 (ACCESS), 8 (benchmark_XML), 10 (contrib), 11 (development), 12 (development_benchmark_xml), 13 (f2c), 14 (metis), 15 (spooles), 16 (spoolesMPI). For the first time, it will take 20-30 minutes to build Tahoe executable. Subsequent building will be faster. 
- Add the following directories as your system path.
```bash
echo export TAHOE_MOD=$TAHOE_MAIN/tahoe >> ~/.zshrc
echo export TAHOE_DIR=$TAHOE_MOD/tahoe >> ~/.zshrc
echo export PATH=$PATH:$TAHOE_MOD/bin >> ~/.zshrc
source ~/.zshrc 
```
- On macOS with Apple Silicon processors, then you may have to add the following to your .zshrc file. We have not found a workaround for it yet. Make sure to choose your `user-name`.
```bash
install_name_tool -add_rpath “$ACCESS/lib” /Users/user-name/tahoe-install/tahoe/tahoe/tahoe
```

## Executing Tahoe

### Serial executation

### Parallel execution

