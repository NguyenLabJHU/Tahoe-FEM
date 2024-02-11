#!/bin/sh

## Please note this file was created in Windows Subsystem for Linux so it might have DOS format
## If UNIX conversion is required, then install (Ubuntu): sudo apt install dos2unix


## Before running the script, open the bash terminal and create a directory like the following
## mkdir /home/user-name/TAHOE
## Copy this file to the above directory and change the directory name below as you created
## Convert to UNIX format: dos2unix install_pre_ubuntu.sh
## To give permission: chmod a+x install_pre_ubuntu.sh
## To run this script: ./install_pre_ubuntu.sh


## Tahoe main directory: make sure to change this here with your directory

TAHOE_MAIN="/home/user-name/TAHOE"
cd $TAHOE_MAIN


## Installing pre-requisite packages 
## Ubuntu already comes with some of these packages installed
sudo apt update && sudo apt upgrade
# install compilers, perl
sudo apt install gcc -y
sudo apt install g++ -y
sudo apt install gfortran -y
sudo apt install perl
# intall python and python packages
sudo apt install python3
sudo apt install python-pip -y
sudo apt install software-properties-common -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt update
# install FEniCS (python)
sudo add-apt-repository ppa:fenics-packages/fenics 
sudo apt update
sudo apt install fenics -y
# mpi, ssh, and other development tools
sudo apt install openmpi-bin openmpi-common libopenmpi-dev openssh-client openssh-server -y
sudo apt install build-essential -y
sudo apt install xutils-dev -y
# install git and cvs
sudo spt install git
sudo apt install cvs -y
sudo apt dist-upgrade


## Installing SEACAS

echo MAIN TAHOE INSTALLATION DIRECTORY IS: $PWD
git clone https://github.com/sandialabs/seacas.git
cd seacas
./install-tpl.sh
mkdir build && cd build
# use cmake-exodus instead of cmake-config for the exodus file format to be used in TAHOE
../cmake-exodus
make && make install



## Add environment variables to bash profile
echo export CVS_RSH=ssh >> ~/.bashrc
echo export TAHOE_MAIN=$TAHOE_MAIN >> ~/.bashrc
echo export ACCESS=$TAHOE_MAIN/seacas >> ~/.bashrc
echo export LD_LIBRARY_PATH=$ACCESS/lib >> ~/.bashrc

source ~/.bashrc

# Soft link needs to be created (only for first-time compilation)
# It's needed for the modules that depend on SEACAS
# If you make changes to any of those modules then this will be required
cd $ACCESS
## this might have been already done
ln -s include/ inc
cd lib
ln -s libexodus.dylib libexoIIv2c.dylib
ln -s libexodus.a libexoIIv2c.a
