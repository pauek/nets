#!/bin/bash

# Copyright (c) 2007, Pau Fernández

# Check dependencies
if [ -z /usr/bin/g++ ]; then
   echo "You have to install g++!"
   echo 'In Ubuntu/Debian, type "aptitude install g++"'
   exit 1
fi

if ! [ -r /usr/include/random/uniform.h ]; then
   echo "You have to install Blitz++!" 
   echo 'In Ubuntu, type "apt-get install libblitz0-dev"'
   exit 1
fi

echo "Compiling programs..."
make -k -C graphics
make -k -C networks
echo "Done."

[ -d bin ] || mkdir bin
cp -f networks/G* bin/.

echo "+----------------------------------------------------------+"
echo '|   You will find the executables in the "bin" directory   |'
echo "+----------------------------------------------------------+"
