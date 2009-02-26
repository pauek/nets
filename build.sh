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
   echo 'In Ubuntu/Debian, type "aptitude install blitz++"'
   exit 1
fi

echo -e "Compiling programs...\n"
make -C graphics
make -C networks

files=$(find ! -type d ! -name "*.sh" -perm '/u+x')
mkdir bin
for f in $files; do
  cp $f bin/.
done
echo
echo "+----------------------------------------------------------+"
echo "|                                                          |"
echo '|   You will find the executables in the "bin" directory   |'
echo "|                                                          |"
echo "+----------------------------------------------------------+"
echo
