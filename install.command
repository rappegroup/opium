#!/bin/sh
cd "$(dirname "$0")"
./configure
make
#sudo cp opium /usr/local/bin/opium
echo "*****************************************************************************"
echo "Installation finished! Run opium using command line in a new terminal window!"
echo "The usual location for the executable is /usr/local/bin"
echo "*****************************************************************************"
#make distcheck
# make install
