# Path to the shellscript interpreter
# /bin/bash
#! C:/msys64/usr/bin/bash
#
gfortran -c -Wall -std=legacy odepack.f
ar rcs libodepack.a odepack.o

rm -f odepack.o

echo "Normal end of execution."
