# Windows bat
#
gfortran -c -Wall -std=legacy odepack.f
ar rcs libodepack.a odepack.o

del "odepack.o"

echo "Normal end of execution."

pause