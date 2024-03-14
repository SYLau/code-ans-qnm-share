# code-ans-qnm-share
Fortran code for computing quasi-normal modes of anisotropic neutron stars
This code has only been tested in Windows 10 OS with the gFortran compiler. It contains Fortran 2008 features and is incompatible with older versions of Fortran.

1. Before compiling the main program, compile the odepack library by running lib\odepack\install.sh.
2. Compile the program with projects\proj_aniso\Sid_ld\Makefile
3. cd to projects\proj_aniso\Sid_ld and run bin\proj_aniso_ld_single\proj_aniso_ld_single.exe
4. Modify projects\proj_aniso\Sid_ld\inlist_ld_single.nml to adjust the anisotropic neutron star parameters and search frequency for the modes
