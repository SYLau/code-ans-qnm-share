# code-ans-qnm-share
Fortran code for computing quasi-normal modes of anisotropic neutron stars.
This code has only been tested in Windows 10 OS with the gFortran compiler. It contains Fortran 2008 features and is incompatible with older versions of Fortran.

1. Before compiling the main program, compile the odepack library by running "install.sh" at "lib\odepack". This creates the compiled library file "libodepack.a".
2. Compile the program with "projects\proj_aniso\Sid_ld\Makefile" or use the CodeBlocks project file "projects\proj_aniso_ld_single.cbp" as an alternative for Windows users
3. To run the compiled program, cd to "projects\proj_aniso\Sid_ld" and execute "bin\proj_aniso_ld_single\proj_aniso_ld_single.exe"
4. Modify "projects\proj_aniso\Sid_ld\inlist_ld_single.nml" to adjust the anisotropic neutron star parameters and search frequency for the modes

It is recommended to have gnuplot installed to plot the results generated. Two gnuplot scripts are provided in "projects\proj_aniso\Sid_ld\results_single\" to plot the ingoing wave amplitude at the far-field region at a given range of frequencies.
