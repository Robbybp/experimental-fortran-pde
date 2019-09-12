# Code used in PSE seminar

Source code written for this seminar is in two files: pde.f90 and pde_direct.f90. They (attempt to) solve the PDE with iterative and direct solvers respectively. 

pde.f90 depends on mgmres.f90, an implementation of the GMRES algorithm by Lili Ju and John Burkardt. It can be obtained by:

```sh
wget https://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.f90
```

and compiled:

```sh
gfortran -o pde pde.f90 mgmres.f90
```

pde_direct.f90 depends on HSL's MA38 code and a BLAS library, e.g. https://github.com/xianyi/OpenBLAS. MA38 must be requested from HSL's website. This code can be compiled, for example, by:

```sh
gfortran -o pde_direct pde_direct.f90 libma38.a libopenblas.a
```

Command line arguments for both codes are ```which_problem``` (i.e. ```steady``` or ```dynamic```), ```nt```, and ```nx``` (the number of time and space discretization points). ```nt``` must be provided even if the problem is ```steady```, but it can be anything. 
For example:

```sh
./pde dynamic 100 100
```
