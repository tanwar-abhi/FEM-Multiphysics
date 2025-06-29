# Multiphysics Product Development
Simulation of solid and structure mechanics problems.
The code uses open source [**Eigen3**](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for matrix multiplications, linear algebra and solving linear system of equations.
[Eigen-Documentation](https://eigen.tuxfamily.org/dox/group__QuickRefPage.html).

## Version 0.1.3
* Switched solvers from iterative to direct solver (LU based) to maintain stability in code.


## Compile and Run Code
To compile the code, invoke the call to Makefile i.e. by simply typing “make -jN” in command line.
* The `-j` is the jobs tag i.e. give directive to compiler for parallel build.
* `N` is the number of parallel jobs to run, this depends on the total threads and cores of user’s CPU.
* To Remove all object and executable files from main directory : *make clean*
* For more inforamation about make options, use syntax : *make help*

To run the code, go to directory contining mesh and JSON files then run excutable : **PATH-TO-MPS-Product-Dev/RunCode**


### Comments from Development Team

* The code needs a total of 3 JSON input files i.e. `solver` , `boundary` and `parts`. 
  * `solver.json` contains the inputs for solver, equations and material properties of the domain.
  * `boundary.json` file contains the boundary conditions {Dirichlet, Neumann, Initial} defined as per the problem by user.
  * `parts.json` contain all the subsequent parts of the whole domain {Multibody problems}. `parts` input must be defined for each problem even for non multibody problems.

* Read the [*Documentation*](https://github.com/bosonqpsi/MPS-Product-Dev/tree/master/documentation) for inputs JSON file format.
