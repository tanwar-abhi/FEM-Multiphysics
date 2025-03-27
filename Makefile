# Makefile to compile the solver
#
#
# For compiling using eigen3 headers.
# g++ -I /path/to/eigen/ my_program.cpp -o my_program
#
#
# Another option on Linux/MacOs is to symlink or copy the Eigen folder into /usr/local/include/. So we can compile the program with::
# g++ my_program.cpp -o my_program 

# C++ compiler type
CXX = g++
# CXX = clang++-12

# Libraries path flag
# LIB = -I Libraries/eigen3.4.0 -I Libraries/cppJson/
LIB = -I Libraries/eigen3.3 -I Libraries/cppJson/

# Path to source code
PATHSC = src/

# C++ standard
CXXFLAGS := -std=c++14

# This is the debug flag; "-g" for debugging using GDB
# DFLAG := -g -O2 -Wall
DFLAG = -g


all:Solver

readMesh.o: $(PATHSC)readMesh.hpp $(PATHSC)readMesh.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)readMesh.cpp

Quadratures.o: $(PATHSC)Quadratures.hpp $(PATHSC)Quadratures.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Quadratures.cpp

BasisFunction.o: $(PATHSC)BasisFunction.hpp $(PATHSC)BasisFunction.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)BasisFunction.cpp

Material.o: $(PATHSC)Material.hpp $(PATHSC)Material.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Material.cpp

DirichletBC.o: $(PATHSC)DirichletBC.hpp $(PATHSC)DirichletBC.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)DirichletBC.cpp

Boundary_HT.o: $(PATHSC)Boundary_HT.hpp $(PATHSC)Boundary_HT.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Boundary_HT.cpp

PostProcess.o: $(PATHSC)PostProcess.hpp $(PATHSC)PostProcess.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)PostProcess.cpp

Matrix_Assemble.o: $(PATHSC)Matrix_Assemble.hpp $(PATHSC)Matrix_Assemble.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Matrix_Assemble.cpp


# For transient problem
# Time_Discretization.o: $(PATHSC)Time_Discretization.hpp $(PATHSC)Time_Discretization.cpp
#	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Time_Discretization.cpp

structureBoundary.o: $(PATHSC)structureBoundary.hpp $(PATHSC)structureBoundary.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)structureBoundary.cpp

neumanBC.o: $(PATHSC)neumanBC.hpp $(PATHSC)neumanBC.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)neumanBC.cpp

SolverInputs.o: $(PATHSC)SolverInputs.hpp $(PATHSC)SolverInputs.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)SolverInputs.cpp

PartInputs.o: $(PATHSC)PartInputs.hpp $(PATHSC)PartInputs.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)PartInputs.cpp

BoundaryInputs.o: $(PATHSC)BoundaryInputs.hpp $(PATHSC)BoundaryInputs.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)BoundaryInputs.cpp

Element.o: $(PATHSC)Element.hpp $(PATHSC)Element.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c $(PATHSC)Element.cpp

main_HT.o: main_HT.cpp
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) -c main_HT.cpp

# Solver: PartInputs.o SolverInputs.o BoundaryInputs.o BasisFunction.o Quadratures.o Element.o readMesh.o Material.o DirichletBC.o Boundary_HT.o neumanBC.o Time_Discretization.o PostProcess.o Matrix_Assemble.o main_HT.o
# 	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) *.o -o RunCode -ljsoncpp

Solver: PartInputs.o SolverInputs.o BoundaryInputs.o BasisFunction.o Quadratures.o Element.o readMesh.o Material.o DirichletBC.o Boundary_HT.o neumanBC.o structureBoundary.o PostProcess.o Matrix_Assemble.o main_HT.o
	$(CXX) $(CXXFLAGS) $(DFLAG) $(LIB) *.o -o RunCode -ljsoncpp

###############################################################################################################################
.PHONY : cleanobj
cleanobj:
	-@rm -rv *.o


.PHONY : clean
clean : cleanobj
	-@rm -rv RunCode
	-@echo "Removed C++ object files and executables"


# Help Target
.PHONY : help
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo " "
	@echo "... all = The default if no target is provided."
	@echo "... clean = Removes all object(.o) and executable files."
	@echo "... cleanobj = Removes all object{.o} files only."
	@echo " "
	@echo "Build Options:"
	@echo "-j[N], --jobs[=N]  = For parallel build's, j is for setting parallel build jobs, "
	@echo "			N is the number of parallel build jobs."
	@echo "			For example parallel build using 2 threads : make -j2"

