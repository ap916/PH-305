all: gaussElimination.f90 gaussPivot.f90
	gfortran gaussElimination.f90 -o elim.out
	gfortran gaussPivot.f90 -o pivot.out
	g++ matmul-assign.cpp -o matmul.out
clean:
	rm *.out
