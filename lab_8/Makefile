# make - Compiles the file euler.f90
# make clean - removes the data files
# make plot - plots the data files into output.ps file

all: euler.f90
	gfortran euler.f90 -o euler
clean:
	rm *.txt
plot: data_plot.gp
	gnuplot data_plot.gp
