# make - Compiles the file euler.f90
# make clean - removes the data files
# make plot - plots the data files into output.ps file

all: heun_euler.f90
	gfortran heun_euler.f90 -o heun
clean:
	rm *.txt
plot: data_plot.gp
	gnuplot data_plot.gp
