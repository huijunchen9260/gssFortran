all:
	gfortran -c goldenSectionSearch.f90
	gfortran -o main *.o main.f90
