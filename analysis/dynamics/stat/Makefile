#FC=gfortran
#FFLAGS=-fcheck=all -fbacktrace 
FC=ifx
FFLAGS=-traceback -check noarg_temp_created

all:
	${FC} ${FFLAGS} stat.f90 voxels.f90 main.f90 

clean:
	rm -f *.o *.mod a.out *.dat *.xyz *.png
