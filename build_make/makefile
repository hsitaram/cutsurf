SOURCEDIR=../src
ofiles = ${SOURCEDIR}/MathOps.o ${SOURCEDIR}/searchADT.o ${SOURCEDIR}/PointTriangleDistance.o \
	 ${SOURCEDIR}/tribfdata.o ${SOURCEDIR}/buildADT.o  ${SOURCEDIR}/DistCalc.o ${SOURCEDIR}/cartwrite.o 

FLAGS=-O3

%.o: %.f90
	        gfortran -c ${FLAGS} $< -o $@


cutsurf.exe:${ofiles}
	gfortran ${ofiles} ${SOURCEDIR}/driver.f90 -o cutsurf.exe
clean:
	rm -rf ${SOURCEDIR}/*.o ${SOURCEDIR}/*.mod *.mod *.o *.exe
