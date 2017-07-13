PROJECT(adt_dist)
ENABLE_LANGUAGE(Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_C_FLAGS "-fPIC -O2")
#SET(CMAKE_Fortran_FLAGS "-fbounds-check -fbacktrace -fdefault-real-8 -fconvert=big-endian")
#SET(CMAKE_Fortran_FLAGS "-O2 -fdefault-real-8 -fconvert=big-endian")
SET(CMAKE_Fortran_FLAGS "-O2")
#SET(CMAKE_Fortran_FLAGS "-CB -traceback")

set(DIST_SRCS
buildADT.f90
driver.f90
PointTriangleDistance.f90
searchADT.f90
tribfdata.f90   
MathOps.f90
DistCalc.f90
cartwrite.f90
)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
ADD_EXECUTABLE(cutsurf.exe ${DIST_SRCS})

