#!/bin/bash

#SNAP=../snapshot_4/snap_004
#SNAP=/import/cage3/kohji/tng50/snapshot_4/snap_004
SNAP=/import/cage1/kohji/tng50/snapshot_6/snap_006

rm -rf `basename ${SNAP}`.log

#./HI_forest ${SNAP} 0.0 0.0 0.0 0.1 0.1 1.0 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 0.0 0.0 0.1 0.1 1.0 
mv ${SNAP}_los.dat `basename ${SNAP}`_1Z.dat
./read_spectrum `basename ${SNAP}`_1Z.dat `basename ${SNAP}`_1Z

#./HI_forest ${SNAP} 35000.0 0.0 0.0 -0.1 0.1 1.0 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 35000.0 0.0 0.0 -0.1 0.1 1.0 
mv ${SNAP}_los.dat `basename ${SNAP}`_2Z.dat
./read_spectrum `basename ${SNAP}`_2Z.dat `basename ${SNAP}`_2Z

#./HI_forest ${SNAP} 0.0 35000.0 0.0 0.1 -0.1 1.0 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 35000.0 0.0 0.1 -0.1 1.0 
mv ${SNAP}_los.dat `basename ${SNAP}`_3Z.dat
./read_spectrum `basename ${SNAP}`_3Z.dat `basename ${SNAP}`_3Z

#./HI_forest ${SNAP} 35000.0 35000.0 0.0 -0.1 -0.1 1.0 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 35000.0 35000.0 0.0 -0.1 -0.1 1.0 
mv ${SNAP}_los.dat `basename ${SNAP}`_4Z.dat
./read_spectrum `basename ${SNAP}`_4Z.dat `basename ${SNAP}`_4Z

# ------------------------------------------------------------

#./HI_forest ${SNAP} 0.0 0.0 0.0 0.1 1.0 0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 0.0 0.0 0.1 1.0 0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_1Y.dat
./read_spectrum `basename ${SNAP}`_1Y.dat `basename ${SNAP}`_1Y

#./HI_forest ${SNAP} 35000.0 0.0 0.0 -0.1 1.0 0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 35000.0 0.0 0.0 -0.1 1.0 0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_2Y.dat
./read_spectrum `basename ${SNAP}`_2Y.dat `basename ${SNAP}`_2Y

#./HI_forest ${SNAP} 0.0 0.0 35000.0 0.1 1.0 -0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 0.0 35000.0 0.1 1.0 -0.1
mv ${SNAP}_los.dat `basename ${SNAP}`_3Y.dat
./read_spectrum `basename ${SNAP}`_3Y.dat `basename ${SNAP}`_3Y

#./HI_forest ${SNAP} 35000.0 0.0 35000.0 -0.1 1.0 -0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 35000.0 0.0 35000.0 -0.1 1.0 -0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_4Y.dat
./read_spectrum `basename ${SNAP}`_4Y.dat `basename ${SNAP}`_4Y


# ------------------------------------------------------------

#./HI_forest ${SNAP} 0.0 0.0 0.0 1.0 0.1 0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 0.0 0.0 1.0 0.1 0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_1X.dat
./read_spectrum `basename ${SNAP}`_1X.dat `basename ${SNAP}`_1X

#./HI_forest ${SNAP} 0.0 35000.0 0.0 1.0 -0.1 0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 35000.0 0.0 1.0 -0.1 0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_2X.dat
./read_spectrum `basename ${SNAP}`_2X.dat `basename ${SNAP}`_2X

#./HI_forest ${SNAP} 0.0 0.0 35000.0 1.0 0.1 -0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 0.0 35000.0 1.0 0.1 -0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_3X.dat
./read_spectrum `basename ${SNAP}`_3X.dat `basename ${SNAP}`_3X

#./HI_forest ${SNAP} 0.0 35000.0 35000.0 1.0 -0.1 -0.1 &>> `basename ${SNAP}`.log
./HI_forest ${SNAP} 0.0 35000.0 35000.0 1.0 -0.1 -0.1 
mv ${SNAP}_los.dat `basename ${SNAP}`_4X.dat
./read_spectrum `basename ${SNAP}`_4X.dat `basename ${SNAP}`_4X
