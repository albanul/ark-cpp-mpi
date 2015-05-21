//
// Created by Alex on 01.05.2015.
//

#ifndef _ARK_CPP_PARALLEL_H_
#define _ARK_CPP_PARALLEL_H_

#include <mpi.h>

//	-----------------------------------------------------------------------------------
//  left        - previous processor along the X1 axis
//  right       - next processor along the X1 axis
//  down        - previous processor along the X2 axis
//  up          - next processor along the X2 axis
//  bottom      - previous processor along the X3 axis
//  top         - next processor along the X3 axis
//  nproc       - number of processors
//  rank        - processor number
//  newComm     - new communicator
//  coordX      - coordinates in topology
//  coordY      - coordinates in topology
//  coordZ      - coordinates in topology
//  px          - number of processors along the X1 axis
//  py          - number of processors along the X2 axis
//  pz          - number of processors along the X3 axis
//  filetype    - new MPItype

int left, right, up, down, top, bottom, nproc, rank, px, py, pz, coordX, coordY, coordZ;
int coords[3], dims[3];

MPI_Comm newComm;
MPI_Datatype filetype;

#endif //_ARK_CPP_PARALLEL_H_