//
// Created by Alex on 15.03.2015.
//

#ifndef _ARK_CPP_BULK_H_
#define _ARK_CPP_BULK_H_

#include "Arr3d.h"

//	---------------------------------------------------------------------------------
//	l - geometry index: l=2 - cylindrical coordinates, l=1 - cartesian coordinates
//	n1 - number of grid nodes in x1 direction for one processor
//	n2 - number of grid nodes in X2 direction for one processor
//	n3 - number of grid nodes in X3 direction for one processor
//	nStop - total number of time steps
//	nPrint - print interval

int l, n1, n2, n3, nStop, nPrint, nStep, n1_g, n2_g, n3_g;

//	---------------------------------------------------------------------------------
//	x1_e - coordinates of the east border
//	x1_w - coordinates of the west border
//	x2_n - coordinates of the north border
//	x2_s - coordinates of the south border
//	x3_t - coordinates of the top border
//	x3_b - coordinates of the bottom border

double x1_e, x1_w, x2_n, x2_s, x3_t, x3_b;

//	---------------------------------------------------------------------------------
//	dx1 - grid step along x1 axis
//	dx2 - grid step along X2 axis
//	dx3 - grid step along X3 axis

double dx1, dx2, dx3;

//	---------------------------------------------------------------------------------
//	x1(:) - coordinates of the grid nodes along the axis x1
//	X2(:) - coordinates of the grid nodes along the axis X2
//	X3(:) - coordinates of the grid nodes along the axis X3

double *x1, *x2, *x3;

//	----------------------------------------------------------------------------------
//	roCon(:,:,:) - conservative density at the current time step
//	u1Con(:,:,:) - conservative velocity at the current time step along the axis x1
//	u2Con(:,:,:) - conservative velocity at the current time step along the axis X2
//	u3Con(:,:,:) - conservative velocity at the current time step along the axis X3
//	tCon(:,:,:) - conservative temperature at the current time step

Arr3d *roCon, *u1Con, *u2Con, *u3Con, *tCon;

//	----------------------------------------------------------------------------------
//	ronCon(:,:,:) - conservative density on the next time step
//	u1nCon(:,:,:) - conservative velocity on the next time step along the axis x1
//	u2nCon(:,:,:) - conservative velocity on the next time step along the axis X2
//	u3nCon(:,:,:) - conservative velocity on the next time step along the axis X3
//	tnCon(:,:,:) - conservative temperature at the next time step

Arr3d *ronCon, *u1nCon, *u2nCon, *u3nCon, *tnCon;

//	----------------------------------------------------------------------------------
//	p1(:,:,:) -  pressure on the faces perpendicular to the axis x1
//	ro1(:,:,:) - density on the faces perpendicular to the axis x1
//	u11(:,:,:) - velocity along the axis x1 on the faces perpendicular to the axis x1
//	u21(:,:,:) - velocity along the axis X2 on the faces perpendicular to the axis x1
//	u31(:,:,:) - velocity along the axis X3 on the faces perpendicular to the axis x1
//	t1(:,:,:) -  temperature on the faces perpendicular to the axis x1

Arr3d *p1, *ro1, *u11, *u21, *u31, *t1;

//	----------------------------------------------------------------------------------
//	p2(:,:,:) -  pressure on the faces perpendicular to the axis X2
//	ro2(:,:,:) - density on the faces perpendicular to the axis X2
//	u12(:,:,:) - velocity along the axis x1 on the faces perpendicular to the axis X2
//	u22(:,:,:) - velocity along the axis X2 on the faces perpendicular to the axis X2
//	u32(:,:,:) - velocity along the axis X3 on the faces perpendicular to the axis X2
//	t2(:,:,:) -  temperature on the faces perpendicular to the axis X2

Arr3d *p2, *ro2, *u12, *u22, *u32, *t2;

//	----------------------------------------------------------------------------------
//	p3(:,:,:) -  pressure on the faces perpendicular to the axis X3
//	ro3(:,:,:) - density on the faces perpendicular to the axis X3
//	u13(:,:,:) - velocity along the axis x1 on the faces perpendicular to the axis X3
//	u23(:,:,:) - velocity along the axis X2 on the faces perpendicular to the axis X3
//	u33(:,:,:) - velocity along the axis X3 on the faces perpendicular to the axis X3
//	t3(:,:,:) -  temperature on the faces perpendicular to the axis X3

Arr3d *p3, *ro3, *u13, *u23, *u33, *t3;

//	-----------------------------------------------------------------------------------
//	R(:)- one-dimensional buffer array for storing Riemann R invariant calculated values
//	Q(:)- one-dimensional buffer array for storing Riemann Q invariant calculated values

double *rBuf, *qBuf, *tfBuf, *tbBuf, *u2fBuf, *u2bBuf, *u3fBuf, *u3bBuf;

//	-----------------------------------------------------------------------------------
//  condition(:,:,:) - 1 if solid medium - 0 if liquid medium

Arr3d *condition;

//	-----------------------------------------------------------------------------------
// averageRot(:,:,:) - average value of rotor

Arr3d *averageU1;
Arr3d *averageU2;
Arr3d *averageU3;

#endif //_ARK_CPP_BULK_H_
