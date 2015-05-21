//
// Created by Alex on 15.03.2015.
//

#ifndef _ARK_CPP_CONSTANTS_H_
#define _ARK_CPP_CONSTANTS_H_

//	-----------------------------------------------------------------------------------
//	sound	-	sound speed
//	ro0_g	-	unperturbed density of the liquid
//	dt		-	time step
//	VIS		-	kinematic viscosity
//	u10		-	initial speed along the axis x1
//	u20		-	initial speed along the axis X2
//	u30		-	initial speed along the axis X3
//	t0		-	initial temperature
//	TIME	-	current time
//	CFL		-	Courant number

double sound, ro0_g, dt, VIS, u10, u20, u30, t0, TIME, CFL;
bool x1Period, x2Period, x3Period, needSwap, useTecplot;
const char* dirPath;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#endif //_ARK_CPP_CONSTANTS_H_
