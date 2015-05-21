//
// Created by Alex on 15.03.2015.
//

#ifndef _ARK_CPP_STRESS_H_
#define _ARK_CPP_STRESS_H_

#include "../headers/Arr3d.h"

//	--------------------------------------------------------------------------------------
//	SIGMA11(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis x1
//	SIGMA21(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis x1
//	SIGMA31(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis x1
//	SIGMA12(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis X2
//	SIGMA22(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis X2
//	SIGMA32(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis X2
//	SIGMA13(:,:,:) - friction stress in the direction X1 to the faces perpendicular to the axis X3
//	SIGMA23(:,:,:) - friction stress in the direction X2 to the faces perpendicular to the axis X3
//	SIGMA33(:,:,:) - friction stress in the direction X3 to the faces perpendicular to the axis X3

Arr3d   *sigm11, *sigm21, *sigm31,
        *sigm12, *sigm22, *sigm32,
        *sigm13, *sigm23, *sigm33;

#endif //_ARK_CPP_STRESS_H_
