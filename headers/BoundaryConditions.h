//
// Created by Alex on 20.03.2015.
//

#ifndef _ARK_CPP_BOUNDARYCONDITIONS_H_
#define _ARK_CPP_BOUNDARYCONDITIONS_H_


#include "Arr3d.h"

class BoundaryConditions {
private:
    Arr3d *_u1, *_u2, *_u3, *_ro, *_t;
    bool _x1Period, _x2Period, _x3Period;
    int _n1, _n2, _n3;
    double _t0;
public:
    BoundaryConditions(Arr3d *u1, Arr3d *u2, Arr3d *u3, Arr3d *ro, Arr3d *t, bool x1Period, bool x2Period, bool x3Period,
            int n1, int n2, int n3, double t0);
    void alongX1();
    void alongX2();
    void alongX3();
};


#endif //_ARK_CPP_BOUNDARYCONDITIONS_H_
