//
// Created by Alex on 20.03.2015.
//

#include "../headers/BoundaryConditions.h"

BoundaryConditions::BoundaryConditions(Arr3d *u1, Arr3d *u2, Arr3d *u3, Arr3d *ro, Arr3d *t, bool x1Period, bool x2Period, bool x3Period, int n1, int n2, int n3, double t0) {
    _u1 = u1;
    _u2 = u2;
    _u3 = u3;
    _ro = ro;
    _t = t;
    _x1Period = x1Period;
    _x2Period = x2Period;
    _x3Period = x3Period;
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _t0 = t0;
}

void BoundaryConditions::alongX1() {
    // along the X1 axis
    for (int j = 1; j < _n2; ++j) {
        for (int k = 1; k < _n3; ++k) {
            // on the west plane
            _u1->elem(_n1, j, k) = _u1->elem(_n1 - 1, j, k);
            _u2->elem(_n1, j, k) = _u2->elem(_n1 - 1, j, k);
            _u3->elem(_n1, j, k) = _u3->elem(_n1 - 1, j, k);
            _ro->elem(_n1, j, k) = _ro->elem(_n1 - 1, j, k);
            _t->elem(_n1, j, k) = _t->elem(_n1 - 1, j, k);

            // on the east plane
            _u1->elem(0, j, k) = _u1->elem(1, j, k);
            _u2->elem(0, j, k) = _u2->elem(1, j, k);
            _u3->elem(0, j, k) = _u3->elem(1, j, k);
            _ro->elem(0, j, k) = _ro->elem(1, j, k);
            _t->elem(0, j, k) = _t->elem(1, j, k);

        }
    }
}

void BoundaryConditions::alongX2() {
    // along the X2 axis
    for (int i = 1; i < _n1; i++) {
        for (int k = 1; k < _n3; k++) {
            // on the north plane
            _u1->elem(i, _n2, k) = _x2Period ? _u1->elem(i, 1, k) : 0.;
            _u2->elem(i, _n2, k) = _x2Period ? _u2->elem(i, 1, k) : 0.;
            _u3->elem(i, _n2, k) = _x2Period ? _u3->elem(i, 1, k) : 0.;
            _ro->elem(i, _n2, k) = _x2Period ? _ro->elem(i, 1, k) : _ro->elem(i, _n2 - 1, k);
            _t->elem(i, _n2, k) = _x2Period ? _t->elem(i, 1, k) : _t0;

            // on the south plane
            _u1->elem(i, 0, k) = _x2Period ? _u1->elem(i, _n2 - 1, k) : 0.;
            _u2->elem(i, 0, k) = _x2Period ? _u2->elem(i, _n2 - 1, k) : 0.;
            _u3->elem(i, 0, k) = _x2Period ? _u3->elem(i, _n2 - 1, k) : 0.;
            _ro->elem(i, 0, k) = _x2Period ? _ro->elem(i, _n2 - 1, k) : _ro->elem(i, 1, k);
            _t->elem(i, 0, k) = _x2Period ? _t->elem(i, _n2 - 1, k) : _t0;
        }
    }
}

void BoundaryConditions::alongX3() {
    for (int i = 1; i < _n1; ++i) {
        for (int j = 1; j < _n2; ++j) {
            // on the top plane
            _u1->elem(i, j, _n3) = _x3Period ? _u1->elem(i, j, 1) : 0.;
            _u2->elem(i, j, _n3) = _x3Period ? _u2->elem(i, j, 1) : 0.;
            _u3->elem(i, j, _n3) = _x3Period ? _u3->elem(i, j, 1) : 0.;
            _ro->elem(i, j, _n3) = _x3Period ? _ro->elem(i, j, 1) : _ro->elem(i, j, _n3 - 1);
            _t->elem(i, j, _n3) = _x3Period ? _ro->elem(i, j, 1) : _t0;

            // on the bottom plane
            _u1->elem(i, j, 0) = _x3Period ? _u1->elem(i, j, _n3 - 1) : 0.;
            _u2->elem(i, j, 0) = _x3Period ? _u2->elem(i, j, _n3 - 1) : 0.;
            _u3->elem(i, j, 0) = _x3Period ? _u3->elem(i, j, _n3 - 1) : 0.;
            _ro->elem(i, j, 0) = _x3Period ? _ro->elem(i, j, _n3 - 1) : _ro->elem(i, j, 1);
            _t->elem(i, j, 0) = _x3Period ? _t->elem(i, j, _n3 - 1) : _t0;
        }
    }
}


