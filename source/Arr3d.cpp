//
// Created by Alex on 15.03.2015.
//

#include "../headers/Arr3d.h"

Arr3d::Arr3d(int n1, int n2, int n3) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;

    _arr = new double[_n1*_n2*_n3];
}

Arr3d::~Arr3d() {
    delete[] _arr;
}

double& Arr3d::elem(int i, int j, int k) {
    return _arr[k *_n1*_n2 + j*_n1 + i];
}

double* Arr3d::ref() {
	return _arr;
}