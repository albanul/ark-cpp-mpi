//
// Created by Alex on 15.03.2015.
//

#ifndef _ARK_CPP_ARR3D_H_
#define _ARK_CPP_ARR3D_H_


class Arr3d {
private:
    int _n1, _n2, _n3;
    double* _arr;
public:
    Arr3d(int n1, int n2, int n3);
    double& elem(int i, int j, int k);
    double* ref();
    ~Arr3d();
};


#endif //_ARK_CPP_ARR3D_H_
