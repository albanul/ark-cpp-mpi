//
// Created by Alex on 15.03.2015.
//

#ifndef _ARK_CPP_ARR2D_H_
#define _ARK_CPP_ARR2D_H_


class Arr2d {
private:
    int _n1, _n2;
    double* _arr;
public:
    Arr2d(int n1, int n2);
    double& elem(int i, int j);
    double* ref();
    ~Arr2d();
};


#endif //_ARK_CPP_ARR2D_H_
