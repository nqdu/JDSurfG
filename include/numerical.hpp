#ifndef _JDSURFG_NUMERICAL_H
#define _JDSURFG_NUMERICAL_H
#include <unsupported/Eigen/CXX11/Tensor>

typedef Eigen::Array<float,-1,-1> fmat2;
typedef Eigen::Array<double,-1,-1> dmat2;
typedef Eigen::Tensor<float,3> fmat3;
typedef Eigen::Tensor<double,3> dmat3;
typedef Eigen::Array<int,-1,-1> imat2;
typedef Eigen::Array<float,-1,1> fvec;
typedef Eigen::Array<double,-1,1> dvec;
typedef Eigen::Array<int,-1,1> ivec;

#endif