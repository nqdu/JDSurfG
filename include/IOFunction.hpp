#pragma once
#include<iostream>
#include<fstream>
#include"tomography.hpp"

void inline skipread(std::ifstream &infile,std::string &line)
{
    // skip all comment lines
    line[0] = '#';
    while(line[0] == '#' || line.length() == 0){
        getline(infile,line);
    }
}

template<typename... Args>
void skipread(std::ifstream &infile,std::string &line,
                const char *__restrict__ __format, Args ...args)
{
    // skip all comment lines
    line[0] = '#';
    while(line[0] == '#' || line.length() == 0){
        getline(infile,line);
    }

    // unpack input variables
    sscanf(line.c_str(),__format,args...);
}

void __read_traveltime(std::string &surfdata,SurfTime &surf);
int __read_parfile(std::ifstream &infile,std::string &paramfile,MOD3d &mod,
                SurfTime &surf,InverseParamsBase &param);
void __read_InputModel(std::string &modfile,Eigen::Tensor<float,3> &vs,
                        float *dep,bool print_depth=false);