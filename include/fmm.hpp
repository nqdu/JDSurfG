#include<iostream>
#include<vector>
#include<Eigen/Core>
#define earth 6371.0
#define pi 3.1415926535898
class backpointer
{
    public:
    int px,pz;
};

class FastMarching
{
    private:
    int nvx,nvz,nnx,nnz,fom,gdx,gdz;
    int vnl,vnr,vnt,vnb,nrnx,nrnz,sgdl,rbint;
    int nnxr,nnzr,asgr;
    Eigen::MatrixXi nsts,nstsr,srs;
    float gox,goz,dnx,dnz,dvx,dvz,snb;
    float goxd,gozd,dvxd,dvzd,dnxd,dnzd;
    float drnx,drnz,gorx,gorz,dnxr,dnzr,goxr,gozr;
    Eigen::MatrixXf velv,veln,velnb;
    Eigen::MatrixXf ttn,ttnr;
    int ntr;
    std::vector<backpointer> btg;

    int fouds1(int iz,int ix);
    int fouds2(int iz,int ix);
    int addtree(int iz,int ix);

    public:
    int travel();
};