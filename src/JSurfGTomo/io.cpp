#include "JSurfGTomo/JSurfGTomo.hpp"

void JSurfGTomo:: 
write_syn(const fvec &dsyn,const std::string &swdfile,
         const std::string &gravfile) const
{
    int m1 = surf.obst.size(), m2 = obsg.size();
    fvec data1 = dsyn.segment(0,m1), data2 = dsyn.segment(m1,m2);
    
    // write dispersion data
    this->surf.write_syn(data1,swdfile);

    // write gravity data
    FILE *fp = fopen(gravfile.c_str(),"w");
    for(int i = 0; i < m2; i ++){
        fprintf(fp,"%f %f %g %g\n",lon_grav[i],lat_grav[i],obsg[i],data2[i]);
    }
    fclose(fp);

}