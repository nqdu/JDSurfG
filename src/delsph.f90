subroutine delsph(flat1,flon1,flat2,flon2,del)bind(c,name="delsph")
    use,intrinsic :: iso_c_binding
    implicit none
    real(c_float),parameter:: R=6371.0
    REAL(c_float),parameter:: pi=3.1415926535898
    real(c_float),value,INTENT(In) :: flat1,flat2,flon1,flon2
    real(c_float),INTENT(OUT) :: del
  
    real(c_float) dlat,dlon,lat1,lat2,a,c
  
    dlat=flat2-flat1
    dlon=flon2-flon1
    lat1=pi/2-flat1
    lat2=pi/2-flat2
    a=sin(dlat/2)*sin(dlat/2)+sin(dlon/2)*sin(dlon/2)*cos(lat1)*cos(lat2)
    c=2*atan2(sqrt(a),sqrt(1-a))
    del=R*c
  end subroutine
  