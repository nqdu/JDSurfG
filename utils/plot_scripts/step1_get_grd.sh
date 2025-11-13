#!/bin/bash
####################################
# This script is used to get different depth Vs.
#################################### 
set -e
. parameters.sh

# generate slice
rm -rf slice grd 
mkdir -p slice grd 
for idx in $run_idx;
do
  MODEL=M`echo $idx | awk '{printf("%03d",$1)}'`
  outmodel=$RESULT_DIR/mod_iter${idx}.dat

  for ((i=0;i<$NSLICE_H;i++));
  do
    dep=${DEPTH_H[$i]}
    echo "Generating horizontal slice $i for model $MODEL at depth = $dep km ..."
    python ${PY_SCRIPTS}/generate_plane.py $LON0_H $LON1_H $LAT0_H $LAT1_H $dep $NGRD profile.txt

    # interpolate model values onto the slice
    python ${PY_SCRIPTS}/interp3d.py profile.txt $outmodel slice/horiz.$MODEL.$i.txt
    \rm profile.txt

    # make it to grd
    bounds=-R${LON0_H}/${LON1_H}/${LAT0_H}/${LAT1_H}
    dx=`echo "$LON1_H $LON0_H" | awk -v n=$NGRD '{print ($1-$2)/(n-1)}'`
    dy=`echo "$LAT1_H $LAT0_H" | awk -v n=$NGRD '{print ($1-$2)/(n-1)}'`
    inc=-I${dx}/${dy}
    awk '{print $1,$2,$4}' slice/horiz.$MODEL.$i.txt | gmt surface  $bounds $inc -Ggrd/horiz.$MODEL.${i}.grd -Vq
  done 

  # vertical slices
  for ((i=0;i<$NSLICE_V;i++));
  do
    echo "Generating vertical slice $i for model $MODEL ..."
    python ${PY_SCRIPTS}/generate_gc.py ${LON0_V[$i]} ${LAT0_V[$i]} ${LON1_V[$i]} ${LAT1_V[$i]} $MIN_DEPTH $MAX_DEPTH $NGRD profile.txt
    # interpolate model values onto the slice
    python ${PY_SCRIPTS}/interp3d.py profile.txt $outmodel slice/vert.$MODEL.$i.txt
    \rm profile.txt

    # make it to grd 
    info=`gmt gmtinfo -C slice/vert.$MODEL.$i.txt`
    dmin=`echo $info | awk '{print $7}'`
    dmax=`echo $info | awk '{print $8}'`
    bounds=-R${dmin}/${dmax}/${MIN_DEPTH}/${MAX_DEPTH}
    dx=`echo "$dmax $dmin" | awk -v n=$NGRD '{print ($1-$2)/(n-1)}'`
    dy=`echo "$MAX_DEPTH $MIN_DEPTH" | awk -v n=$NGRD '{print ($1-$2)/(n-1)}'`
    inc=-I${dx}/${dy}
    awk '{print $4,$3,$5}' slice/vert.$MODEL.$i.txt | gmt surface $bounds $inc -Ggrd/vert.$MODEL.${i}.grd -Vq
  done
done

if [ $HAS_TOPO -eq 1 ]; then
  echo "Generating topography/gravity/bathymetry grids ..."

  path=$RESULT_DIR/../topography.dat
  info=`head -3 $path`
  nlat=`echo $info |head -1| awk '{print $1}'`
  nlon=`echo $info |head -1| awk '{print $2}'`
  lat0=`echo $info |head -2 |tail -1 | awk '{print $1}'`
  lon0=`echo $info |head -2 |tail -1 | awk '{print $2}'`
  dlat=`echo $info | tail -1 | awk '{print $1}'`
  dlon=`echo $info | tail -1 | awk '{print $2}'`

  bounds=-R${lon0}/`echo "$lon0 $dlon $nlon" | awk '{printf("%.4f",$1+$2*($3-1))}'`/${lat0}/`echo "$lat0 $dlat $nlat" | awk '{printf("%.4f",$1+$2*($3-1))}'`
  sed -n '4,$p' $path | gmt xyz2grd  -Ggrd/topo.grd $bounds -I${dlon}/${dlat} -ZTL

  # generate slice
  for ((i=0;i<$NSLICE_V;i++));
  do 
    echo "Generating topography slice $i ..."
    MODEL=M`echo $run_idx | awk '{printf "%03d",$1}'`
    awk '{print $1,$2}' slice/vert.$MODEL.$i.txt | gmt grdtrack -Ggrd/topo.grd  > slice/vert.topo.$i.txt
  done 
fi