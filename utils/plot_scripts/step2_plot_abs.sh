#!/bin/bash

. parameters.sh

mkdir -p pics

for idx in $run_idx;
do 
  MODEL=M`echo $idx | awk '{printf("%03d",$1)}'`
  
  # horizontal slice
  for ((i=0;i<$NSLICE_H;i++));
  do
    dep=${DEPTH_H[$i]}
    echo "Plotting horizontal slice $i for model $MODEL at depth = $dep km ..."
    # get grd info 
    info=`gmt grdinfo -C grd/horiz.$MODEL.${i}.grd`
    echo $info
    zmin=`echo $info | awk '{print $6}'`
    zmax=`echo $info | awk '{print $7}'`
    gmt makecpt -Cseis -T$zmin/$zmax/100+n -D -Z > out.cpt
    bounds=-R${LON0_H}/${LON1_H}/${LAT0_H}/${LAT1_H}
    proj=-JM12c

    gmt begin pics/${MODEL}_horiz_slice${i} jpg
    gmt basemap $bounds $proj -Bxaf -Byaf+l"Latitude (deg)" -BWSne+t"Horizontal Slice at Depth = ${dep} km"
    gmt grdimage grd/horiz.$MODEL.${i}.grd -Cout.cpt -E200

    # plot stations
    grep -v '^#' $DISP_FILE | awk '{print $2,$1}' | gmt plot -St0.15c -Gred

    gmt colorbar -Cout.cpt -Bx0.5f0.1+l"Vs,km/s" -By+lkm/s
    gmt end 
  done 

  # vertical slice
  for ((i=0;i<$NSLICE_V;i++));
  do
    echo "Plotting vertical slice $i for model $MODEL ..."

    # get grd info 
    info=`gmt grdinfo -C grd/vert.$MODEL.${i}.grd`
    echo $info
    xmin=`echo $info | awk '{print $2}'`
    xmax=`echo $info | awk '{print $3}'`
    zmin=`echo $info | awk '{print $6}'`
    zmax=`echo $info | awk '{print $7}'`
    gmt makecpt -Cseis -T$zmin/$zmax/100+n -D -Z > out.cpt
    bounds=-R${xmin}/${xmax}/${MIN_DEPTH}/${MAX_DEPTH}
    proj=-JX12c/-6c

    gmt begin pics/${MODEL}_vert_slice${i} jpg
    gmt basemap $bounds $proj -Bxaf -Byaf+l"Depth (km)" -BWSne+t"Vertical Slice $i"
    gmt grdimage grd/vert.$MODEL.${i}.grd -Cout.cpt -E200
    gmt colorbar -Cout.cpt -Bx0.5f0.1+l"Vs,km/s" -By+lkm/s
    gmt end
  done

  \rm out.cpt
done 