#!/bin/bash 

set -e 

. parameters.sh

for idx in $run_idx;
do 
  MODEL=M`echo $idx | awk '{printf("%03d",$1)}'`

  if [ $idx -eq 0 ]; then
    continue
  fi

  # horizontal slice
  for ((i=0;i<$NSLICE_H;i++));
  do
    dep=${DEPTH_H[$i]}
    echo "Plotting variation of horiz slice $i for model $MODEL at depth = $dep km ..."

    # grd difference 
    initmod=grd/horiz.M000.${i}.grd
    curmod=grd/horiz.$MODEL.${i}.grd
    gmt grdmath $curmod $initmod SUB $initmod DIV 100 MUL = grd/diff_horiz.$MODEL.${i}.grd
    #gmt grdmath grd/depth.$dep.true.grd $initgrd  SUB $initgrd DIV 100 MUL = grd/diff.$dep.true.grd

    # get info
    info=`gmt grdinfo -C grd/diff_horiz.$MODEL.${i}.grd`
    echo $info
    zmin=`echo $info | awk '{print $6}'`
    zmax=`echo $info | awk '{print $7}'`
    # make zmin/zmax symmetric
    absmax=`echo $zmin $zmax | awk '{if(-$1>$2) print -$1; else print $2}'`
    zmin=`echo $absmax | awk '{print -$1}'`
    zmax=`echo $absmax | awk '{print $1}'`
    gmt makecpt -Cvik -T$zmin/$zmax/100+n -D -Z -I > out.cpt
    bounds=-R${LON0_H}/${LON1_H}/${LAT0_H}/${LAT1_H}
    proj=-JM12c

    gmt begin pics/diff_${MODEL}_horiz_slice${i} jpg
    gmt basemap $bounds $proj -Bxaf -Byaf+l"Latitude (deg)" -BWSne+t"Horizontal Slice at Depth = ${dep} km"
    gmt grdimage grd/diff_horiz.$MODEL.${i}.grd -Cout.cpt -E200
    grep -v '^#' $DISP_FILE | awk '{print $2,$1}' | gmt plot -St0.15c -Gred
    gmt colorbar -Cout.cpt -Bxaf+l"@[\delta \ln(V_s), \% @["
    gmt end 
  done 

  # vertical slice
  for ((i=0;i<$NSLICE_V;i++));
  do
    echo "Plotting variation of vert slice $i for model $MODEL ..."
    # grd difference
    initmod=grd/vert.M000.${i}.grd
    curmod=grd/vert.$MODEL.${i}.grd
    gmt grdmath $curmod $initmod SUB $initmod DIV 100 MUL = grd/diff_vert.$MODEL.${i}.grd
    #gmt grdmath grd/depth.$dep.true.grd $initgrd  SUB $initgrd DIV 100 MUL = grd/diff.$dep.true.grd    
    # get info
    info=`gmt grdinfo -C grd/diff_vert.$MODEL.${i}.grd`
    echo $info
    xmin=`echo $info | awk '{print $2}'`
    xmax=`echo $info | awk '{print $3}'`
    zmin=`echo $info | awk '{print $6}'`
    zmax=`echo $info | awk '{print $7}'`
    # make zmin/zmax symmetric
    absmax=`echo $zmin $zmax | awk '{if(-$1>$2) print -$1; else print $2}'`
    zmin=`echo $absmax | awk '{print -$1}'`
    zmax=`echo $absmax | awk '{print $1}'`
    gmt makecpt -Cvik -T$zmin/$zmax/100+n -D -Z -I > out.cpt
    bounds=-R${xmin}/${xmax}/${MIN_DEPTH}/${MAX_DEPTH}
    proj=-JX12c/-6c

    gmt begin pics/diff_${MODEL}_vert_slice${i} jpg
    gmt basemap $bounds $proj -Bxaf+l"Distance along slice (km)" -Byaf+l"Depth (km)" -BWSne+t"Vertical Slice ${i}"
    gmt grdimage grd/diff_vert.$MODEL.${i}.grd -Cout.cpt -E200
    gmt colorbar -Cout.cpt -Bxaf+l"@[\delta \ln(V_s), \% @["
    gmt end 
  done

  \rm out.cpt
done 