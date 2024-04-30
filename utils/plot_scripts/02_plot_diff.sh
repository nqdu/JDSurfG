#!/bin/bash
. parameters.sh

# set your profile if required
gmt project -C105/35.3 -E105/25.7 -G1.0 -Q> prof.dat

# interpolate
depth=(`awk '{print $3}' $RESULT_DIR//mod_iter$MODEL.dat |sort -n |uniq`)
nz=`echo "${#depth[@]} - 1"|bc -l`
args=""
:> out.dat
:> out.true.dat
for ((i=0;i<$nz;i++));
do  
    dep=${depth[$i]}
    awk '{print $1,$2}' prof.dat | gmt grdtrack -Ggrd/diff.$dep.surf.grd -donodata > tmp.dat
    n=`wc -l tmp.dat |awk '{print $1}'`
    
    head -$n prof.dat > tmp1.dat
    #echo $n `wc -l tmp1.dat`
    paste tmp1.dat tmp.dat > tmp2.dat 
    awk -v a=$dep '{print $3,$2*0+a,$6}' tmp2.dat >> out.dat  

    awk '{print $1,$2}' prof.dat | gmt grdtrack -Ggrd/diff.$dep.true.grd -donodata > tmp.dat
    n=`wc -l tmp.dat |awk '{print $1}'`
    
    head -$n prof.dat > tmp1.dat
    #echo $n `wc -l tmp1.dat`
    paste tmp1.dat tmp.dat > tmp2.dat 
    awk -v a=$dep '{print $3,$2*0+a,$6}' tmp2.dat >> out.true.dat  
done 
#exit

info=`gmt gmtinfo -C out.dat`
xmin=`echo $info |awk '{print $1}'`
xmax=`echo $info |awk '{print $2}'`
ymin=`echo $info |awk '{print $3}'`
ymax=`echo $info |awk '{print $4}'`
region="$xmin/$xmax/$ymin/$ymax"
dx=`echo "$xmin $xmax" |awk '{print ($2-$1) / 127.}'`
dy=`echo "$ymin $ymax" |awk '{print ($2-$1) / 127.}'`

# to grd
gmt surface out.dat -Gout.grd -I$dx/$dy  -R$region -Vq
gmt surface out.true.dat -Gout.true.grd -I$dx/$dy  -R$region -Vq
vmin=`gmt grdinfo -C out.true.grd |awk '{print $6}'`
vmax=`gmt grdinfo -C out.true.grd  |awk '{print $7}'`
echo $vmin $vmax
gmt makecpt -Cturbo -T$vmin/$vmax/100+n -D -Z -I > out.cpt

proj=-JX12c/-6c

gmt begin profA jpg 
gmt basemap -R$region $proj -Bxaf -Byaf+l"Depth (km)" -BWSne+t"Inverted Model"
gmt grdimage  out.grd -Cout.cpt -E200
#gmt grdcontour out.grd -R$region $proj -Gd6c -A0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6, -Gd6c
#gmt colorbar -DjMR+w3c/.25c+o-0.8c/0c+m -Cout.cpt -Bx0.5f0.1+l"Vs,km/s" -By+lkm/s

gmt basemap -R$region $proj -Bxaf -Byaf -BwSne+t"True model" -X13c
gmt grdimage  out.true.grd -Cout.cpt -E200
#gmt grdcontour out.grd -R$region $proj -Gd6c -A0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6, -Gd6c
gmt colorbar -Bxaf -X-6.5c -Cout.cpt -Bx0.5f0.1+l"Vs,km/s" -By+lkm/s
gmt end

rm *.grd *.cpt tmp* out.dat out.true.dat
