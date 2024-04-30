#!/bin/bash
. parameters.sh 

info=`gmt gmtinfo -C $RESULT_DIR/res_grav$MODEL.dat`
xmin=`echo $info |awk '{print $1}'`
xmax=`echo $info |awk '{print $2}'`
ymin=`echo $info |awk '{print $3}'`
ymax=`echo $info |awk '{print $4}'`
region="$xmin/$xmax/$ymin/$ymax"
dx=`echo "$xmin $xmax" |awk '{print ($2-$1) / 127.}'`
dy=`echo "$ymin $ymax" |awk '{print ($2-$1) / 127.}'`
awk '{print $1,$2,$3}' $RESULT_DIR/res_grav$MODEL.dat |gmt surface -R$region -I$dx/$dy -Gobsg.grd -Vq
awk '{print $1,$2,$4}' $RESULT_DIR/res_grav$MODEL.dat |gmt surface -R$region -I$dx/$dy -Gsyng.grd -Vq

# colorbar
vmin=`gmt grdinfo -C obsg.grd |awk '{print $6}'`
vmax=`gmt grdinfo -C obsg.grd  |awk '{print $7}'`
echo $vmin $vmax
gmt makecpt -Cturbo -T$vmin/$vmax/100+n -D -Z -I > out.cpt
proj=-JM12c

gmt begin gravity jpg 
gmt basemap -R$region $proj -Bxaf -Byaf+l"Depth (km)" -BWSne+t"Inverted Gravity"
gmt grdimage  syng.grd -Cout.cpt -E200
#gmt grdcontour out.grd -R$region $proj -Gd6c -A0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6, -Gd6c
#gmt colorbar -DjMR+w3c/.25c+o-0.8c/0c+m -Cout.cpt -Bx0.5f0.1+l"Vs,km/s" -By+lkm/s

gmt basemap -R$region $proj -Bxaf -Byaf -BwSne+t"Observed Gravity" -X13c
gmt grdimage  obsg.grd -Cout.cpt -E200
#gmt grdcontour out.grd -R$region $proj -Gd6c -A0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6, -Gd6c
gmt colorbar -Bxaf -X-6.5c -Cout.cpt -Bx0.5f0.1+l"dG,mGal" -By+lkm/s
gmt end

\rm *.grd *.cpt 
