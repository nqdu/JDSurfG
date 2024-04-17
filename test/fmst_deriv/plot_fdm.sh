#!/bin/bash
info=`gmt gmtinfo -C fdm_phase.dat`
xmin=`echo $info |awk '{print $1}'`
xmax=`echo $info |awk '{print $2}'`
ymin=`echo $info |awk '{print $3}'`
ymax=`echo $info |awk '{print $4}'`
bounds=-R$xmin/$xmax/$ymin/$ymax 
proj=-JM12c

dx=`echo $xmin $xmax | awk '{print ($2-$1)/127.}'`
dy=`echo $ymin $ymax | awk '{print ($2-$1)/127.}'`
echo $bounds

# surface 
gmt surface fdm_phase.dat -Gfdm_phase.grd -I$dx/$dy $bounds 
gmt surface fdm_group.dat -Gfdm_group.grd -I$dx/$dy $bounds 

# colorbar
info=`gmt grdinfo fdm_group.grd -C`
vmin=`echo $info |awk '{print $6}'`
#vmax=`echo $info |awk '{print $7}'`
vmax=`echo "0. - $vmin" |bc -l`
gmt makecpt -Cvik -I -T$vmin/$vmax/100+n -D -Z > grad.cpt 

gmt begin frechet jpg 
gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"Phase Velocity Frechet"
gmt grdimage fdm_phase.grd -Cgrad.cpt -E200 $bounds $proj 
gmt plot sr.dat -W1p,black
head -1 sr.dat | gmt plot -Sa0.5c -Gyellow
sed -n '2,$p' sr.dat | gmt plot -St0.5c -Gpurple

gmt basemap $bounds $proj -Bxaf -Byaf -BwSen+t"Group Velocity Frechet" -X13c
gmt grdimage fdm_group.grd -Cgrad.cpt -E200 $bounds $proj 
gmt plot sr.dat -W1p,black
head -1 sr.dat | gmt plot -Sa0.5c -Gyellow
sed -n '2,$p' sr.dat | gmt plot -St0.5c -Gpurple

gmt colorbar -Cgrad.cpt -G$vmin/0 -Bxaf+l"frechet derivative" -X-6.5c
gmt end 

\rm *.grd *.cpt 