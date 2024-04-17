#!/bin/bash
info=`gmt gmtinfo -C grad.org.dat`
xmin=`echo $info |awk '{print $1}'`
xmax=`echo $info |awk '{print $2}'`
ymin=`echo $info |awk '{print $3}'`
ymax=`echo $info |awk '{print $4}'`
bounds=-R$xmin/$xmax/$ymin/$ymax 
proj=-JX12c/5c

dx=`echo $xmin $xmax | awk '{print ($2-$1)/127.}'`
dy=`echo $ymin $ymax | awk '{print ($2-$1)/127.}'`
echo $bounds

# surface 
gmt surface grad.org.dat -Ggrad.org.grd -I$dx/$dy $bounds 
gmt surface grad.smooth.conv.dat  -Ggrad.smooth.conv.grd -I$dx/$dy $bounds
gmt surface grad.smooth.pde.dat  -Ggrad.smooth.pde.grd -I$dx/$dy $bounds

# compute diff
gmt grdmath grad.smooth.conv.grd grad.smooth.pde.grd SUB = diff.grd

# colorbar
gmt makecpt -Cpolar -I -T-1/1/100+n -D -Z > grad.cpt 

gmt begin smooth jpg 

gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"org" -X13c -Y13c
gmt grdimage grad.org.grd -Cgrad.cpt -E200

gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"smooth_conv" -X13c
gmt grdimage grad.smooth.conv.grd -Cgrad.cpt -E200

gmt basemap $bounds $proj -Bxaf -Byaf -BwSen+t"smooth_pde" -X-13c -Y-8c
gmt grdimage grad.smooth.pde.grd -Cgrad.cpt -E200

gmt basemap $bounds $proj -Bxaf -Byaf -BwSen+t"smooth_diff" -X13c
gmt grdimage diff.grd -Cgrad.cpt -E200

gmt colorbar -Cgrad.cpt -Bxaf -X-6.5c 
gmt end 

\rm *.grd *.cpt 