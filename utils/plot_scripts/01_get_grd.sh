#!/bin/bash
####################################
# This script is used to get different depth Vs.
#################################### 
. parameters.sh
outmodel=$RESULT_DIR/mod_iter$MODEL.dat
do_syn_test=true

#
rm -rf slice
mkdir -p slice

depth=(`awk '{print $3}' $RESULT_DIR/mod_iter$MODEL.dat |sort -n |uniq`)
nz=`echo "${#depth[@]} - 1"|bc -l`
for ((i=0;i<$nz;i++));
do
    dep=${depth[$i]}
    echo $dep
    awk -v a=$dep '($3 - a <= 1.0e-4 && $3-a >=-1.0e-4) {print $1,$2,$4}' $outmodel >slice/depth.${dep}.surf.dat
    if [ $do_syn_test ];then 
        awk -v a=$dep '($3 - a <= 1.0e-4 && $3-a >=-1.0e-4) {print $1,$2,$4}' $RESULT_DIR/mod_true.dat >slice/depth.${dep}.true.dat
        awk -v a=$dep '($3 - a <= 1.0e-4 && $3-a >=-1.0e-4) {print $1,$2,$4}' $RESULT_DIR/mod_iter0.dat >slice/depth.${dep}.init.dat
    fi
done 

rm -rf grd
mkdir -p grd

filename=`ls slice |head -1`
info=`gmt gmtinfo -C slice/$filename`
xmin=`echo $info |awk '{print $1}'`
xmax=`echo $info |awk '{print $2}'`
ymin=`echo $info |awk '{print $3}'`
ymax=`echo $info |awk '{print $4}'`
region="$xmin/$xmax/$ymin/$ymax"
dx=`echo "$xmin $xmax" |awk '{print ($2-$1) / 127.}'`
dy=`echo "$ymin $ymax" |awk '{print ($2-$1) / 127.}'`

for ((i=0;i<$nz;i++));
do
    dep=${depth[$i]}
    echo $dep
    gmt surface slice/depth.$dep.surf.dat -Ggrd/depth.$dep.surf.grd  -I$dx/$dy -R$region
    if [ $do_syn_test ] ;then 
        gmt surface slice/depth.$dep.true.dat -Ggrd/depth.$dep.true.grd  -I$dx/$dy -R$region -Vq
        gmt surface slice/depth.$dep.init.dat -Ggrd/depth.$dep.init.grd  -I$dx/$dy -R$region -Vq

        # relative variation
        initgrd=grd/depth.$dep.init.grd
        gmt grdmath grd/depth.$dep.true.grd $initgrd  SUB $initgrd DIV 100 MUL = grd/diff.$dep.true.grd
        gmt grdmath grd/depth.$dep.surf.grd $initgrd  SUB $initgrd DIV 100 MUL = grd/diff.$dep.surf.grd

        # colorbar 
        vmin=`gmt grdinfo -C grd/diff.$dep.surf.grd|awk '{print $6}'`
        vmax=`gmt grdinfo -C grd/diff.$dep.surf.grd|awk '{print $7}'`
        gmt makecpt -Cpolar -T-10/10/100+n -D -Z -I > out.cpt


        # plot checkerbaord
        gmt begin depth_$dep jpg 
            gmt basemap -R$region -JM12c -Bxaf -Byaf 
            gmt grdimage grd/diff.$dep.surf.grd -E200 -Cout.cpt 

            # stations
            cat $RESULT_DIR/res_swd0.dat |awk '{print $5,$4}' |gmt plot -St0.2c -Gblack
            cat $RESULT_DIR/res_swd0.dat |awk '{print $3,$2}' |gmt plot -St0.2c -Gblack
            

            gmt basemap -R$region -JM12c -Bxaf -Byaf  -X14c
            gmt grdimage grd/diff.$dep.true.grd -E200 -Cout.cpt 
            cat $RESULT_DIR/res_swd0.dat |awk '{print $5,$4}' |gmt plot -St0.2c -Gblack
            cat $RESULT_DIR/res_swd0.dat |awk '{print $3,$2}' |gmt plot -St0.2c -Gblack

            gmt colorbar -Bxaf -Cout.cpt -X-7c
        gmt end 
    fi 
done