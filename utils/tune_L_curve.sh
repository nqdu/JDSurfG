#/usr/bin/bash
# this file is for joint inversion
# if it is applied to surface wave tomography, remember to change some commands!
inputfile=JointSG.in 

# clear previous files
mkdir -p storage

# loop around each smoothing factor
for i in 10 15 20 25 30 35 50 70 85 100 150 180 200 350 500 600 800 1000 2000 5000;
do
    # replace smooth factor
    newline="$i"" 0.01                          # smooth damp"
    line=`grep "smooth damp" $inputfile`
    sed "s/$line/$newline/g" $inputfile > tmp.out 

    # set maxiteration = 1 
    # on L-curve analysis, you only need one iteration to see the relative
    # change of ||d-dobs|| and ||Lm||
    newline="1                                  # maximum iteration"
    line=`grep "maximum iteration" $inputfile`
    sed "s/$line/$newline/g" tmp.out >  storage/results${i}/$inputfile
    rm tmp.out
    
    # inversion 
    ../bin/JointTomo storage/results${i}/$inputfile surfdataSC.dat obsgrav.dat gravmat.dat MOD.surf MOD.ref >joint.out
    
    # backup results
    mkdir -p storage/results${i}
    mv results/* storage/results${i}/
    cp joint.out storage/results${i}/
done
