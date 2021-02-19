#/usr/bin/bash
mkdir -p storage
for i in 10 15 20 25 30 35 50 70 85 100 150 180 200 350 500 600 800 1000 2000 5000;do
    echo $i
    mkdir -p storage/results${i}
    sed "4c $i 0.0                        c: weight damp " JointSG.in > tmp.out
    sed "7c 1                                c: maximum iteration" tmp.out > storage/results${i}/JointSG.in
    rm tmp.out
    
    # inversion 
    ../bin/JointTomo storage/results${i}/JointSG.in surfdataSC.dat obsgrav.dat ../sichuan_test/gravmat.dat MOD.surf MOD.ref >joint.out
    
    # copy results
    mv  results/* storage/results${i}/
    cp results.bak/*.dat storage/results${i}/
    cp joint.out storage/results${i}/
done
