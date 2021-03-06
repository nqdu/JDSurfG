# first run direct surface wave tomography
# use ../bin/DSurfTomo -h for help
../bin/DSurfTomo DSurfTomo.in surfdataSC.dat MOD MOD.true

# generate gravity matrix
# use ../bin/mkmat -h for help
../bin/mkmat DSurfTomo.in obsgrav.dat MOD 

# then compute gravity for results of DSurfTomo
#../bin/syngrav 

# run joint inversion
# use ../bin/JointTomo -h for help
../bin/JointTomo JointSG.in surfdataSC.dat obsgrav.dat gravmat.dat MOD MOD MOD.true
