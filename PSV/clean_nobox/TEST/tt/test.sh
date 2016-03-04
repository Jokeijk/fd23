if [ $# -ne 5 ] ;then
  echo "$0 h dt vp vs alpha"
  exit -1
fi


h=$1
dt=$2
vp=$3
vs=$4
alpha=$5

cat > homo.mdl <<END
BACK $vp $vs 3.0
END
cat > station <<END
0 0
END

cat > par <<END
usebox=0
strike=0
azimuth=45
dip=90
rake=0
output=0
U_sta_file=station
W_sta_file=station
xs=512
zs=512
model=homo.mdl
usepunch=0
nx=1024
nz=1024

gpuid=1

h=$h
dt=$dt
nt=1500
ntsnap=200

#source
srctime=Gaussian
alpha=$alpha

srctype=1

#PML
pml_r=1E-11
pml_dt=0.005
pml_v=30.0
pml_fc=2
npml=32
END

nbpsv2d par=par
