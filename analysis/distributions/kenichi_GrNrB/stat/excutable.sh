module purge
make clean 
module load legacy/CentOS7
module load intel-oneapi/2021.3

echo making

make
