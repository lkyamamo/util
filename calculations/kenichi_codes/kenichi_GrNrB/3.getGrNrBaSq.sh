
# changeHToDerterium and getStatExcutable and getGrNrBaSq

cd stat
sh excutable.sh

echo "created executable"

module purge
module load legacy/CentOS7
module load intel-oneapi/2021.3

cd ../
stat/a.out pos_vel
cd ../..

