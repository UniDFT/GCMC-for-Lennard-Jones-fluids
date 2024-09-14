#!/bin/csh
cd Run${MPIRUN_RANK}
#Replace following line with desired sequential program
../run
exit
