#### This script was launched with grid engine qsub. to execute a fit with MPI
#### use the bash shell
#!/bin/bash
#$ -S /bin/bash
#### replace $job_name with a meaningful name.
#$ -N fitneuron.bin
#### IMPORTANT: $num_CPUs must match with that below
#$ -pe mpich-f 72
#### Email notification at the beginning (b), the end (e),
#### and the suspension (s) of the job.
#$ -M brookings@brandeis.edu
#$ -m bes
## Call mpirun to run the parallel job. IMPORTRANT: the $num_CPUs variable 
## must match with the one in the -pe mpich $num_CPUs option specified above.
## Replace $job_command with your own. To run CHARMM, replace $job_command 
## with 'charmm < input > output'. To run NAMD, replace $job_command with
## 'namd2 control > output'. Be sure that either charmm or namd2 is in the
## path to the executables.
/home/tbrookin/OpenMPI/bin/mpirun -mca btl self,tcp -np 72 /home/tbrookin/cpp/fitneuron/bin/fitneuron.bin -startup startFit.txt -accuracy 1.0 -populationSize 943
