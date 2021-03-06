#!/bin/bash
#SBATCH --job-name="sleepyABM"
#SBATCH --time=1-0:15:00
#### Make array * tasks = 1000
### use array 0 to n - 1
#   SBATCH --array=0-0
#   SBATCH --ntasks=10
#SBATCH --array=0-0
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32768

module load java/1.8.0_25
module load R-gcc/3.2.3
export JAVA_HOME=/usr/local/java/1.8.0_25/jdk1.8.0_25/ 
export JAVA=/usr/local/java/1.8.0_25/jdk1.8.0_25/bin/java 
export JAVAC=/usr/local/java/1.8.0_25/jdk1.8.0_25/bin/javac 
export JAVAH=/usr/local/java/1.8.0_25/jdk1.8.0_25/bin/javah 
export JAR=/usr/local/java/1.8.0_25/jdk1.8.0_25/bin/jar 
export JAVA_LD_LIBRARY_PATH="/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/lib/amd64/server:/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/lib/amd64:/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/../lib/amd64:/usr/java/packages/lib/amd64:/usr/lib64:/lib64:/lib:/usr/lib" 
export JAVA_LIBS="-L/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/lib/amd64/server -L/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/lib/amd64 -L/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/../lib/amd64 -L/usr/java/packages/lib/amd64 -L/usr/lib64 -L/lib64 -L/lib -L/usr/lib -ljvm"
export JAVA_CPPFLAGS="-I/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/../include -I/usr/local/java/1.8.0_25/jdk1.8.0_25/jre/../include/linux"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$JAVA_LD_LIBRARY_PATH

cd /vlsci/VR0212/mrke/achaves/abm

for i in `seq 1 $SLURM_NTASKS`
do
    srun --nodes=1 --ntasks=1 --cpus-per-task=1 R --no-save --args $(($SLURM_ARRAY_TASK_ID*$SLURM_NTASKS+$i)) 0 0 < RNL_twostrategies.R &
done
# IMPORTANT must wait for all to finish, or all get killed
wait