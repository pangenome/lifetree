#!/bin/bash
#SBATCH --job-name=pylauncher_example
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=2-00:00:00

# note: must specify --nodes=N on command line

cmdfile=$1

module load python
module load pylauncher

python -c "from pylauncher import ClassicLauncher; ClassicLauncher('"$cmdfile"')"
