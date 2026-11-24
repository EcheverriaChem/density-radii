#!/bin/bash -l
#SBATCH -J bondfixB
#SBATCH -o logs/bondfixB.o%A_%a
#SBATCH -e logs/bondfixB.e%A_%a
#SBATCH -p agustina_thin
#SBATCH --account=noncovnano
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=24:00:00
#SBATCH --array=1-3%50

# One ORCA job per line of joblist.txt. Writes status.txt in each job folder.
echo "host: $(hostname)   task: $SLURM_ARRAY_TASK_ID"
cd "$SLURM_SUBMIT_DIR" || exit 4
JOBDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" joblist.txt)
echo "jobdir: $JOBDIR"
[ -d "$JOBDIR" ] || { echo "ERROR: cannot find $JOBDIR"; exit 5; }

module load Orca/6.1.1 || { echo "ERROR: module load failed" | tee "$JOBDIR/status.txt"; exit 2; }

module list 2>&1

ORCA_EXE=$(command -v orca)
[ -z "$ORCA_EXE" ] && { echo "ERROR: orca not in PATH after module load" | tee "$JOBDIR/status.txt"; exit 3; }
echo "orca: $ORCA_EXE"
ORCADIR=$(dirname "$ORCA_EXE")

cd "$JOBDIR"
echo "running on $(hostname) with $ORCA_EXE" > status.txt
"$ORCA_EXE" job.inp > job.out 2> job.err
rc=$?
echo "orca exit code: $rc" | tee -a status.txt
"$ORCADIR/orca_2mkl" job -molden > /dev/null 2>&1 || echo "WARNING: orca_2mkl failed" | tee -a status.txt
rm -f job.tmp* job*.tmp 2>/dev/null
exit $rc
