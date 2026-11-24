#!/bin/bash -l
#SBATCH -J {jobname}
#SBATCH -o logs/{jobname}.o%A_%a
#SBATCH -e logs/{jobname}.e%A_%a
#SBATCH -p {partition}
#SBATCH --account={account}
#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={mem_per_cpu}M
#SBATCH --time={walltime}
#SBATCH --array=1-{njobs}%50

# One ORCA job per line of joblist.txt. Writes status.txt in each job folder.
echo "host: $(hostname)   task: $SLURM_ARRAY_TASK_ID"
cd "$SLURM_SUBMIT_DIR" || exit 4
JOBDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" joblist.txt)
echo "jobdir: $JOBDIR"
[ -d "$JOBDIR" ] || {{ echo "ERROR: cannot find $JOBDIR"; exit 5; }}

module load {orca_module} || {{ echo "ERROR: module load failed" | tee "$JOBDIR/status.txt"; exit 2; }}
{extra_module_line}
module list 2>&1

ORCA_EXE=$(command -v orca)
[ -z "$ORCA_EXE" ] && {{ echo "ERROR: orca not in PATH after module load" | tee "$JOBDIR/status.txt"; exit 3; }}
echo "orca: $ORCA_EXE"
ORCADIR=$(dirname "$ORCA_EXE")

# private OpenMPI/ORCA temporary folder: several ORCA jobs of the same user on one node
# otherwise share /tmp session files, and the first job to finish kills the others
export TMPDIR=/tmp/orca_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}
export OMPI_MCA_orte_tmpdir_base=$TMPDIR
mkdir -p "$TMPDIR"
cd "$JOBDIR"
echo "running on $(hostname) with $ORCA_EXE" > status.txt
"$ORCA_EXE" job.inp > job.out 2> job.err
rc=$?
echo "orca exit code: $rc" | tee -a status.txt
"$ORCADIR/orca_2mkl" job -molden > /dev/null 2>&1 || echo "WARNING: orca_2mkl failed" | tee -a status.txt
rm -f job.tmp* job*.tmp 2>/dev/null
rm -rf "$TMPDIR"
exit $rc
