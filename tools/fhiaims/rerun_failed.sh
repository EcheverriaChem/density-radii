#!/bin/bash
# Rerun the FHI-aims EOS points that stopped without "Have a nice day".
# Run it from the eos_aims folder ONLY after the aims_eos array has finished
# (squeue -u $USER shows nothing), otherwise running points look failed.
#
#   bash rerun_failed.sh          # list failed points and submit them
#   bash rerun_failed.sh --dry    # only list them
#
# Each failed point gets the keywords that handle an ill-conditioned overlap
# matrix (near-linear dependence of the basis at short interatomic distances)
# and a longer SCF limit. Smearing, basis and k-grid are unchanged, so the
# EOS stays consistent. Each point runs as its own array task.

cd "$(dirname "$0")" 2>/dev/null
: > failed.txt
for d in */; do
  d=${d%/}
  [ -f "$d/control.in" ] || continue
  if [ ! -f "$d/aims.out" ] || ! grep -q "Have a nice day" "$d/aims.out"; then
    echo "$d" >> failed.txt
  fi
done
n=$(wc -l < failed.txt)
echo "$n unfinished points (listed in failed.txt)"
[ "$n" -eq 0 ] && exit 0
[ "$1" == "--dry" ] && { cat failed.txt; exit 0; }

# Extra settings, added once per point (skipped if already present).
while read d; do
  c="$d/control.in"
  grep -q "# rerun settings\|override_illconditioning" "$c" && { [ -f "$d/aims.out" ] && mv "$d/aims.out" "$d/aims.out.failed"; continue; }
  sed -i '/^sc_iter_limit/d' "$c"
  # insert before the species block (first line starting with "  species" or "species")
  awk 'BEGIN{done=0}
       /^[[:space:]]*species/ && !done {print "# rerun settings\nsc_iter_limit       400\nbasis_threshold     1.0E-5\noverride_illconditioning .true.\n"; done=1}
       {print}' "$c" > "$c.tmp" && mv "$c.tmp" "$c"
  [ -f "$d/aims.out" ] && mv "$d/aims.out" "$d/aims.out.failed"
done < failed.txt

mkdir -p logs
cat > rerun_eos.sh <<EOF
#!/bin/bash -l
#SBATCH -J aims_rerun
#SBATCH -o logs/rerun.o%A_%a
#SBATCH -e logs/rerun.e%A_%a
#SBATCH -p agustina_thin
#SBATCH --account=noncovnano
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --time=24:00:00
#SBATCH --array=1-$n%20
cd "\$SLURM_SUBMIT_DIR" || exit 4
module load gnu12 openmpi4/4.1.4 hwloc/2.7.0 fhi-aims/260331 || exit 2
ulimit -s unlimited
export OMP_NUM_THREADS=1
export TMPDIR=/tmp/aimsr_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}; export OMPI_MCA_orte_tmpdir_base=\$TMPDIR; mkdir -p \$TMPDIR
d=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" failed.txt)
( cd "\$d" && mpirun -np \$SLURM_NTASKS \$FHIAIMS_BINARY > aims.out 2> aims.err < /dev/null )
rm -rf \$TMPDIR
EOF
sbatch rerun_eos.sh
