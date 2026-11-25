#!/bin/bash
# Add the ill-conditioning keywords to every EOS point that has not finished yet.
# Safe while the array is running: FHI-aims reads control.in only at start-up,
# so running points are unaffected and pending points pick up the change.
# Run from the eos_aims folder.
n=0
for d in */; do
  d=${d%/}; c="$d/control.in"
  [ -f "$c" ] || continue
  [ -f "$d/aims.out" ] && grep -q "Have a nice day" "$d/aims.out" && continue
  grep -q "override_illconditioning" "$c" && continue
  awk 'BEGIN{done=0}
       /^[[:space:]]*species/ && !done {print "basis_threshold     1.0E-5\noverride_illconditioning .true.\n"; done=1}
       {print}' "$c" > "$c.tmp" && mv "$c.tmp" "$c"
  n=$((n+1))
done
echo "keywords added to $n unfinished points"
