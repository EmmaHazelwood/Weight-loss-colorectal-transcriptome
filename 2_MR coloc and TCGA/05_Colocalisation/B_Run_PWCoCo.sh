
reference_panel=data/1000GenomesReferenceFiles/EUR

# PWCoCo binary; override with e.g. PWCOCO_BIN=~/src/pwcoco/build/pwcoco bash B_Run_PWCoCo.sh
pwcoco_bin=${PWCOCO_BIN:-./pwcoco/build/pwcoco}

for exp in *_exposure.csv; do
  out=${exp/_exposure.csv/_outcome.csv}
  base=${exp/_exposure.csv/}
  # Skip trios already run; delete ${base}_coloc.coloc to force a re-run
  [ -f "${base}_coloc.coloc" ] && continue
  "$pwcoco_bin" --bfile $reference_panel \
           --sum_stats1 $exp \
           --sum_stats2 $out \
           --coloc_pp 5e-5 5e-5 1e-6 \
           --out ${base}_coloc
done
