
reference_panel=data/1000GenomesReferenceFiles/EUR

for exp in *_exposure.csv; do
  out=${exp/_exposure.csv/_outcome.csv}
  base=${exp/_exposure.csv/}
  ./pwcoco/build/pwcoco --bfile $reference_panel \
           --sum_stats1 $exp \
           --sum_stats2 $out \
           --coloc_pp 5e-5 5e-5 1e-6 \
           --out ${base}_coloc
done
