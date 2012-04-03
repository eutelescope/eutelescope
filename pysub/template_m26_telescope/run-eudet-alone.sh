#!/bin/sh
#
# run1 =
# run2 =
# 
# seq run1 run2 | awk '{print "./submit-converter.py --config=config/ibl/config-04-2fei4-0deg-geoid0.cfg -r --hot "$1" "$1 }' | sh -x
# seq run1 run2 | awk '{print "./submit-clusearch.py --config=config/ibl/config-04-2fei4-0deg-geoid0.cfg -r --hot "$1" "$1 }' | sh -x
# seq run1 run2 | awk '{print "./submit-hitmaker.py  --config=config/ibl/config-04-2fei4-0deg-geoid0.cfg -r -o "$1" -n "$1"  run0"$1"-clu-p.slcio "}' | sh -x
# seq run1 run2 | awk '{print "./submit-align.py     --config=config/ibl/config-04-2fei4-0deg-geoid0.cfg -r -o "$1" -iPreAlignedHit  "$1"-hit.slcio "}' | sh -x
# seq run1 run2 | awk '{print "./submit-fitter.py    --config=config/ibl/config-04-2fei4-0deg-geoid0.cfg -r -o "$1" -a "$1"-align-db.slcio "$1"-hit.slcio "}' | sh -x


run1=2016 
run2=2016  
config="config-desy2012-02-anemone.cfg"

for run in $(seq $run1 $run2)
do
  echo "./submit-converter.py --config=config/"$config"  --hot "$run" "$run
  echo "./submit-clusearch.py --config=config/"$config"  --hot "$run" "$run
  echo "./submit-hitmaker.py  --config=config/"$config"  -o "$run"  run0"$run"-clu-p.slcio "
  echo "./submit-align.py     --config=config/"$config"  -o "$run" -f \"0\" -e \"\" -iPreAlignedHit  "$run"-hit.slcio "
  echo "./submit-fitter.py    --config=config/"$config"   -o "$run" -a "$run"-align-db.slcio "$run"-hit.slcio "
done
