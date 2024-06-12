#!/bin/sh

# clean prior results
rm -rf amp_sub*

# prepare directories
AMP_SUB_DIRS=$(ls configs/ | sed 's/config/amp_sub/g' | sed 's/.csv//g' )
mkdir $AMP_SUB_DIRS

# execute notebooks from templates and configs

for d in amp_sub*
do
#echo $d
cd $d

CONTRAST=$(echo $d | sed 's/amp_sub_//g')

echo "computing $CONTRAST in $d
"

cp ../configs/config_$CONTRAST\.csv ./config.csv
cp ../template/farrell_preprocess.ipynb $CONTRAST\_preprocess.ipynb
cp ../template/farrell_results.ipynb $CONTRAST\_results.ipynb

cat config.csv

ls

jupyter nbconvert --execute --to html $CONTRAST\_preprocess.ipynb
jupyter nbconvert --execute --to html $CONTRAST\_results.ipynb

rm config.csv

ls

cd ..

done

echo "--------

ALL FINISHED

---------"