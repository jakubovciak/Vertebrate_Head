#!/bin/sh

cd amp_G4toN2_zebra_5to9
jupyter nbconvert --execute --to html amp_G4_N2_zf_f_5_9_preprocess.ipynb
jupyter nbconvert --execute --to html amp_G4_N2_zf_f_5_9_results.ipynb
cd ..

cd amp_G4toN5_zebra_5to12
jupyter nbconvert --execute --to html amp_all_zf_f_50to12_preprocess.ipynb
jupyter nbconvert --execute --to html amp_all_zf_f_50to12_results.ipynb
cd ..

cd amp_N5_zebra_10to12
jupyter nbconvert --execute --to html amp_N5_zf_10_12_preprocess.ipynb
jupyter nbconvert --execute --to html amp_N5_zf_10_12_results.ipynb
cd ..


echo "--------

ALL FINISHED

---------"