# ./histogram_ensemble -bin 0 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin0 -hist 100 439000 445000 -update
# ./histogram_ensemble -bin 0 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin0 -hist 100 439000 445000 -update
# ./histogram_ensemble -bin 0 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin0 -hist 100 439000 445000 -update

# ./histogram_ensemble -bin 1 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin1 -hist 100 356000 362000 -update
# ./histogram_ensemble -bin 1 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin1 -hist 100 356000 362000 -update
# ./histogram_ensemble -bin 1 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin1 -hist 100 356000 362000 -update

# ./histogram_ensemble -bin 2 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin2 -hist 100 136000 140000 -update
# ./histogram_ensemble -bin 2 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin2 -hist 100 136000 140000 -update
# ./histogram_ensemble -bin 2 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin2 -hist 100 136000 140000 -update

# ./histogram_ensemble -bin 3 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin3 -hist 100 42500 45500 -update
# ./histogram_ensemble -bin 3 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin3 -hist 100 42500 45500 -update
# ./histogram_ensemble -bin 3 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin3 -hist 100 42500 45500 -update

# ./histogram_ensemble -bin 4 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin4 -hist 100 11750 13500 -update
# ./histogram_ensemble -bin 4 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin4 -hist 100 11750 13500 -update
# ./histogram_ensemble -bin 4 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin4 -hist 100 11750 13500 -update

# ./histogram_ensemble -bin 5 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin5 -hist 100 2900 3900 -update
# ./histogram_ensemble -bin 5 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin5 -hist 100 2900 3900 -update
# ./histogram_ensemble -bin 5 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin5 -hist 100 2900 3900 -update

# ./histogram_ensemble -bin 6 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin6 -hist 100 500 1300 -update
# ./histogram_ensemble -bin 6 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin6 -hist 100 500 1300 -update
# ./histogram_ensemble -bin 6 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin6 -hist 100 500 1300 -update

# ./histogram_ensemble -bin 7 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin7 -hist 100 0 450 -update
# ./histogram_ensemble -bin 7 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin7 -hist 100 0 450 -update
# ./histogram_ensemble -bin 7 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin7 -hist 100 0 450 -update

# ./histogram_ensemble -bin 8 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin8 -hist 100 10 110 -update
# ./histogram_ensemble -bin 8 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin8 -hist 100 10 110 -update
# ./histogram_ensemble -bin 8 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin8 -hist 100 -100 200 -update

# ./histogram_ensemble -bin 9 -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin9 -hist 100 5 17 -update
# ./histogram_ensemble -bin 9 -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin9 -hist 100 5 17 -update
# ./histogram_ensemble -bin 9 -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin9 -hist 100 -75 125 -update

exe="./histogram_ensemble"

bin="0"
n="100"
xmin="1762000"
xmax="1774000"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="1"
n="100"
xmin="1430000"
xmax="1441000"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="2"
n="100"
xmin="549000"
xmax="556500"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="3"
n="100"
xmin="173000"
xmax="177000"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="4"
n="100"
xmin="49500"
xmax="51750"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="5"
n="100"
xmin="13250"
xmax="14600"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="6"
n="100"
xmin="3300"
xmax="4000"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="7"
n="100"
xmin="800"
xmax="1150"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="8"
n="100"
xmin="125"
xmax="325"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="9"
n="100"
xmin="0"
xmax="120"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

bin="10"
n="100"
xmin="-15"
xmax="50"
${exe} -bin ${bin} -i ../data/elisa_analysis_\*.root -o ensemble.root elisa_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/bayes_analysis_\*.root -o ensemble.root bayes_bin${bin} -hist ${n} ${xmin} ${xmax} -update
${exe} -bin ${bin} -i ../data/mle_analysis_\*.root -o ensemble.root mle_bin${bin} -hist ${n} ${xmin} ${xmax} -update

