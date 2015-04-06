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

