algo=${1}
seed=${RANDOM}
if [ "${algo}" == "elisa" ] ; then 
	# elisa 
	# example/unfold -algorithm 0 -statistical_analysis 100 example/elisa_analysis_test.root
	example/unfold -algorithm 0 -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_${seed}.root -covariance
	# for((i=0;i<8;++i)) ; do
		# seed=${RANDOM}
		# echo "example/unfold -algorithm 0 -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_${seed}.root"
	# done
elif [ "${algo}" == "mle" ] ; then 
	# maximum likelihood
	example/unfold -algorithm 3 -seed ${seed} -statistical_analysis 1000000 example/mle_analysis_${seed}.root
elif [ "${algo}" == "bayes" ] ; then 
	# bayesian iteration
	example/unfold -algorithm 2 -seed ${seed} -iterations 5 -statistical_analysis 1000000 example/bayes_analysis_${seed}.root
fi
