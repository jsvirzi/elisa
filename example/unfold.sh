# example/unfold -algorithm 0 -error_analysis 1000 example/elisa_error_analysis.root -prior example/ensemble_2.root elisa_bin

algo=${1}
n=${2}
if [ "${algo}" == "elisa" ] ; then 
	# elisa 
	# example/unfold -algorithm 0 -statistical_analysis 100 example/elisa_analysis_test.root
	# dfile="example/measurement_1.root"
	# dname="meas_x"
	dfile="example/measurement_1_meas.dat"
	tfile="example/measurement_1_true.dat"
	nstat="1"
	if [ "${n}" == "" ] ; then
		seed=${RANDOM}
		# example/unfold -algorithm 0 -seed ${seed} -meas ${dfile} ${dname} -statistical_analysis ${nstat} example/elisa_analysis_${seed}.root -covariance
		# example/unfold -algorithm 0 -meas ${dfile} ${dname} -statistical_analysis ${nstat} example/elisa_analysis.root
		example/unfold -algorithm 0 -format 0 -meas ${dfile} -statistical_analysis ${nstat} example/elisa_analysis_${seed}.root
	else
		for((i=0;i<${n};++i)) ; do
			seed=${RANDOM}
			example/unfold -algorithm 0 -seed ${seed} -meas ${dfile} ${dname} -statistical_analysis ${nstat} example/elisa_analysis_${seed}.root -covariance
		done
	fi
elif [ "${algo}" == "elisa1" ] ; then 
	# elisa 
	# example/unfold -algorithm 0 -statistical_analysis 100 example/elisa_analysis_test.root
	dfile="example/measurement_1.root"
	dname="meas_x"
	if [ "${n}" == "" ] ; then
		seed=${RANDOM}
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_1_${seed}.root -covariance
	else
		for((i=0;i<${n};++i)) ; do
			seed=${RANDOM}
			example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_1_${seed}.root -covariance
		done
	fi
elif [ "${algo}" == "combined" ] ; then 
	# elisa 
	dfile="example/measurement_1_2.root"
	dname="meas_x"
	if [ "${n}" == "" ] ; then
		seed=${RANDOM}
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_1_2_${seed}.root -covariance
	else
		for((i=0;i<${n};++i)) ; do
			seed=${RANDOM}
			example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_1_2_${seed}.root -covariance
		done
	fi
elif [ "${algo}" == "elisa2" ] ; then 
	# elisa 
	# example/unfold -algorithm 0 -statistical_analysis 100 example/elisa_analysis_test.root
	dfile="example/measurement_2.root"
	dname="meas_x"
	if [ "${n}" == "" ] ; then
		seed=${RANDOM}
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_2_${seed}.root -covariance
	else
		for((i=0;i<${n};++i)) ; do
			seed=${RANDOM}
			example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_2_${seed}.root -covariance
		done
	fi
elif [ "${algo}" == "elisa12" ] ; then 
	# elisa. do both at the same time
	for((i=0;i<${n};++i)) ; do
		seed=${RANDOM}
		dfile="example/measurement_1.root"
		dname="meas_x"
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_1_${seed}.root -covariance
		seed=${RANDOM}
		dfile="example/measurement_2.root"
		dname="meas_x"
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 example/elisa_analysis_experiment_2_${seed}.root -covariance
	done
elif [ "${algo}" == "mle" ] ; then 
	# maximum likelihood
	seed=${RANDOM}
	dfile="example/measurement_1_meas.dat"
	tfile="example/measurement_1_true.dat"
	example/unfold -algorithm 3 -seed ${seed} -format 0 -meas ${dfile} -statistical_analysis 1 example/mle_analysis_${seed}.root
elif [ "${algo}" == "bayes" ] ; then 
	# bayesian iteration
	# dfile="example/measurement_1.root"
	# dname="meas_x"
	dfile="example/measurement_1_meas.dat"
	tfile="example/measurement_1_true.dat"
	seed=${RANDOM}
	nstat="1"
	# example/unfold -algorithm 2 -seed ${seed} -meas ${dfile} ${dname} -iterations 5 -statistical_analysis 1 example/bayes_analysis_${seed}.root
	example/unfold -algorithm 2 -format 0 -meas ${dfile} -iterations 5 -statistical_analysis ${nstat} example/bayes_analysis_${seed}.root
	# example/unfold -algorithm 2 -seed ${seed} -meas ${dfile} ${dname} -iterations 5 -statistical_analysis 1000000 example/bayes_analysis_${seed}.root
elif [ "${algo}" == "1prior2" ] ; then
	for((i=0;i<${n};++i)) ; do
		seed=${RANDOM}
		dfile="example/measurement_1.root"
		dname="meas_x"
		pfile="example/ensemble_2.root"
		pfile="example/scripts/prior2.root"
		pname="elisa_bin"
		echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname}"
		echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname} -statistical_analysis 1000 example/elisa_analysis_experiment_1_prior_2.root"
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -prior ${pfile} ${pname} -statistical_analysis 10000 example/elisa_analysis_experiment_1_prior_2_${seed}.root -covariance
	done
elif [ "${algo}" == "2prior1" ] ; then
	for((i=0;i<${n};++i)) ; do
		seed=${RANDOM}
		dfile="example/measurement_2.root"
		dname="meas_x"
		pfile="example/ensemble_1.root"
		pfile="example/scripts/prior1.root"
		pname="elisa_bin"
		# example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -prior ${pfile} ${pname} 
		example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -prior ${pfile} ${pname} -statistical_analysis 10000 example/elisa_analysis_experiment_2_prior_1_${seed}.root -covariance
	done
fi
