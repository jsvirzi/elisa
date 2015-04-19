shopt -s expand_aliases

source ~/utils/setup_root.sh
export UNFOLDDIR=/eliza18/atlas/jsvirzi/physics/elisa
export LD_LIBRARY_PATH=${UNFOLDDIR}/lib:${LD_LIBRARY_PATH}
pushd ${UNFOLDDIR} 

exe="example/unfold -algorithm 0"

n=1
for((i=0;i<${n};++i)) ; do
	seed=${RANDOM}

	dfile="example/measurement_2.root"
	dname="meas_x"
	${exe} -seed ${seed} -meas ${dfile} ${dname} -create_pdfs 2000000000 example/pdf_2/pdf_2.seed_${seed}.root pdf_bin

	continue

	dfile="example/measurement_1.root"
	dname="meas_x"
	${exe} -seed ${seed} -meas ${dfile} ${dname} -create_pdfs 2000000000 example/pdf_1/pdf_1.seed_${seed}.root pdf_bin

	continue

	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/scripts/prior2.root"
	pname="elisa_bin"
	ofile="example/pdf_1_prior_2/pdf_1_prior_2.seed_${seed}.root"
	oname="pdf_bin"
	${exe} -algorithm 0 -seed ${seed} -meas ${dfile} ${dname} -prior ${pfile} ${pname} -create_pdfs 2000000000 ${ofile} ${oname} 

	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/scripts/prior1.root"
	pname="elisa_bin"
	ofile="example/pdf_2_prior_1/pdf_2_prior_1.seed_${seed}.root"
	oname="pdf_bin"
	${exe} -algorithm 0 -seed ${seed} -meas ${dfile} ${dname} -prior ${pfile} ${pname} -create_pdfs 2000000000 ${ofile} ${oname} 
	continue

	dfile="example/measurement_1.root"
	dname="meas_x"
	ifile="example/experiment_1/elisa_intermediate_results_experiment_1_seed${seed}.root"
	sfile="example/experiment_1/elisa_statistical_analysis_experiment_1_${seed}.root"
	# example/unfold -algorithm 0 -save_intermediate ${ifile} -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 ${sfile} -covariance

	dfile="example/measurement_2.root"
	dname="meas_x"
	ifile="example/experiment_2/elisa_intermediate_results_experiment_2_seed${seed}.root"
	sfile="example/experiment_2/elisa_statistical_analysis_experiment_2_${seed}.root"
	# example/unfold -algorithm 0 -save_intermediate ${ifile} -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 ${sfile} -covariance

	dfile="example/measurement_1_2.root"
	dname="meas_x"
	ifile="example/experiment_1_2/elisa_intermediate_results_experiment_1_2_seed${seed}.root"
	sfile="example/experiment_1_2/elisa_statistical_analysis_experiment_1_2_${seed}.root"
	# example/unfold -algorithm 0 -save_intermediate ${ifile} -meas ${dfile} ${dname} -seed ${seed} -statistical_analysis 10000 ${sfile} -covariance

	# 2 as prior for 1
	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/ensemble_2.root"
	pfile="example/scripts/prior2.root"
	pname="elisa_bin"

	# echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname}"
	# echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname} -statistical_analysis 1000 example/elisa_analysis_experiment_1_prior_2.root"
	# example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -prior ${pfile} ${pname} -statistical_analysis 10000 example/elisa_analysis_experiment_1_prior_2_${seed}.root -covariance

	example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -error_analysis 1000000 example/experiment_1_prior_2/elisa_error_analysis_seed_${seed}.1_prior_2.root -prior ${pfile} ${pname} -statistical_analysis 100000 example/experiment_1_prior_2/elisa_statistical_analysis_experiment_seed_${seed}.1_prior_2.root

	# 1 as a prior for 2
	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/ensemble_1.root"
	pfile="example/scripts/prior1.root"
	pname="elisa_bin"

	example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -error_analysis 1000000 example/experiment_2_prior_1/elisa_error_analysis_seed_${seed}.2_prior_1.root -prior ${pfile} ${pname} -statistical_analysis 100000 example/experiment_2_prior_1/elisa_statistical_analysis_experiment_seed_${seed}.2_prior_1.root
	done
