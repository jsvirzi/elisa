shopt -s expand_aliases

source ~/utils/setup_root.sh
export UNFOLDDIR=/eliza18/atlas/jsvirzi/physics/elisa
export LD_LIBRARY_PATH=${UNFOLDDIR}/lib:${LD_LIBRARY_PATH}
pushd ${UNFOLDDIR} 

exe="example/unfold -algorithm 0"

n=2
for((i=0;i<${n};++i)) ; do
	seed=${RANDOM}

# make pdf for experiment #1
	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/pdf_1/pdf_1.seed_${seed}.root"
	pname="pdf_bin"
	# ${exe} -seed ${seed} -meas ${dfile} ${dname} -create_pdfs 2000000000 ${pfile} ${pname} 

# make pdf for experiment #2
	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/pdf_2/pdf_2.seed_${seed}.root"
	pname="pdf_bin"
	# ${exe} -seed ${seed} -meas ${dfile} ${dname} -create_pdfs 2000000000 ${pfile} ${pname} 

# make pdf for experiment #1 + #2
	dfile="example/measurement_1_2.root"
	dname="meas_x"
	pfile="example/pdf_1_2/pdf_1_2.seed_${seed}.root"
	pname="pdf_bin"
	${exe} -seed ${seed} -meas ${dfile} ${dname} -create_pdfs 2000000000 ${pfile} ${pname} 

# make pdf for experiment #1 using the posterior from #2
	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/scripts/prior2.root"
	pname="elisa_bin"
	ofile="example/pdf_1_prior_2/pdf_1_prior_2.seed_${seed}.root"
	oname="pdf_bin"
	# ${exe} -seed ${seed} -meas ${dfile} ${dname} -prior ${pfile} ${pname} -create_pdfs 2000000000 ${ofile} ${oname} 

# make pdf for experiment #2 using the posterior from #1
	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/scripts/prior1.root"
	pname="elisa_bin"
	ofile="example/pdf_2_prior_1/pdf_2_prior_1.seed_${seed}.root"
	oname="pdf_bin"
	# ${exe} -seed ${seed} -meas ${dfile} ${dname} -prior ${pfile} ${pname} -create_pdfs 2000000000 ${ofile} ${oname} 

# unfold #1 using pdf_1.root. store intermediate solutions to form the posterior
	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/pdf_1.root"
	pname="pdf_bin"
# the ntuple with the intermediate solutions
	nfile="example/experiment_1/elisa_intermediate_solutions_experiment_1.seed_${seed}.root"
# perform statistical analysis
	# sfile="example/experiment_1/elisa_statistical_analysis_experiment_1.seed_${seed}.root"
	# ${exe} -seed ${seed} -o ${nfile} -meas ${dfile} ${dname} -statistical_analysis 10000 ${sfile} -covariance
# no statistical analysis. just the measurement
	# ${exe} -seed ${seed} -o ${nfile} -meas ${dfile} ${dname} -pdf ${pfile} ${pname} -covariance -throws 1000000 -expert

# unfold #2 using pdf_1.root. store intermediate solutions to form the posterior
	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/pdf_2.root"
	pname="pdf_bin"
# the ntuple with the intermediate solutions
	nfile="example/experiment_2/elisa_intermediate_solutions_experiment_2.seed_${seed}.root"
# perform statistical analysis
	# sfile="example/experiment_2/elisa_statistical_analysis_experiment_2.seed_${seed}.root"
	# ${exe} -seed ${seed} -o ${nfile} -meas ${dfile} ${dname} -statistical_analysis 10000 ${sfile} -covariance
# no statistical analysis. just the measurement
	# ${exe} -seed ${seed} -o ${nfile} -meas ${dfile} ${dname} -pdf ${pfile} ${pname} -covariance -throws 1000000 -expert

	continue

	dfile="example/measurement_1.root"
	dname="meas_x"
	ifile="example/experiment_1/elisa_intermediate_results_experiment_1_seed${seed}.root"
	sfile="example/experiment_1/elisa_statistical_analysis_experiment_1_${seed}.root"
	# ${exe} -seed ${seed} -save_intermediate ${ifile} -meas ${dfile} ${dname} -statistical_analysis 10000 ${sfile} -covariance

	dfile="example/measurement_2.root"
	dname="meas_x"
	ifile="example/experiment_2/elisa_intermediate_results_experiment_2_seed${seed}.root"
	sfile="example/experiment_2/elisa_statistical_analysis_experiment_2_${seed}.root"
	# ${exe} -seed ${seed} -save_intermediate ${ifile} -meas ${dfile} ${dname} -statistical_analysis 10000 ${sfile} -covariance

	dfile="example/measurement_1_2.root"
	dname="meas_x"
	ifile="example/experiment_1_2/elisa_intermediate_results_experiment_1_2_seed${seed}.root"
	sfile="example/experiment_1_2/elisa_statistical_analysis_experiment_1_2_${seed}.root"
	# ${exe} -seed ${seed} -save_intermediate ${ifile} -meas ${dfile} ${dname} -statistical_analysis 10000 ${sfile} -covariance

	# 2 as prior for 1
	dfile="example/measurement_1.root"
	dname="meas_x"
	pfile="example/ensemble_2.root"
	pfile="example/scripts/prior2.root"
	pname="elisa_bin"

	# echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname}"
	# echo "example/unfold -algorithm 0 -meas ${dfile} ${dname} -prior ${pfile} ${pname} -statistical_analysis 1000 example/elisa_analysis_experiment_1_prior_2.root"
	# example/unfold -algorithm 0 -meas ${dfile} ${dname} -seed ${seed} -prior ${pfile} ${pname} -statistical_analysis 10000 example/elisa_analysis_experiment_1_prior_2_${seed}.root -covariance

	${exe} -seed ${seed} -meas ${dfile} ${dname} -error_analysis 1000000 example/experiment_1_prior_2/elisa_error_analysis_seed_${seed}.1_prior_2.root -prior ${pfile} ${pname} -statistical_analysis 100000 example/experiment_1_prior_2/elisa_statistical_analysis_experiment_seed_${seed}.1_prior_2.root

	# 1 as a prior for 2
	dfile="example/measurement_2.root"
	dname="meas_x"
	pfile="example/ensemble_1.root"
	pfile="example/scripts/prior1.root"
	pname="elisa_bin"

	${exe} -seed ${seed} -meas ${dfile} ${dname} -error_analysis 1000000 example/experiment_2_prior_1/elisa_error_analysis_seed_${seed}.2_prior_1.root -prior ${pfile} ${pname} -statistical_analysis 100000 example/experiment_2_prior_1/elisa_statistical_analysis_experiment_seed_${seed}.2_prior_1.root
	done
