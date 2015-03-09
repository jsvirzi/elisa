#include "unfold.h"

#include "stdlib.h"

bool debug = false, verbose = false;

int main(int argc, char **argv) {

	std::string name, dfile, dname, rfile, rname, tfile, tname, ofile;
	int algorithm = -1, nstat = 0, seed = 0, max_trials = 0;
	double epsilon = 0.001;
	bool cov = false, bootstrap = false;

/* defaults */
	rfile = "example/response_matrix.root";
	rname = "x";
	name = "unfold"; 
	dfile = "example/measurement_1.root";
	dname = "meas_x";
	tfile = "example/measurement_1.root";
	tname = "true_x";

	algorithm = Unfold::MaximumLikelihood;
	// algorithm = -1;

	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
		} else if(strcmp("-r", argv[i]) == 0) { rfile = argv[++i]; rname = argv[++i]; /* response matrix */
		} else if(strcmp("-meas", argv[i]) == 0) { dfile = argv[++i]; dname = argv[++i]; /* the data */
		} else if(strcmp("-true", argv[i]) == 0) { tfile = argv[++i]; tname = argv[++i]; /* truth */
		} else if(strcmp("-trials", argv[i]) == 0) { max_trials = atoi(argv[++i]); 
		} else if(strcmp("-stat", argv[i]) == 0) { nstat = atoi(argv[++i]); 
		} else if(strcmp("-seed", argv[i]) == 0) { seed = atoi(argv[++i]); 
		} else if(strcmp("-algorithm", argv[i]) == 0) { algorithm = atoi(argv[++i]); 
		} else if(strcmp("-o", argv[i]) == 0) { ofile = argv[++i];
		} else if(strcmp("-name", argv[i]) == 0) { name = argv[++i];
		} else if(strcmp("-bootstrap", argv[i]) == 0) { bootstrap = true;
		} else if(strcmp("-cov", argv[i]) == 0) { cov = true;
		} else if(strcmp("-epsilon", argv[i]) == 0) { epsilon = atof(argv[++i]); 
		} else if(strcmp("-convergence", argv[i]) == 0) { 
			max_trials = atoi(argv[++i]); 
			epsilon = atof(argv[++i]); 
		} else { printf("unrecognized argument [%s]\n", argv[i]); exit(1); 
		}
	}

	printf("unfolding [%s]\n", dname.c_str());
	printf("convergence criteria = %f\n", epsilon);
	// getchar();

	Unfold *unfold = 0;
	if(algorithm < 0) { unfold = new Unfold(Unfold::Elisa, name.c_str()); } 
	else { unfold = new Unfold(algorithm, name.c_str()); }
	if(seed) unfold->set_seed(seed);
	unfold->initialize_response_matrix(rfile.c_str(), rname.c_str());
	unfold->set_meas(dfile.c_str(), dname.c_str());
	int i, j, k, nt = unfold->get_n_true(), nr = unfold->get_n_meas();
	if(bootstrap) unfold->bootstrap(); /* create new sample bootstrapped from input sample */
	if(tname.length()) unfold->set_true(tfile.c_str(), tname.c_str());
	double *n = unfold->get_meas();
	double *y = unfold->get_true();
	for(i=0;i<nr;++i) { printf("MEAS(%d) = %f. TRUE = %f\n", i, n[i], y[i]); }

	// unfold->set_output_file(ofile.c_str());
	// unfold->set_epsilon(epsilon);
	// unfold->enable_autosave();
	// if(max_trials) unfold->set_max_trials(max_trials);

	unfold->run();

	for(i=0;i<nr;++i) { printf("UNFOLDED(%d) = %f\n", i, y[i]); }

	double **R = unfold->get_response_matrix();

	for(int j=0;j<nr;++j) {
		double acc = 0.0;
		for(int i=0;i<nt;++i) {
			acc += y[i] * R[i][j];
		}
		printf("blah(%d) = %f\n", j, acc);
	}

return 0;

#if 0

	if(nstat) {
		printf("evaluating statistical uncertainty with %d trials\n", max_trials);
		// unfold->evaluate_statistical_uncertainty(0, nstat);
		if(bootstrap) unfold->calculate_bias(nstat, Unfold::UseTruth, "bias_ntuple_elisa.root");
		else unfold->calculate_bias(nstat, Unfold::UseUnfolded, "bias_ntuple_elisa.root");
	}
	if(tname.length()) unfold->closure_test(tfile.c_str(), tname.c_str());
	unfold->write_info();
	unfold->write();

#endif

#if 0
	int nt = unfold->get_nt(), nr = unfold->get_nr();
	double *y = new double [ nt ];
	double *y_true = new double [ nt ];
	double *data = new double [ nr ];

	if(cov) {
		printf("calculating covariance matrix...\n"); getchar();
		unfold->jacobian(0);
	}

	unfold->get_data(data);
	unfold->get_solution(y);
	unfold->get_truth(y_true);
	for(int i=0;i<nt;++i) printf("final comparison. bin(%d) (M)WL=%f. TRUE=%f. DIFF=%f. DATA=%f\n", 
		i, y[i], y_true[i], y[i] / y_true[i], data[i]);

	delete [] y_true;
	delete [] y;
	delete [] data;

#endif

	delete unfold;
	return 0;

#if 0

	for(int i=0;i<nt;++i) printf("final comparison. bin(%d) ML=%f. DATA=%f\n", 
		i, y[i], data[i]);

	unfold->get_data(data);
	unfold->get_solution(y_ml);
	unfold->get_truth(y_true);
	return 0;

	unfold->echo_parameters();

	unfold->write(); /* treat the initial guess as iteration #0 */
	for(int iter=0;iter<iters;++iter) {
		unfold->iterate(1);
		unfold->write();
		if(tname.length()) unfold->closure_test(tfile.c_str(), tname.c_str());
	}

	double *y_true = new double [ nt ];
	double *y_bayes = new double [ nt ];
	double *y_elisa = new double [ nt ];
	unfold->get_solution(y_bayes);
	unfold->get_maximum_likelihood_solution(y_ml);
	for(int i=0;i<nt;++i) printf("final comparison. bin(%d) BAYES=%f ML=%f. DIFF=%f. DATA=%f\n", 
		i, y_bayes[i], y_ml[i], y_bayes[i] / y_ml[i], data[i]);

	unfold->get_weighted_likelihood_solution(y_elisa);
	for(int i=0;i<nt;++i) printf("final comparison. bin(%d) WL=%f. DATA=%f\n", 
		i, y_elisa[i], data[i]);

	if(nstat) {
		printf("evaluating statistical uncertainty\n");
		unfold->evaluate_statistical_uncertainty(0, nstat);
		unfold->write();
	}

	delete [] y_true;

	return 0;

#endif

}