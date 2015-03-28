#include "TFile.h"
#include "TH1D.h"

#include "unfold.h"

#include "stdlib.h"

bool debug = false, verbose = false;

int main(int argc, char **argv) {

	std::string name, dfile, dname, rfile, tfile, tname, ofile, sfile, pfile, pname;
	int algorithm = -1, nstat = 0, seed = 0, max_trials = 0, iterations = 5, option = 0;
	double epsilon = 0.001;
	bool covariance = false, bootstrap = false, truth = false, use_prior = false;

/* defaults */
	rfile = "example/response_matrix.dat";
	name = "unfold"; 
	dfile = "example/measurement.root";
	dname = "meas_x";
	// tfile = "example/measurement.root";
	// tname = "true_x";
	// truth = true;

	algorithm = Unfold::MaximumLikelihood;
	algorithm = Unfold::BayesianIteration;
	// algorithm = -1;

	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
		} else if(strcmp("-r", argv[i]) == 0) { rfile = argv[++i]; /* response matrix */
		} else if(strcmp("-meas", argv[i]) == 0) { dfile = argv[++i]; dname = argv[++i]; /* the data */
		} else if(strcmp("-true", argv[i]) == 0) { tfile = argv[++i]; tname = argv[++i]; truth = true; 
		} else if(strcmp("-prior", argv[i]) == 0) { pfile = argv[++i]; pname = argv[++i]; use_prior = true; 
		} else if(strcmp("-trials", argv[i]) == 0) { max_trials = atoi(argv[++i]); 
		} else if(strcmp("-statistical_analysis", argv[i]) == 0) { 
			nstat = atoi(argv[++i]); /* number of pseudo-experiments */
			sfile = argv[++i]; /* output file for pseudo-experiments */
		} else if(strcmp("-seed", argv[i]) == 0) { seed = atoi(argv[++i]); 
		} else if(strcmp("-option", argv[i]) == 0) { option = atoi(argv[++i]); 
		} else if(strcmp("-iterations", argv[i]) == 0) { iterations = atoi(argv[++i]); 
		} else if(strcmp("-algorithm", argv[i]) == 0) { algorithm = atoi(argv[++i]); 
		} else if(strcmp("-o", argv[i]) == 0) { ofile = argv[++i];
		} else if(strcmp("-name", argv[i]) == 0) { name = argv[++i];
		} else if(strcmp("-bootstrap", argv[i]) == 0) { bootstrap = true;
		} else if(strcmp("-covariance", argv[i]) == 0) { covariance = true;
		} else if(strcmp("-epsilon", argv[i]) == 0) { epsilon = atof(argv[++i]); 
		} else if(strcmp("-convergence", argv[i]) == 0) { 
			max_trials = atoi(argv[++i]); 
			epsilon = atof(argv[++i]); 
		} else { printf("unrecognized argument [%s]\n", argv[i]); exit(1); 
		}
	}

	printf("unfolding [%s]\n", dname.c_str());
	printf("convergence criteria = %f\n", epsilon);

	Unfold *unfold = 0;
	if(algorithm < 0) { unfold = new Unfold(Unfold::Elisa, name.c_str()); } 
	else { unfold = new Unfold(algorithm, name.c_str()); }
	if(algorithm == Unfold::BayesianIteration) {
		unfold->set_iterations(iterations);
		unfold->set_progress_report_frequency(10000); /* avoid excessive output to screen */
	} else if(algorithm == Unfold::MaximumLikelihood) {
		unfold->set_progress_report_frequency(10000); /* avoid excessive output to screen */
	} else if(algorithm == Unfold::Elisa) {
		// unfold->set_progress_report_frequency(10000); /* avoid excessive output to screen */
	}
	if(seed) {
		unfold->set_seed(seed);
		printf("seed = %d\n", seed);
	}
	unfold->initialize_response_matrix(rfile.c_str());
	unfold->set_meas(dfile.c_str(), dname.c_str());
	int i, j, k, nt = unfold->get_n_true(), nr = unfold->get_n_meas();
	if(bootstrap) unfold->bootstrap(); /* create new sample bootstrapped from input sample */
	if(tname.length()) unfold->set_true(tfile.c_str(), tname.c_str());
	double *n = unfold->get_meas();
	double *y_true = unfold->get_true();
	double *y = unfold->get_solution();
	for(i=0;i<nr;++i) { printf("MEAS(%d) = %f. TRUE = %f\n", i, n[i], y_true[i]); }

	TH1D **prior = 0;
	if(use_prior) {
		char *str = new char [ pname.length() + 32 ];
		TFile fp(pfile.c_str(), "read");
		prior = new TH1D * [ nt ];
		for(i=0;i<nt;++i) {
			sprintf(str, "%s%d", pname.c_str(), i);
			prior[i] = (TH1D *)fp.Get(str);
			prior[i]->SetDirectory(0);
		}
		fp.Close();
		delete [] str;

		unfold->set_prior(prior);
	}

	unfold->run(option);

	for(i=0;i<nr;++i) { printf("UNFOLDED(%d) = %f\n", i, y[i]); }

	double **R = unfold->get_response_matrix();

	for(int j=0;j<nr;++j) {
		double acc = 0.0;
		for(int i=0;i<nt;++i) {
			acc += y[i] * R[i][j];
		}
		printf("REFOLDED(%d) = %f\n", j, acc);
	}

	if(nstat) {
		printf("evaluating statistical uncertainty with %d trials. saving to %s\n", max_trials, sfile.c_str());
		int source = truth ? Unfold::UseTruth : Unfold::UseUnfolded;
		printf("using %s for statistical analysis\n", truth ? "truth" : "unfolded");
		unfold->statistical_analysis(nstat, source, sfile.c_str(), covariance);
		// unfold->statistical_analysis(nstat, source, sfile.c_str(), true, Unfold::ResponseMatrixVariationUniform, 0.05);
	}
	if(truth) unfold->closure_test();
#if 0
	unfold->write_info();
	unfold->write();

#endif

return 0;

#if 0
	int nt = unfold->get_nt(), nr = unfold->get_nr();
	double *y = new double [ nt ];
	double *y_true = new double [ nt ];
	double *data = new double [ nr ];

	if(covariance) {
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
