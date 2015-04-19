#include "TFile.h"
#include "TH1D.h"

#include "unfold.h"

#include "stdlib.h"

bool debug = false, verbose = false;

/* jsv. read from file */
#define NBINS 12
struct {
	int bin, nbins;
	double x_min, x_max;
} pdf_info[NBINS+1] = {
	{0, 250, 1759000, 1778000},
	{1, 250, 1426000, 1442000},
	{2, 250, 546200, 556000},
	{3, 250, 172300, 179000},
	{4, 250, 49000, 53000},
	{5, 250, 12700, 14900},
	{6, 250, 3000, 4200},
	{7, 250, 550, 1250},
	{8, 250, 90, 400},
	{9, 250, 0, 140},
	{10, 250, 0, 40},
	{11, 250, 0, 25},
	{0, 0, 0, 0}
};

int main(int argc, char **argv) {

	std::string name, dfile, dname, rfile, tfile, tname, ofile, sfile, pfile, pname, efile, ifile,
		pdf_ofile, pdf_oname;
	int algorithm = -1, nstat = 0, seed = 0, max_trials = 0, iterations = 5, option = 0, nerrm = 0;
	int pdf_throws = 0;
	double epsilon = 0.001;
	bool covariance = false, bootstrap = false, truth = false, use_prior = false, use_pdf = false,
		save_intermediate = false, create_pdfs = false, limits_known = false, expert = false;

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
		} else if(strcmp("-expert", argv[i]) == 0) { expert = true;
		} else if(strcmp("-r", argv[i]) == 0) { rfile = argv[++i]; /* response matrix */
		} else if(strcmp("-meas", argv[i]) == 0) { dfile = argv[++i]; dname = argv[++i]; /* the data */
		} else if(strcmp("-true", argv[i]) == 0) { tfile = argv[++i]; tname = argv[++i]; truth = true; 
		} else if(strcmp("-save_intermediate", argv[i]) == 0) { ifile = argv[++i]; save_intermediate = true; 
		} else if(strcmp("-prior", argv[i]) == 0) { pfile = argv[++i]; pname = argv[++i]; use_prior = true; 
		} else if(strcmp("-pdf", argv[i]) == 0) { pfile = argv[++i]; pname = argv[++i]; use_pdf = true; 
		} else if(strcmp("-create_pdfs", argv[i]) == 0) { 
			pdf_throws = atoi(argv[++i]);
			pdf_ofile = argv[++i]; 
			pdf_oname = argv[++i];
			create_pdfs = true; 
		} else if(strcmp("-trials", argv[i]) == 0) { max_trials = atoi(argv[++i]); 
		} else if(strcmp("-error_analysis", argv[i]) == 0) { 
			nerrm= atoi(argv[++i]); /* number of pseudo-experiments */
			efile = argv[++i]; /* output file for pseudo-experiments */
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

	if(use_pdf && use_prior) {
		printf("conflicting options chosen. pdf = true and prior = true\n");
		return 1;
	}

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
		unfold->save_intermediate(save_intermediate, ifile.c_str());
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
	double *ytrue = unfold->get_true();
	for(i=0;i<nr;++i) { printf("MEAS(%d) = %f. TRUE = %f\n", i, n[i], ytrue[i]); }

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

		// unfold->set_prior(prior);
	}

	if(create_pdfs) {
		limits_known = true;
		double *x_min = 0, *x_max = 0;
		int *nbins = 0;
		bool bias_removal = true; /* subtract 1 from measured spectrum */
		printf("creating pdfs with %d throws. file = %s. name = %s\n", pdf_throws, pdf_ofile.c_str(), pdf_oname.c_str());
		if(limits_known) {
			x_min = new double [ NBINS ];
			x_max = new double [ NBINS ];
			nbins = new int [ NBINS];
			for(i=0;i<NBINS;++i) {
				x_min[i] = pdf_info[i].x_min;
				x_max[i] = pdf_info[i].x_max;
				nbins[i] = pdf_info[i].nbins;
				printf("LIMIT(%d) = %d [%f, %f]\n", i, nbins[i], x_min[i], x_max[i]);
			}
		}
		unfold->create_theta_pdfs(bias_removal, prior, pdf_throws, pdf_ofile.c_str(), pdf_oname.c_str(),
			nbins, x_min, x_max);
		return 0;
	}

	TH1D **pdf = 0;
	if(use_pdf) {
		char *str = new char [ pname.length() + 32 ];
		TFile fp(pfile.c_str(), "read");
		pdf= new TH1D * [ nt ];
		for(i=0;i<nt;++i) {
			sprintf(str, "%s%d", pname.c_str(), i);
			pdf[i] = (TH1D *)fp.Get(str);
			pdf[i]->SetDirectory(0);
		}
		fp.Close();
		delete [] str;
	}

	if(expert) {
		double *ntemp = new double [ nr ];
		double *ytemp = new double [ nt ];
		for(int i=0;i<nt;++i) ntemp[i] = (n[i] > 1.0) ? (n[i] - 1.0) : 0.0;
		bool detail = true; /* want A, B and C */
		bool require_convergence = false; /* we will do our own compilation later */
		printf("running in expert mode\n");
		unfold->get_weighted_likelihood_solution(ytemp, ntemp, detail, require_convergence, pdf, ofile.c_str());
		delete [] ntemp;
		delete [] ytemp;
	} else {
		unfold->run(option);
	}

	double *y = unfold->get_solution();
	unfold->save_intermediate(false, 0); /* turn off saving intermediate results */

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

	if(nerrm) {
		printf("evaluating error analysis with %d trials. saving to %s\n", nerrm, efile.c_str());
		unfold->error_analysis(nerrm, efile.c_str());
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
	double *ytrue = new double [ nt ];
	double *data = new double [ nr ];

	if(covariance) {
		printf("calculating covariance matrix...\n"); getchar();
		unfold->jacobian(0);
	}

	unfold->get_data(data);
	unfold->get_solution(y);
	unfold->get_truth(ytrue);
	for(int i=0;i<nt;++i) printf("final comparison. bin(%d) (M)WL=%f. TRUE=%f. DIFF=%f. DATA=%f\n", 
		i, y[i], ytrue[i], y[i] / ytrue[i], data[i]);

	delete [] ytrue;
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
	unfold->get_truth(ytrue);
	return 0;

	unfold->echo_parameters();

	unfold->write(); /* treat the initial guess as iteration #0 */
	for(int iter=0;iter<iters;++iter) {
		unfold->iterate(1);
		unfold->write();
		if(tname.length()) unfold->closure_test(tfile.c_str(), tname.c_str());
	}

	double *ytrue = new double [ nt ];
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

	delete [] ytrue;

	return 0;

#endif

}
