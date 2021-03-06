#ifndef UNFOLD_H
#define UNFOLD_H

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

bool get_maximum_likelihood_solution(double **R, double *y, int nt, double *n, int nr);

class Unfold {

	public:

	Unfold(int algorithm = NAlgorithms, const char *name = 0);
	~Unfold();

	enum {
		UseUnfolded,
		UseTruth
	};

/* the different algorithm for unfolding */
	enum {
		Elisa,
		FullBayesian,
		BayesianIteration,
		MaximumLikelihood,
		NAlgorithms
	} UnfoldingAlgorithm;

	enum {
		CentralValue,
		InputValue,
		PoissonExtraction,
		PoissonExtractionUsingCentralValue,
		PoissonExtractionUsingTruth,
		Intermediate,
		NBiasNtupleOptions
	} BiasNtupleOptions;

	enum {
		ResponseMatrixVariationNone,
		ResponseMatrixVariationFromHistogram,
		ResponseMatrixVariationUniform
	} ResponseMatrixVariationOptions;

/* file formats */
	enum {
		Binary = 0,
		Text,
		Root,
		NFormats
	};

	void set_seed(int seed) { this->seed = seed; rndm->SetSeed(seed); gRandom->SetSeed(seed); };
	int get_seed() { return seed; };
	// void set_dimension(int N) { this->N = N; };
	bool set_true(double *true_bins, int n = 0);
	bool set_meas(double *meas_bins, int n = 0);
	bool set_true(TH1D *h);
	bool set_true(TH2D *h);
	bool set_true(TH3D *h);
	bool set_true(const char *file, const char *name);
	bool set_true(const char *file);
	bool set_meas(TH1D *h);
	bool set_meas(TH2D *h);
	bool set_meas(TH3D *h);
	bool set_meas(const char *file, const char *name);
	bool set_meas(const char *file);
	bool invert_matrix(double **M, int n = 0);
	bool initialize();
	bool fill(double xtrue, double xmeas);
	bool miss(double xtrue);
	bool save(const char *filename, int mode = 0);
	bool add_response_matrix(const char *file, const char *name, double weight = 1.0);
	// bool set_prior(TH1D **prior);
	// bool initialize_response_matrix(const char *file, const char *name, double weight = 1.0);
	bool initialize_response_matrix(const char *file);
	bool run(double *y = 0, double *n = 0, int option = 0, bool detail = false);
//	bool run(double *y, double *n, TH1D **prior, int option = 0, bool detail = false);
//	bool run(double *y, double *n, double *prior_mean, double *prior_width, int option = 0, bool detail = false);
//	bool run(int option = 0, bool detail = false);
	// bool run(TH1D **prior = 0, int option = 0, bool detail = false);
	// jsv clashes with other run() bool run(double *prior_mean, double *prior_width, int option = 0, bool detail = false);
	bool get_weighted_likelihood_solution(double *y, double *n, bool detail);
	bool get_weighted_likelihood_solution(double *y, double *n, bool detail, TH1D **prior);
	bool get_weighted_likelihood_solution(double *y, double *n, bool detail, bool require_convergence, TH1D **pdf, const char *file, int max_trials);
	// bool get_weighted_likelihood_solution(double *y, double *n, bool detail, double *prior_mean, double *prior_width);
	bool get_bayesian_iterative_solution(double *y, double *n, int niters, double *guess);
	// TH1D **create_pdfs(double *n, int nr);
	TH1D **create_mu_pdfs(bool bias_removal, double *n = 0, const char *file = 0,
		const char *name = 0, int *nbins = 0, double *x_min = 0, double *x_max = 0);
	TH1D **create_theta_pdfs(bool bias_removal, TH1D **prior = 0, int nthrows = 2000000000, const char *file = 0, 
		const char *name = 0, int *nbins = 0, double *x_min = 0, double *x_max = 0); 
	// TH1D **create_pdfs(double **Rinv, TH1D **pdf0, int nr);
	// bool get_maximum_likelihood_solution(double *y);
	bool get_maximum_likelihood_solution(double *y = 0, double *n = 0);
	bool bootstrap(double *n = 0);
	double *get_true();
	double *get_meas();
	double *get_solution();
	double **get_response_matrix();
	int get_n_true() { return nt; };
	int get_n_meas() { return nr; };
	bool statistical_analysis(double *y0, int ntrials, const char *file, 
		bool detail = false, int dR_options = ResponseMatrixVariationNone, double dR_nom = 0.0);
	bool statistical_analysis(int ntrials, int option, const char *file,
		bool detail = false, int dR_options = ResponseMatrixVariationNone, double dR_nom = 0.0);
	bool error_analysis(int ntrials, const char *file);
	bool save_intermediate(const char *file) { intermediate_results_file = file; return true; }
	bool save_intermediate(std::string file) { intermediate_results_file = file; return true; }
	bool save_intermediate(bool flag, const char *file) { 
		if(flag == false) intermediate_results_file = ""; 
		else intermediate_results_file = file; 
		return true;
	}
	bool calculate_response(double *y, double *mu, double **R = 0);
	double closure_test();
	bool set_iterations(int iterations) { this->iterations = iterations; };
	bool get_efficiency(double *eff);
	bool set_progress_report_frequency(int frequency) { progress_report_frequency = frequency; };

	private:

	std::string intermediate_results_file;
	double *make_guess(int option);
	bool write_basic_info(const char *file);
	bool true_uf, true_ov, meas_uf, meas_ov;
	int N;
	int progress_report_frequency; /* how often to report progress */
	bool initialized, verbose, debug;
	int seed; /* random number seed */
	double **R, **dR, **Rinv; /* the response matrix (and its inverse) */
	double **M, **Minv; /* jsv */
	double *true_bins; /* the edges of the bins containing truth occupancies */
	double *meas_bins; /* the edges of the bins containing measured value occupancies */ 

	bool cleanup();

	double *n, *y, *z, *p, *accr, *acct, *mean, *rms, *closure_ratio, *y_true;
	double *A, *B, **C, **cov, **icov, **J;
	double *guess, *bias;
	// double *eff, *deff;
	TH1D **prior;
	// double *prior, *dprior;
	double epsilon; /* convergence criteria */
	int counter0, max_trials, trials;
	TH2D *response;
	TRandom3 *rndm;
	int nt, nr, iterations, algorithm;
	int dimensions_true, dimensions_meas, nbinsx_true, nbinsy_true, nbinsz_true, nbinsx_meas, nbinsy_meas, nbinsz_meas;
	char *name, *str;
	std::string ofile;
	// TFile *fp;
	bool autosave, M_init; /* jsv what is M_init */
	int n_response; /* counts number of files containing response matrices */
};

#endif

