#ifndef UNFOLD_H
#define UNFOLD_H

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

class Unfold {

	public:

	Unfold(int algorithm, const char *name = 0);
	~Unfold();

/* the different algorithm for unfolding */
	enum {
		Elisa,
		FullBayesian,
		BayesianIteration,
		MaximumLikelihood,
	} UnfoldingAlgorithm;

	void set_seed(int seed) { this->seed = seed; };
	int get_seed() { return seed; };
	// void set_dimension(int N) { this->N = N; };
	bool set_true(double *true_bins, int n = 0);
	bool set_meas(double *meas_bins, int n = 0);
	bool set_true(TH1D *h);
	bool set_true(TH2D *h);
	bool set_true(TH3D *h);
	bool set_true(const char *file, const char *name);
	bool set_meas(TH1D *h);
	bool set_meas(TH2D *h);
	bool set_meas(TH3D *h);
	bool set_meas(const char *file, const char *name);
	bool invert_matrix(double **M, int n = 0);
	bool initialize();
	bool fill(double xtrue, double xmeas);
	bool miss(double xtrue);
	bool save(const char *filename, int mode = 0);
	bool add_response_matrix(const char *file, const char *name, double weight = 1.0);
	bool initialize_response_matrix(const char *file, const char *name, double weight = 1.0);
	bool initialize_response_matrix();
	bool run(double *y, double *n, int option = 0);
	bool run(int option = 0);
	bool get_weighted_likelihood_solution(double *y, double *n);
	bool get_bayesian_iterative_solution(double *y, double *n, int niters, double *guess);
	TH1D **create_pdfs(double *n, int nr);
	TH1D **create_pdfs(double **Rinv, TH1D **pdf0, int nr);
	bool get_maximum_likelihood_solution(double *y);
	bool get_maximum_likelihood_solution(double *y, double *n);
	bool bootstrap(double *n = 0);

	private:

	bool ov, uf;
	TH1D *h_x_true, *h_x_meas;
	TH2D *h_x_y_true, *h_x_y_meas;
	TH3D *h_x_y_z_true, *h_x_y_z_meas;
	TH1D *h_efficiency, *h_efficiency_numer, *h_efficiency_denom;
	int N;
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
	double *eff, *deff;
	double *prior, *dprior;
	double epsilon; /* convergence criteria */
	int counter0, max_trials, trials;
	TH2D *response;
	TRandom3 *rndm;
	int nt, nr, iterations, algorithm;
	int dimensions_true, dimensions_meas, nbinsx_true, nbinsy_true, nbinsz_true, nbinsx_meas, nbinsy_meas, nbinsz_meas;
	char *name, *str;
	std::string ofile;
	// TFile *fp;
	bool autosave, M_init;
	int n_response; /* counts number of files containing response matrices */
};

#endif

