#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "response_matrix.h"
#include "unfold.h"

#include "math.h"

static bool find_inverse(double **A, double **Ainv, int n) {
	int i, j, k, m;
	double det;
	TMatrixD M(n, n);
	for(j=0;j<n;++j) { for(k=0;k<n;++k) M[j][k] = A[j][k]; }
	TMatrixD Minv = M.Invert(&det);
	for(i=0;i<n;++i) { for(j=0;j<n;++j) { Ainv[i][j] = Minv[i][j]; } }
	if(det == 0.0) return false;

#if 0 
	for(i=0;i<n;++i) {
		printf("ROW(%d) ", i);
		for(j=0;j<n;++j) {
			double acc = 0.0;
			for(k=0;k<n;++k) acc += A[i][k] * Ainv[k][j];
			printf("%6.4f ", acc);
		}
		printf("\n");
	}
#endif
#if 1 
	printf("matrix to invert: A=\n");
	for(i=0;i<n;++i) {
		printf("ROW(%2d) ", i);
		for(j=0;j<n;++j) {
			printf("%.3f ", A[i][j]);
		}
		printf("\n");
	}
#endif
#if 1 
	printf("A * A^{-1} = \n");
	for(i=0;i<n;++i) {
		printf("ROW(%2d) ", i);
		for(j=0;j<n;++j) {
			double acc = 0.0;
			for(k=0;k<n;++k) acc += A[i][k] * Ainv[k][j];
			printf("%.f ", fabs(acc));
		}
		printf("\n");
	}
#endif

	return true;
}

/* find solution to R * y = n using Cramer's rule. jsv not tested */
bool get_maximum_likelihood_solution(double **R, double *y, int nt, double *n, int nr) { 
/* this is written as if nr could be different than nt, but nr = nt is required. 
	written this way for transparency */
	int i, j, k;
	TMatrixD Rcopy(nr, nt);
	for(j=0;j<nr;++j) { for(k=0;k<nt;++k) Rcopy[j][k] = R[j][k]; }
	double det, det0 = Rcopy.Determinant();
	if(det0 == 0.0) { for(i=0;i<nt;++i) { y[i] = 0.0; } return false; }
	for(i=0;i<nt;++i) {
		TMatrixD M(nr, nt);
		for(j=0;j<nr;++j) { for(k=0;k<nt;++k) M[j][k] = R[j][k]; }
		for(j=0;j<nr;++j) M[j][i] = n[j];
		det = M.Determinant();
		y[i] = det / det0;
	}

/*
	printf("solution: ");
	for(i=0;i<nt;++i) printf("%6.4f ", y[i]); 
	printf("\n");
	getchar();
 */

	return true;
}

bool Unfold::bootstrap(double *n) {
	int i;
	if(n == 0) for(i=0;i<nr;++i) this->n[i] = rndm->Poisson(n[i]);
	else n[i] = rndm->Poisson(n[i]);
	return true;
}

Unfold::Unfold(int algorithm, const char *name) : 
	true_uf(false), true_ov(false), meas_uf(false), meas_ov(false),
	verbose(false), debug(false), 
	rndm(0), nt(0), nr(0), iterations(0), 
	dimensions_true(0), dimensions_meas(0),
	nbinsx_true(0), nbinsy_true(0), nbinsz_true(0),
	nbinsx_meas(0), nbinsy_meas(0), nbinsz_meas(0),
	R(0), M(0), dR(0), Rinv(0), Minv(0), response(0), 
	n(0), y(0), z(0), p(0), 
	guess(0), 
	accr(0), acct(0), mean(0), rms(0), closure_ratio(0), y_true(0), A(0), B(0), C(0), cov(0), icov(0), J(0), 
	bias(0),
	autosave(false), n_response(0), seed(4357), prior(0),
	progress_report_frequency(0)
	{
	this->algorithm = algorithm;

	rndm = new TRandom3;

	int len = 0;
	if(name) len = strlen(name);
	this->name = new char [ len + 1 ];
	sprintf(this->name, name);
	str = new char [ len + 256 ]; /* general purpose usage when manipulating name */

/* convergence criteria for estimator */
	max_trials = 1000000;
	counter0 = 100;
	epsilon = 0.00001;
}

Unfold::~Unfold() {
	cleanup();
}

bool Unfold::add_response_matrix(const char *file, const char *name, double weight) {
	TFile fp(file, "read");

	if(n_response == 0) {
		sprintf(str, "response_%s", name);
		response = (TH2D *)fp.Get(str);
		response->Scale(weight);
		response->SetDirectory(0);
	} else {
		sprintf(str, "response_%s", name);
		TH2D *h = (TH2D *)fp.Get(str);
		response->Add(h, weight);
	}

	++n_response;

#if 0
	sprintf(str, "efficiency_numer_%s", name);
	h_efficiency_numer = (TH1D *)fp.Get(str);
	h_efficiency_numer->Scale(weight);
	h_efficiency_numer->SetDirectory(0);

	sprintf(str, "efficiency_denom_%s", name);
	h_efficiency_denom = (TH1D *)fp.Get(str);
	h_efficiency_denom->Scale(weight);
	h_efficiency_denom->SetDirectory(0);
#endif

/* jsv. need to deal with overflow and underflow */

#if 0

	sprintf(str, "dimensions_%s", name);
	TH2D *h_dim = (TH2D *)fp.Get(str);
	dimensions_true = h_dim->GetNbinsX();
	dimensions_meas = h_dim->GetNbinsY();

	if(dimensions_true == 1) {
		sprintf(str, "response_distribution_true_%s", name);
		h_x_true = (TH1D *)fp.Get(str);
		h_x_true->SetDirectory(0);
		nbinsx_true = h_x_true->GetNbinsX();
	} else if(dimensions_true == 2) {
		sprintf(str, "response_distribution_true_%s", name);
		h_x_y_true = (TH2D *)fp.Get(str);
		h_x_y_true->SetDirectory(0);
		nbinsx_true = h_x_y_true->GetNbinsX();
		nbinsy_true = h_x_y_true->GetNbinsY();
	} else if(dimensions_true == 3) {
		sprintf(str, "response_distribution_true_%s", name);
		h_x_y_z_true = (TH3D *)fp.Get(str);
		h_x_y_z_true->SetDirectory(0);
		nbinsx_true = h_x_y_z_true->GetNbinsX();
		nbinsy_true = h_x_y_z_true->GetNbinsY();
		nbinsz_true = h_x_y_z_true->GetNbinsZ();
	}

	if(dimensions_meas == 1) {
		sprintf(str, "response_distribution_meas_%s", name);
		h_x_meas = (TH1D *)fp.Get(str);
		h_x_meas->SetDirectory(0);
		nbinsx_meas = h_x_meas->GetNbinsX();
	} else if(dimensions_meas == 2) {
		sprintf(str, "response_distribution_meas_%s", name);
		h_x_y_meas = (TH2D *)fp.Get(str);
		h_x_y_meas->SetDirectory(0);
		nbinsx_meas = h_x_y_meas->GetNbinsX();
		nbinsy_meas = h_x_y_meas->GetNbinsY();
	} else if(dimensions_meas == 3) {
		sprintf(str, "response_distribution_meas_%s", name);
		h_x_y_z_meas = (TH3D *)fp.Get(str);
		h_x_y_z_meas->SetDirectory(0);
		nbinsx_meas = h_x_y_z_meas->GetNbinsX();
		nbinsy_meas = h_x_y_z_meas->GetNbinsY();
		nbinsz_meas = h_x_y_z_meas->GetNbinsZ();
	}

#endif

	return true;
}

#if 0

bool Unfold::initialize_response_matrix(const char *file, const char *name, double weight) {
	if(add_response_matrix(file, name, weight) == false) return false;
	return initialize_response_matrix();
}

#endif

/* after reading in and adding all the response matrices, we need to perform a few more steps before
 * we are ready */
bool Unfold::initialize_response_matrix(const char *file) {

	int i, j;

	ResponseMatrix *response = new ResponseMatrix(file);

	nt = response->get_nt();
	nr = response->get_nr();
	R = new double * [ nt ];
	for(i=0;i<nt;++i) R[i] = new double [ nr ];

	response->get_response_matrix(R);

	M = new double * [ nt ];
	dR = new double * [ nt ];
	// eff = new double [ nt ];
	// deff = new double [ nt ];
	for(i=0;i<nt;++i) M[i] = new double [ nr ];
	for(i=0;i<nt;++i) dR[i] = new double [ nr ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) M[i][j] = 0.0; }

/* the inverse response matrix */
	Rinv = new double * [ nr ];
	for(i=0;i<nr;++i) Rinv[i] = new double [ nt ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) Rinv[i][j] = 0.0; }

	Minv = new double * [ nr ];
	for(i=0;i<nr;++i) Minv[i] = new double [ nt ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) Minv[i][j] = 0.0; }

	y = new double [ nt ]; /* the solution */
	n = new double [ nr ]; /* the data */
	mean = new double [ nt ];
	rms = new double [ nt ];
	bias = new double [ nt ];

/* intermediate variables */
	p = new double [ nt ];
	accr = new double [ nr ];
	acct = new double [ nt ];
	z = new double [ nt ];
	guess = new double [ nt ]; /* the initial guess */
	// prior = new double [ nt ]; /* the initial guess */
	// dprior = new double [ nt ]; /* the initial guess */
	y_true = new double [ nt ]; 
	closure_ratio = new double [ nt ];
	A = new double [ nr ];
	B = new double [ nr ];
	C = new double * [ nr ];
	J = new double * [ nr ];
	cov = new double * [ nr ];
	icov = new double * [ nr ];

	for(i=0;i<nr;++i) C[i] = new double [ nr ];
	for(i=0;i<nr;++i) J[i] = new double [ nr ];
	for(i=0;i<nr;++i) cov[i] = new double [ nr ];
	for(i=0;i<nr;++i) icov[i] = new double [ nr ];

/*** initialize arrays ***/
	for(i=0;i<nr;++i) { for(j=0;j<nr;++j) J[i][j] = C[i][j] = cov[i][j] = icov[i][j] = 0.0; }
	for(i=0;i<nr;++i) { A[i] = B[i] = 0.0; icov[i][i] = 1.0; } /* default is diag(1,1,1,...,1) for covariance matrix */

	for(i=0;i<nt;++i) {
		y[i] = acct[i] = z[i] = mean[i] = rms[i] = closure_ratio[i] = y_true[i] = bias[i] = 0.0;
	}

	for(j=0;j<nr;++j) n[j] = accr[i] = 0.0; 

	printf("Response Matrix = ");
	for(i=0;i<nt;++i) {
		printf("\nROW(%d) = ", i);
		for(j=0;j<nr;++j) {
			printf("%6.3f ", R[i][j]);
		}
	}
	printf("\n");

	printf("inverting %d X %d matrix...\n", nt, nt);
	bool stat = find_inverse(R, Rinv, nt);
	printf("matrix inversion complete...\n");

	delete response;

	return stat;
}

bool Unfold::set_prior(TH1D **prior) {
	if(this->prior) { /* if priors were previously set, clean up */
		for(int i=0;i<nt;++i) delete this->prior[i]; 
	} else {
		this->prior = new TH1D * [ nt ];
	}
	for(int i=0;i<nt;++i) {
		this->prior[i] = new TH1D(*prior[i]);
		this->prior[i]->SetDirectory(0);
	}
	return true;
}

bool Unfold::cleanup() {

	if(str) delete [] str;

	if(M) { for(int i=0;i<nt;++i) delete [] M[i]; delete [] M; M = 0; }

	if(R) { for(int i=0;i<nt;++i) delete [] R[i]; delete [] R; R = 0; }

	if(Rinv) {
		for(int i=0;i<nt;++i) delete [] Rinv[i];
		delete [] Rinv;
		Rinv = 0;
	}

	if(dR) {
		for(int i=0;i<nt;++i) delete [] dR[i];
		delete [] dR;
		dR = 0;
	}

	if(prior) {
		for(int i=0;i<nt;++i) delete prior[i]; 
		delete [] prior; 
	}

	// if(eff) { delete [] eff; eff = 0; }
	// if(deff) { delete [] deff; deff = 0; }
	if(y) { delete [] y; y = 0; }
	if(guess) { delete [] guess; guess = 0; }
	// if(prior) { delete [] prior; prior = 0; }
	// if(dprior) { delete [] dprior; dprior = 0; }
	if(y_true) { delete [] y_true; y_true = 0; }
	if(n) { delete [] n; n = 0; }
	if(z) { delete [] z; z = 0; }
	if(p) { delete [] p; p = 0; }
	if(accr) { delete [] accr; accr = 0; }
	if(acct) { delete [] acct; acct = 0; }
	if(mean) { delete [] mean; mean = 0; }
	if(rms) { delete [] rms; rms = 0; }
	if(bias) { delete [] bias; bias = 0; }
	if(closure_ratio) { delete [] closure_ratio; closure_ratio = 0; }
	if(A) { delete [] A; A = 0; }
	if(B) { delete [] B; B = 0; }
	if(C) {
		for(int i=0;i<nr;++i) delete [] C[i];
		delete [] C;
		C = 0;
	}
	if(cov) {
		for(int i=0;i<nr;++i) delete [] cov[i];
		delete [] cov;
		cov = 0;
	}
	if(icov) {
		for(int i=0;i<nr;++i) delete [] icov[i];
		delete [] icov;
		icov = 0;
	}

	return true;

}

double *Unfold::get_true() { return y_true; }

double *Unfold::get_solution() { return y; }

double **Unfold::get_response_matrix() { return R; }

bool Unfold::set_true(TH1D *h) {
	// h_x_true = new TH1D(*h);
	// sprintf(str, "%s", str);
	// h_x_true->SetName(str);
	// h_x_true->SetDirectory(0);
	int nbins = h->GetNbinsX();
	if(h->GetBinContent(0) != 0.0) true_uf = true;
	if(h->GetBinContent(nbins+1) != 0.0) true_ov = true;
	int start = true_uf ? 0 : 1;
	int stop = true_ov ? (nbins + 1) : nbins;
printf("stop = %d. start = %d. true uf/ov = %s/%s\n", stop, start, true_uf ? "true" : "false", true_ov ? "true" : "false");
	int nread = 1 + stop - start;
	if(nt != nread) {
		printf("mismatch number of bins. read %d. expected %d.\n", nread, nr); 
		return false;
	}
	for(int i=0;i<=stop;++i) y_true[i] = h->GetBinContent(i + start);
	// printf("setting true = \n");
	// for(int i=0;i<nt;++i) printf("%.1f\n", y_true[i]);
	return true;
}

#if 0
bool Unfold::set_true(TH1D *h) {
	sprintf(str, "%s", str);
	h_x_true = new TH1D(*h);
	h_x_true->SetName(str);
	h_x_true->SetDirectory(0);
	int start = true_uf ? 0 : 1;
	for(int i=0;i<nt;++i) y[i] = h_x_true->GetBinContent(i + start);
	return true;
}
#endif

bool Unfold::set_true(TH2D *h) {
#if 0
	sprintf(str, "%s", str);
	h_x_y_true = new TH2D(*h);
	h_x_y_true->SetName(str);
	h_x_y_true->SetDirectory(0);
	for(int i=0;i<nt;++i) y_true[i] = h_x_y_true->GetBinContent(i);
#endif
	for(int i=0;i<nt;++i) y_true[i] = h->GetBinContent(i);
	return true;
}

bool Unfold::set_true(TH3D *h) {
#if 0
	sprintf(str, "%s", str);
	h_x_y_z_true = new TH3D(*h);
	h_x_y_z_true->SetName(str);
	h_x_y_z_true->SetDirectory(0);
	for(int i=0;i<nt;++i) y_true[i] = h_x_y_z_true->GetBinContent(i);
#endif
	for(int i=0;i<nt;++i) y_true[i] = h->GetBinContent(i);
	return true;
}

bool Unfold::set_true(double *y, int N) {
	if(this->N != N) return false;
	for(int i=0;i<nr;++i) this->y_true[i] = y[i];
	return true;
}

bool Unfold::set_true(const char *file, const char *name) {
	TFile fp(file, "read");

#if 0
	TFile fp(file, "read");
	sprintf(str, "dimensions_%s", name);
	TH2D *h_dim = (TH2D *)fp.Get(str);
	if(dimensions_true == 0) {
		dimensions_true = h_dim->GetNbinsX();
		printf("dimensions of response matrix (TRUE=%d. MEAS=%d)\n", dimensions_true, dimensions_meas);
	} else {
		int dimensions_true = h_dim->GetNbinsX();
		if(this->dimensions_true != dimensions_true) {
			printf("inconsistency in response matrix dimensions for TRUE (%d vs %d)\n", 
				this->dimensions_true, dimensions_true);
			return false;
		}
	}
#endif

	dimensions_true = 1; /* jsv */

	if(dimensions_true == 1) {
		TH1D *h = (TH1D *)fp.Get(name); 
		int nbins = h->GetNbinsX();
		bool uf = (h->GetBinContent(0) != 0.0), ov = (h->GetBinContent(nbins+1) != 0.0);
		if(uf == false && true_uf) {
			printf("UF != 0 in TRUE histogram but UF is used in response matrix");
			return false;
		}
		if(ov == false && true_ov) {
			printf("OV != 0 in TRUE histogram but OV is used in response matrix");
			return false;
		}
		nbins += uf ? 1 : 0;
		nbins += ov ? 1 : 0;
		if(nbins != nt) { 
			printf("mismatch in number of true bins. expected %d. found %d\n", nt, nbins);
			return false;
		}
		return set_true(h);

#if 0
	} else if(dimensions_true == 2) {
		TFile fp(file, "read");
		TH2D *h = (TH2D *)fp.Get(name); 
		nbinsx_true = h->GetNbinsX();
		nbinsy_true = h->GetNbinsY();
		if((nbinsx_true != nx) || (nbinsy_true != ny)) { 
			printf("mismatch in number of true bins. expected (%d X %d). found (%d X %d)\n", 
				nbinsx_true, nbinsy_true, nx, ny);
			return false;
		}
		return set_true(h);
	} else if(dimensions_true == 3) {
		TFile fp(file, "read");
		TH3D *h = (TH3D *)fp.Get(name); 
		int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
		if((nbinsx_true != nx) || (nbinsy_true != ny) || (nbinsz_true != nz)) { 
			printf("mismatch in number of true bins. expected (%d X %d X %d). found (%d X %d X %d)\n", 
				nbinsx_true, nbinsy_true, nbinsz_true, nx, ny, nz);
			return false;
		}
		return set_true(h);
#endif
	}

	return false;
}

double *Unfold::get_meas() {
	return n;
}

bool Unfold::set_meas(double *n, int N) {
	if(this->N != N) return false;
	for(int i=0;i<nr;++i) this->n[i] = n[i];
	return true;
}

bool Unfold::set_meas(TH1D *h) {
	// h_x_meas = new TH1D(*h);
	// sprintf(str, "%s", str);
	// h_x_meas->SetName(str);
	// h_x_meas->SetDirectory(0);
	int nbins = h->GetNbinsX();
	if(h->GetBinContent(0) != 0.0) meas_uf = true;
	if(h->GetBinContent(nbins+1) != 0.0) meas_ov = true;
	int start = meas_uf ? 0 : 1;
	int stop = meas_ov ? (nbins + 1) : nbins;
printf("stop = %d. start = %d. meas uf/ov = %s/%s\n", stop, start, meas_uf ? "true" : "false", meas_ov ? "true" : "false");
	int nread = 1 + stop - start;
	if(nr != nread) {
		printf("mismatch number of bins. read %d. expected %d.\n", nread, nr); 
		return false;
	}
	for(int i=0;i<=stop;++i) n[i] = h->GetBinContent(i + start);
	// printf("setting meas = \n");
	// for(int i=0;i<nr;++i) printf("%.1f\n", n[i]);
	return true;
}

bool Unfold::set_meas(TH2D *h) {
	sprintf(str, "%s", str);
#if 0
	h_x_y_meas = new TH2D(*h);
	h_x_y_meas->SetName(str);
	h_x_y_meas->SetDirectory(0);
#endif
	for(int i=0;i<nr;++i) n[i] = h->GetBinContent(i);
	return true;
}

bool Unfold::set_meas(TH3D *h) {
#if 0
	sprintf(str, "%s", str);
	h_x_y_z_meas = new TH3D(*h);
	h_x_y_z_meas->SetName(str);
	h_x_y_z_meas->SetDirectory(0);
#endif
	for(int i=0;i<nr;++i) n[i] = h->GetBinContent(i);
	return true;
}

bool Unfold::set_meas(const char *file, const char *name) {
	TFile fp(file, "read");
#if 0
	sprintf(str, "dimensions_%s", name);
	TH2D *h_dim = (TH2D *)fp.Get(str);
	if(dimensions_meas == 0) {
		dimensions_meas = h_dim->GetNbinsY();
		printf("dimensions of response matrix (TRUE=%d. MEAS=%d)\n", dimensions_true, dimensions_meas);
	} else {
		int dimensions_meas = h_dim->GetNbinsY();
		if(this->dimensions_meas != dimensions_meas) {
			printf("inconsistency in response matrix dimensions for MEAS (%d vs %d)\n", 
				this->dimensions_meas, dimensions_meas);
			return false;
		}
	}
#endif

	dimensions_meas = 1; /* jsv */

	if(dimensions_meas == 1) {
		TH1D *h = (TH1D *)fp.Get(name); 
		int nbins = h->GetNbinsX();
		bool uf = (h->GetBinContent(0) != 0.0), ov = (h->GetBinContent(nbins+1) != 0.0);
		if(uf == false && meas_uf) {
			printf("UF != 0 in MEAS histogram but UF is used in response matrix");
			return false;
		}
		if(ov == false && meas_ov) {
			printf("OV != 0 in MEAS histogram but OV is used in response matrix");
			return false;
		}
		nbins += uf ? 1 : 0;
		nbins += ov ? 1 : 0;
		if(nbins != nr) { 
			printf("mismatch in number of meas bins. expected %d. found %d\n", nr, nbins);
			return false;
		}
		return set_meas(h);
#if 0
/* jsv */
	} else if(dimensions_meas == 2) {
		TH2D *h = (TH2D *)fp.Get(name); 
		int nx = h->GetNbinsX(), ny = h->GetNbinsY();
		if((nbinsx_meas != nx) || (nbinsy_meas != ny)) { 
			printf("mismatch in number of meas bins. expected (%d X %d). found (%d X %d)\n", 
				nbinsx_meas, nbinsy_meas, nx, ny);
			return false;
		}
		return set_meas(h);
	} else if(dimensions_meas == 3) {
		TH3D *h = (TH3D *)fp.Get(name); 
		int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
		if((nbinsx_meas != nx) || (nbinsy_meas != ny) || (nbinsz_meas != nz)) { 
			printf("mismatch in number of meas bins. expected (%d X %d X %d). found (%d X %d X %d)\n", 
				nbinsx_meas, nbinsy_meas, nbinsz_meas, nx, ny, nz);
			return false;
		}
		return set_meas(h);
#endif
	}

	return false;
}

double *Unfold::make_guess(int option) {
	double *guess = new double [ nt ];
	if(option == 0) {
		guess = new double [ nt ];
		for(int i=0;i<nt;++i) guess[i] = 1.0;
	} else if(option == 1) {
		guess = new double [ nt ];
		get_maximum_likelihood_solution(guess, n);
	} else if(option == 2) {
		guess = new double [ nt ];
		for(int i=0;i<nt;++i) guess[i] = 1.0;
	} else if(option == 3) {
		guess = new double [ nt ];
		for(int i=0;i<nt;i+=2) { guess[i] = 1.0; }
		for(int i=1;i<nt;i+=2) { guess[i] = 11.0; }
	}
	return guess;
}

/* jsv TODO lots of overlap between run() methods. fix */

bool Unfold::run(int option, bool detail) {
	bool stat = false;
	if(algorithm == BayesianIteration) {
		double *guess = make_guess(option);
		stat = get_bayesian_iterative_solution(y, n, iterations, guess);
		if(guess) delete [] guess;
	} else if(algorithm == MaximumLikelihood) {
		stat = get_maximum_likelihood_solution(y, n);
	} else if(algorithm == Elisa) {
		double *ntemp = new double [ nr ];
		for(int i=0;i<nt;++i) ntemp[i] = (n[i] > 1.0) ? (n[i] - 1.0) : 0.0;
		stat = get_weighted_likelihood_solution(y, ntemp, detail);
		delete [] ntemp;
	}
	return stat;
}

/* input:
 * 	n is data
 * output:
 * 	y is unfolded result
 */
bool Unfold::run(double *y, double *n, int option, bool detail) {
	bool stat = false;
	if(algorithm == BayesianIteration) {
		double *guess = make_guess(option); 
		stat = get_bayesian_iterative_solution(y, n, iterations, guess);
		if(guess) delete [] guess;
	} else if(algorithm == MaximumLikelihood) {
		stat = get_maximum_likelihood_solution(y, n);
	} else if(algorithm == Elisa) {
		double *ntemp = new double [ nr ];
		for(int i=0;i<nr;++i) ntemp[i] = (n[i] > 1.0) ? (n[i] - 1.0) : 0.0;
		if(prior == 0) stat = get_weighted_likelihood_solution(y, ntemp, detail);
		else stat = get_weighted_likelihood_solution(y, ntemp, detail, prior);
		delete [] ntemp;
	}
	return stat;
}

bool Unfold::get_bayesian_iterative_solution(double *y, double *n, int niters, double *guess) {

	int i, j;
	double *eff = new double [ nt ];
	get_efficiency(eff);

	if(guess) { for(i=0;i<nt;++i) y[i] = guess[i]; } 
	else { for(i=0;i<nt;++i) y[i] = 1.0; }

	if(guess && verbose) { 
		printf("initial distribution for bayesian iteration = ...\n");
		for(i=0;i<nt;++i) printf("bin(%d) guess = %f\n", i, y[i]);
		// if(debug) 
		//	getchar();
	}

	if(verbose) {
		printf("iters = %d\n", niters);
		for(i=0;i<nt;++i) printf("start: y(%d) = %f\n", i, y[i]); 
	}

	for(int iter=0;iter<niters;++iter) {

		bool verbose = false;
#if 0
		if(verbose) {
			printf("begin iteration %d\n", iterations+1);
			for(i=0;i<nt;++i) printf("y[%d] = %f\n", i, y[i]); 
		}
#endif

	/* normalize input y */
		double acc = 0.0;
		for(j=0;j<nt;++j) acc += y[j];
		for(j=0;j<nt;++j) p[j] = y[j] / acc;

		for(i=0;i<nr;++i) {
			double acc = 0.0;
			for(j=0;j<nt;++j) { acc += p[j] * R[j][i]; }
			accr[i] = acc;
		}

		for(i=0;i<nt;++i) {
			double acc = 0.0;
			for(j=0;j<nr;++j) {
				double a = (accr[j] == 0.0) ? 0.0 : (p[i] * R[i][j] * n[j] / accr[j]);
				acc += a; 
			}
			if(eff[i] > 0.0) z[i] = acc / eff[i];
			else z[i] = 0.0;
		}

		for(i=0;i<nt;++i) { y[i] = z[i]; }

#if 0
		if(verbose) {
			for(i=0;i<nt;++i) printf("iteration(%d) y[%d] = %f. eff = %f. n = %f\n", 
				iterations, i, y[i], eff[i], n[i]);
		}
		++iterations;
#endif

	}

//	for(i=0;i<nt;++i) printf("final: y(%d) = %f\n", i, y[i]); 
//	getchar();

	delete [] eff;

	return true;

}

bool Unfold::get_weighted_likelihood_solution(double *y, double *n, bool detail) {
	int i, j, k, trial;
	int counter = counter0;
	double *ytemp = new double [ nt ];
	double *ycand = new double [ nr ];
	double *ysave = new double [ nr ];

	trials = -1; /* nothing good has happened yet */

	TH1D **pdf = create_pdfs(n, nr);
	// create_pdfs(Rinv, pdf, nr); /* jsv */

	// printf("get_elisa(%x, %x) with NT=%d NR=%d\n", y, n, nt, nr);

#define JSVDEBUG 0
// #if JSVDEBUG 
/* write them out for now */
	// TFile fp("pdf_elisa.root", "recreate");
	// for(i=0;i<nr;++i) { pdf[i]->Write(); }
	// fp.Write();
	// fp.Close();
// #endif

	double *mu = new double [ nr ];
	double *v = new double [ nr ];
	double a = 0;
	for(i=0;i<nr;++i) { v[i] = 0.0; }
	for(i=0;i<nt;++i) { ysave[i] = ycand[i] = 0.0; }
	for(i=0;i<nr;++i) { A[i] = B[i] = 0.0; for(j=0;j<nr;++j) C[i][j] = 0.0; } 
	bool converge = false, flag;
	for(trial=0;trial<max_trials;++trial) {

		if(progress_report_frequency && ((trial % progress_report_frequency) == 0)) {
			printf("%d pass. %d total\n", trial, max_trials);
		}

	/* draw random mu from the PDF */
		for(i=0;i<nr;++i) { mu[i] = pdf[i]->GetRandom(); }
		flag = get_maximum_likelihood_solution(ytemp, mu); /* ytemp = mu X Rinv. for Theta-function test */
		if(flag == false) continue;

		// for(i=0;i<nr;++i) printf("BIN(%d) : MU = %f. Y = %f. ACC = %f\n", i, mu[i], ytemp[i], ycand[i]);
		// getchar();

		flag = false; /* assume theta-function yields 1 */
		for(i=0;i<nt;++i) { if(ytemp[i] <= 0.0) { flag = true; break; } } /* Theta-function yields 0 */
		if(flag) continue; /* Theta = 0 */

		for(i=0;i<nr;++i) v[i] += mu[i];
		a += 1.0;
		++trial;

		if(detail) {
			for(i=0;i<nr;++i) {
				double s = mu[i];
				double t = TMath::Log(s);
				A[i] += s; 
				B[i] += t;
				for(j=0;j<nr;++j) C[j][i] += mu[j] * t; 
			}
		}

		for(j=0;j<nr;++j) {
			ysave[j] = ycand[j]; /* save previous state */
			ycand[j] = v[j] / a; /* new state */
		}; 

		converge = true; /* assume convergence */
		for(j=0;j<nr;++j) {
			double y_old = ysave[j];
			double y_new = ycand[j];
			double y_ave = 0.5 * (y_new + y_old);
			if(fabs(y_old - y_new) > epsilon * y_ave) { 
				converge = false; /* no convergence yet */
				counter = counter0; /* reset counter */
				break; /* no need to continue after decision about no convergence */
			}
		}; 

		if(converge) {
			--counter; /* count down */
			if(counter == 0) { /* enough successive trials have converged */ 
				printf("convergence criteria reached with %d trials!\n", trial); 
				get_maximum_likelihood_solution(y, ycand);
				// for(i=0;i<nt;++i) y[i] = ycand[i]; /* jsv */
				if(detail) {
				/* jsv. is assumption about the y-convergence valid for A, B and C? */
					for(j=0;j<nr;++j) { A[j] = A[j] / a; }
					for(j=0;j<nr;++j) { B[j] = B[j] / a; }
					for(i=0;i<nr;++i) { for(j=0;j<nr;++j) { C[i][j] = C[i][j] / a; } }
					for(i=0;i<nr;++i) { for(j=0;j<nr;++j) { J[i][j] = C[i][j] - A[i] * B[j]; } }
				}
				break;
			}
		}

	}

#if 0
	printf("J = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", J[i][j]);
		}
		printf("\n");
	}

	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			cov[i][j] = 0.0;
			for(int k=0;k<nr;++k) {
				for(int m=0;m<nr;++m) {
					cov[i][j] += icov[k][m] * J[i][k] * J[j][m];
				}
			}
		}
	}

	printf("INPUT COVARIANCE = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", icov[i][j]);
		}
		printf("\n");
	}

	printf("OUTPUT COVARIANCE = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", cov[i][j]);
		}
		printf("\n");
	}
#endif

	delete [] mu;
	delete [] v;
	delete [] ytemp;
	delete [] ycand;
	delete [] ysave;
	for(i=0;i<nr;++i) { delete pdf[i]; }
	delete [] pdf;

	trials = trial; /* return the number of trials required for convergence */

	return converge && (trials < max_trials);

}

TH1D **Unfold::create_pdfs(double *n, int nr) {

	// printf("creating PDFs\n");
/* create PDFs for random number generation */
	char s[1024];
	int i, j, bin, nbins = 250;
	double x, acc, width, x_min, x_max, *a = new double [ nbins ];
	TH1D **pdf = new TH1D * [ nr ];
	for(i=0;i<nr;++i) {
		TH1D *h;

		sprintf(s, "pdfn%d", i);
		width = 4.0 * sqrt(n[i] + 1.0);
		x_min = n[i] - width;
		x_max = n[i] + width;
		if(x_min <= 0.0) x_min = 0.01; /* jsv TODO figure if configurable. stay away from 0? */
		h = new TH1D(s, "PDF", nbins, x_min, x_max); 
		h->SetDirectory(0);
		// printf("i=%d. x = (%f, %f, %f)\n", i, n[i], x_min, x_max); 
		acc = 0.0;
		for(j=0;j<nbins;++j) {
			bin = j + 1;
			x = h->GetXaxis()->GetBinCenter(bin);
			a[j] = TMath::Poisson(n[i], x);
			acc += a[j];
		}

		for(j=0;j<nbins;++j) {
			bin = j + 1;
			a[j] /= acc;
			h->SetBinContent(bin, a[j]);
		}

		pdf[i] = h;

	}

	delete [] a;

	// printf("PDF creation complete\n");

	return pdf;

}

/* jsv. Can we get rid of ytemp below? */
TH1D **Unfold::create_pdfs(double **Rinv, TH1D **pdf0, int nr) {
	printf("create_pdfs(%x, %x, %d)\n", Rinv, pdf0, nr);
	char str[1024], file[1024];
	int i, j, nt = nr, trial, max_trials = 100000000;
	double *ntemp = new double [ nr ];
	double *ytemp = new double [ nt ];
	double *sumy0 = new double [ nt ];
	double *sumy1 = new double [ nt ];
	double *sumy2 = new double [ nt ];
	TH1D **pdf = new TH1D * [ nt ]; 
	TFile *fp;

#if 0
	fp = new TFile("pdfy.root", "read");
	for(i=0;i<nt;++i) {
		sprintf(str, "pdfy%d", i);
		pdf[i] = (TH1D *)fp->Get(str); 
		pdf[i]->SetDirectory(0);
		int bin0 = pdf[i]->GetXaxis()->FindBin(0.0);
		for(j=0;j<=bin0;++j) pdf[i]->SetBinContent(j, 0.0);
	}

	return pdf;
#endif

	for(i=0;i<nt;++i) { sumy0[i] = sumy1[i] = sumy2[i] = 0.0; } /* initialize accumulators */

	printf("figuring out mean and rms\n");
	max_trials = 1000000;
	for(trial=0;trial<max_trials;++trial) {
		for(i=0;i<nr;++i) ntemp[i] = pdf0[i]->GetRandom(); 
		for(i=0;i<nt;++i) {
			double acc = 0.0;
			for(j=0;j<nr;++j) acc += Rinv[i][j] * ntemp[j];
			ytemp[i] = acc;
		}
		for(i=0;i<nt;++i) {
			double acc = ytemp[i];	
			sumy0[i] += 1.0;
			sumy1[i] += acc;
			sumy2[i] += (acc * acc);
		}
	}

	for(i=0;i<nt;++i) {
		double mean = sumy1[i] / sumy0[i];
		double rms = sqrt(sumy2[i] / sumy0[i] - mean * mean);
		printf("bin(%d). mean=%f. rms=%f\n", i, mean, rms);
		sprintf(str, "pdfy%d", i);
		pdf[i] = new TH1D(str, "PDF(y)", 400, mean - 4.0 * rms, mean + 4.0 * rms); 
		pdf[i]->SetDirectory(0);
	}

	max_trials = 10000000;
	printf("making the pdfs\n");
	for(trial=0;trial<max_trials;++trial) {
		if(trial && ((trial % 1000000) == 0)) printf("%d / %d\n", trial, max_trials);
		for(i=0;i<nr;++i) ntemp[i] = pdf0[i]->GetRandom(); 
		for(i=0;i<nt;++i) {
			double acc = 0.0;
			for(j=0;j<nr;++j) acc += Rinv[i][j] * ntemp[j];
			ytemp[i] = acc;
		}
		for(i=0;i<nt;++i) {
			pdf[i]->Fill(ytemp[i]);
		}
	}

	fp = new TFile("pdfy.root", "recreate");
	for(i=0;i<nt;++i) pdf[i]->Write();
	fp->Write(); 
	fp->Close(); 

	for(i=0;i<nt;++i) {
		int bin0 = pdf[i]->GetXaxis()->FindBin(0.0);
		for(j=0;j<=bin0;++j) pdf[i]->SetBinContent(j, 0.0);
	}

	delete [] sumy0;
	delete [] sumy1;
	delete [] sumy2;
	delete [] ntemp;
	delete [] ytemp;
}

#if 0
/* find solution to R * y = n using Cramer's rule */
bool Unfold::get_maximum_likelihood_solution(double *y) { 
/* this is written as if nr could be different than nt, but nr = nt is required. 
	written this way for transparency */
	int i, j, k;
	TMatrixD Rcopy(nr, nt);
	for(j=0;j<nt;++j) { for(k=0;k<nr;++k) Rcopy[j][k] = R[j][k]; }
	double det, det0 = Rcopy.Determinant();
	if(det0 == 0.0) { for(i=0;i<nt;++i) { y[i] = 0.0; } return false; }
	for(i=0;i<nt;++i) {
		TMatrixD M(nr, nt);
		for(j=0;j<nr;++j) { for(k=0;k<nt;++k) M[j][k] = R[j][k]; }
		for(j=0;j<nr;++j) M[j][i] = n[j];
		det = M.Determinant();
		y[i] = det / det0;
	}
	return true;
}
#endif

/* find solution to R * y = n using matrix inversion */
bool Unfold::get_maximum_likelihood_solution(double *y, double *n) { 
	trials = 1; /* number of tries is always 1 for MLE. here for consistency with other methods */
	if(n == 0) n = this->n;
/* this is written as if nr could be different than nt, but nr = nt is required. 
	written this way for transparency */
	int i, j;
	for(i=0;i<nt;++i) {
		double acc = 0.0;
		for(j=0;j<nr;++j) acc += n[j] * Rinv[j][i];
		y[i] = acc;
	}
	return true;
}

bool Unfold::statistical_analysis(int ntrials, int option, const char *ntuple, bool detail, int dR_options, double dR_nominal) {
	bool stat = false;
	if(option == UseUnfolded) {
		double *y_temp = new double [ nt ]; /* need to buffer y because it is modified by the algorithm */
		for(int i=0;i<nt;++i) y_temp[i] = y[i];
		stat = statistical_analysis(y_temp, ntrials, ntuple, detail, dR_options, dR_nominal);
		delete [] y_temp;
	} else if(option == UseTruth) {
		stat = statistical_analysis(y_true, ntrials, ntuple, detail, dR_options, dR_nominal);
	}
	return stat;
}

bool Unfold::get_efficiency(double *eff) {
	for(int i=0;i<nt;++i) {
		double acc = 0.0;
		for(int j=0;j<nr;++j) acc += R[i][j];
		eff[i] = acc;
	}
	return true;
}

bool Unfold::write_basic_info(const char *file) {

	TTree *tree = 0;
	TFile *fp = 0;
	int i, j, k, dim = nt, dim2 = nt * nr;
	double *Rinv = new double [ dim2 ]; 
	double *R = new double [ dim2 ]; 
	double *dR = new double [ dim2 ]; 
	double *ntrue = new double [ nr ];
	double *ytrue = new double [ nt ];
	double *n = new double [ nr ];
	double *y = new double [ nt ];
	double *eff = new double [ nt ];

	get_efficiency(eff);

	fp = new TFile(file, "update");
	tree = new TTree("info", "info");
	tree->Branch("nt", &nr, "nt/I");
	tree->Branch("nr", &nr, "nr/I");
	tree->Branch("dim", &dim, "dim/I");
	tree->Branch("dim2", &dim2, "dim2/I");
	tree->Branch("R", R, "R[dim2]/D");
	tree->Branch("dR", dR, "dR[dim2]/D");
	tree->Branch("Rinv", Rinv, "Rinv[dim2]/D");
	tree->Branch("eff", eff, "eff[nt]/D");
	tree->Branch("y", y, "y[nt]/D");
	tree->Branch("n", n, "n[nr]/D");
	tree->Branch("ytrue", ytrue, "ytrue[nt]/D");
	tree->Branch("ntrue", ntrue, "ntrue[nr]/D");
	// for(i=0;i<nt;++i) eff[i] = this->eff[i]; 
	// for(i=0;i<nt;++i) deff[i] = this->deff[i]; 
	for(i=0;i<nt;++i) ytrue[i] = this->y_true[i]; 
	for(i=0;i<nr;++i) { ntrue[i] = 0.0; for(j=0;j<nt;++j) { ntrue[i] += ytrue[j] * this->R[i][j]; } } 
	for(i=0;i<nt;++i) y[i] = this->y[i]; 
	for(i=0;i<nr;++i) n[i] = this->n[i]; 
	for(i=k=0;i<nt;++i) { for(j=0;j<nr;++j) { R[k++] = this->R[i][j]; } }
	for(i=k=0;i<nt;++i) { for(j=0;j<nr;++j) { dR[k++] = this->dR[i][j]; } }
	for(i=k=0;i<nt;++i) { for(j=0;j<nr;++j) { Rinv[k++] = this->Rinv[i][j]; } }
	tree->Fill();
	tree->Write();
	fp->Write();
	fp->Close();
	delete fp;

	delete [] R;
	delete [] dR;
	delete [] Rinv;
	delete [] ntrue;
	delete [] ytrue;
	delete [] n;
	delete [] y;
	delete [] eff;

	return true;
}

/* uses input distribution "y" and response matrix "R" to calculate expected response into "mu".
 * the user can input a different response matrix. If R == 0, default response matrix is used */
bool Unfold::calculate_response(double *y, double *mu, double **R) {
	if(R == 0) R = this->R;
	for(int i=0;i<nr;++i) {
		double acc = 0.0;
		for(int j=0;j<nt;++j) acc += y[j] * R[j][i];
		mu[i] = acc;
	}
}

bool Unfold::statistical_analysis(double *y0, int ntrials, const char *file, bool detail, int dR_options, double dR_nominal) {
	int i, j, k, trial, status, throws, dim2;
	double *mu = new double [ nr ];
	double *ntemp = new double [ nr ];
	double *ytemp = new double [ nt ];
	double *sumx0 = new double [ nt ];
	double *sumx1 = new double [ nt ];
	double *sumx2 = new double [ nt ];
	double *atemp = 0, *btemp = 0, *ctemp = 0, *cov = 0, *J = 0; 
	double **R_save = 0, **Rinv_save = 0, *Rtemp = 0;
	if(dR_options != ResponseMatrixVariationNone) {
		R_save = R;
		Rinv_save = Rinv;
		R = new double * [ nt ];
		Rinv = new double * [ nr ];
		for(i=0;i<nt;++i) R[i] = new double [ nr ];
		for(i=0;i<nr;++i) Rinv[i] = new double [ nt ];
		for(i=0;i<nt;++i) { for(j=0;j<nr;++j) R[i][j] = this->R[i][j]; }
		bool stat = find_inverse(R, Rinv, nt);
		if(stat == false) {
			for(i=0;i<nt;++i) delete [] R[i];
			for(i=0;i<nr;++i) delete [] Rinv[i];
			delete [] R;
			delete [] Rinv;
			R = R_save;
			Rinv = Rinv_save;
			return false;
		}
	}

	if(detail) {
		int size = nt;
		atemp = new double [ size ];
		bzero(atemp, sizeof(double) * size);
		btemp = new double [ size ];
		bzero(btemp, sizeof(double) * size);
		size = nt * nr; /* used as contiguous array */
		ctemp = new double [ size ];
		bzero(ctemp, sizeof(double) * size);
		cov = new double [ size ];
		bzero(cov, sizeof(double) * size);
		J = new double [ size ];
		bzero(J, sizeof(double) * size);
		Rtemp = new double [ size ];
		bzero(Rtemp, sizeof(double) * size);
	}
	dim2 = nt * nr;

	write_basic_info(file);

	TFile *fp = new TFile(file, "update");
	TTree *tree = new TTree("extraction", "extraction");
	tree->Branch("throws", &throws, "throws/I");
	tree->Branch("status", &status, "status/I");
	tree->Branch("dim", &nr, "dim/I");
	tree->Branch("dim2", &dim2, "dim2/I");
	tree->Branch("n", ntemp, "n[dim]/D");
	tree->Branch("y", ytemp, "y[dim]/D");
	if(detail) {
		tree->Branch("a", atemp, "a[dim]/D");
		tree->Branch("b", btemp, "b[dim]/D");
		tree->Branch("c", ctemp, "c[dim2]/D");
		tree->Branch("R", Rtemp, "R[dim2]/D");
		tree->Branch("cov", cov, "cov[dim2]/D");
		tree->Branch("J", J, "J[dim2]/D");
/* jsv TODO save Rinv too */
	}

/* put away the central value. use status = 0 */
	for(i=0;i<nt;++i) ytemp[i] = y[i];
	for(i=0;i<nr;++i) ntemp[i] = n[i];
	status = CentralValue;
	throws = 0;
	tree->Fill();

	calculate_response(y0, mu); /* final working value of mu */

/* the input values used. status = 1 */
	for(i=0;i<nt;++i) ytemp[i] = y0[i];
	for(i=0;i<nr;++i) ntemp[i] = mu[i];
	status = InputValue;
	throws = 0;
	tree->Fill();

	calculate_response(y0, mu); /* working value of mu */

/* check for validity of "input" data. it must be at least 1, definitely not 0 */
	bool valid = false;
	while(!valid) {
/* jsv. TODO loop limit? */
		valid = true;
		for(i=0;i<nr;++i) {
			if(mu[i] < 0.001) { /* some token amount above zero */
				valid = false;
				break;
			}
		}
	}

	for(trial=0;trial<ntrials;++trial) {
		if(progress_report_frequency && ((trial % progress_report_frequency) == 0))
			printf("trial %d / %d\n", trial, ntrials);

	/* Poisson extraction on the expected value in the detector. want at least 2 entries in each bin */
	/* jsv. TODO loop limit? */
		for(i=0;i<nr;++i) {
			ntemp[i] = rndm->Poisson(mu[i]);
			while(ntemp[i] < 1.001) ntemp[i] = rndm->Poisson(mu[i]);
		}

		if(dR_options == ResponseMatrixVariationUniform) {
			for(i=0;i<nt;++i) {
				for(j=0;j<nr;++j) {
					R[i][j] = R_save[i][j] * rndm->Gaus(1.0, dR_nominal);
				}
			}
		}
		bool flag = run(ytemp, ntemp, 0, detail);
		status = flag ? PoissonExtraction : NBiasNtupleOptions; /* success or failure */
		throws = this->trials; /* get the number of times dice was thrown to obtain convergence */
		if(detail) {
			for(i=k=0;i<nr;++i) { for(j=0;j<nr;++j) { ctemp[k++] = C[i][j]; } }
			for(i=k=0;i<nr;++i) { for(j=0;j<nr;++j) { cov[k++] = this->cov[i][j]; } } /* jsv where is this used? */
			for(i=k=0;i<nr;++i) { for(j=0;j<nr;++j) { J[k++] = this->J[i][j]; } }
			for(i=0;i<nr;++i) { atemp[i] = A[i]; btemp[i] = B[i]; } 
		}
		tree->Fill();
		if(autosave) tree->AutoSave();
	}

	fp->cd();
	tree->Write();
	fp->Write();
	fp->Close();
	delete fp;

// jsv not required?	delete tree;

	if(R_save) {
		for(i=0;i<nt;++i) delete [] R[i];
		for(i=0;i<nr;++i) delete [] Rinv[i];
		delete [] R;
		delete [] Rinv;
		R = R_save;
		Rinv = Rinv_save;
	}

	if(detail) {
		delete [] J;
		delete [] cov;
		delete [] atemp;
		delete [] btemp;
		delete [] ctemp;
		delete [] Rtemp;
	}
	delete [] mu;
	delete [] ntemp;
	delete [] ytemp;
	delete [] sumx0;
	delete [] sumx1;
	delete [] sumx2;

}

double Unfold::closure_test() {
	for(int i=0;i<nt;++i) {
		closure_ratio[i] = (y_true[i] != 0.0) ? (y[i] / y_true[i]) : 0.0;
	}
	return true;
}

#if 0

/* after reading in and adding all the response matrices, we need to perform a few more steps before
 * we are ready */
bool Unfold::initialize_response_matrix(void) {

	int i, j, binx, biny;

	int nx = response->GetNbinsX(), ny = response->GetNbinsY();

printf("jsv. nx = %d. ny = %d\n", nx, ny);

#if 0
	int nt_expected = 0;
	if(dimensions_true == 1) {
		nt_expected = (nbinsx_true + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	} else if(dimensions_true == 2) {
		nt_expected = (nbinsx_true + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsy_true + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	} else if(dimensions_true == 3) {
		nt_expected = (nbinsx_true + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsy_true + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsz_true + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	}
	if(nt != nt_expected) { printf("mismatch in expected(%d) bins for TRUE(%d)\n", nt_expected, nt); return false; }

	int nr_expected = 0;
	if(dimensions_meas == 1) {
		nr_expected = (nbinsx_meas + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	} else if(dimensions_meas == 2) {
		nr_expected = (nbinsx_meas + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsy_meas + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	} else if(dimensions_meas == 3) {
		nr_expected = (nbinsx_meas + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsy_meas + ((ov ? 1 : 0) + (uf ? 1 : 0))) * (nbinsz_meas + ((ov ? 1 : 0) + (uf ? 1 : 0)));
	}
	if(nr != nr_expected) { printf("mismatch in expected(%d) bins for RECO(%d)\n", nr_expected, nr); return false; }

printf("jsv. expected %d %d\n", nt_expected, nr_expected);

#endif
	// if(overflow == false) { nt -= 2; nr -= 2; }

/* TRUE: if any underflow bins are detected as != 0.0, then underflow bins are to be used */
	true_uf = false, true_ov = false;
	for(i=0;i<=(nx+1);++i) {
		if(response->GetBinContent(i, 0) != 0.0) { 
			true_uf = true;
			break;
		}
	}

/* if any overflow bins are detected as != 0.0, then overflow bins are to be used */
	for(i=0;i<=(nx+1);++i) {
		if(response->GetBinContent(i, ny+1) != 0.0) {
			true_ov = true;
			break;
		}
	}
	nt = nx + (true_uf ? 1 : 0) + (true_ov ? 1 : 0);
printf("TRUE: UF=%s. OV=%s\n", true_uf ? "true" : "false", true_ov ? "true" : "false");

/* MEAS: if any underflow bins are detected as != 0.0, then underflow bins are to be used */
	meas_uf = false, meas_ov = false;
	for(i=0;i<=(ny+1);++i) {
		if(response->GetBinContent(0, i) != 0.0) { 
			meas_uf = true;
			break;
		}
	}

/* if any overflow bins are detected as != 0.0, then overflow bins are to be used */
	for(i=0;i<=(ny+1);++i) {
		if(response->GetBinContent(nx+1, i) != 0.0) { 
			meas_ov = true;
			break;
		}
	}
	nr = ny + (meas_uf ? 1 : 0) + (meas_ov ? 1 : 0);
printf("MEAS: UF=%s. OV=%s\n", meas_uf ? "true" : "false", meas_ov ? "true" : "false");

	M = new double * [ nt ];
	R = new double * [ nt ];
	dR = new double * [ nt ];
	// eff = new double [ nt ];
	// deff = new double [ nt ];
	for(i=0;i<nt;++i) M[i] = new double [ nr ];
	for(i=0;i<nt;++i) R[i] = new double [ nr ];
	for(i=0;i<nt;++i) dR[i] = new double [ nr ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) R[i][j] = M[i][j] = 0.0; }

	int x_start = true_uf ? 0 : 1;
	int y_start = meas_uf ? 0 : 1;

	if(n_response == 1) {

		for(i=0,binx=x_start;i<nt;++i,++binx) {
			for(j=0,biny=y_start;j<nr;++j,++biny) {
				R[i][j] = response->GetBinContent(binx, biny);
				// printf("R[%d][%d] = %f\n", i, j, R[i][j]);
			}
		}

#if 0
/* jsv add support for multiple response matrices */
	} else if(n_response != 1) {

		// int start = uf ? 2 : 1;
		for(i=0,binx=start;i<nt;++i,++binx) {
			double ad = h_efficiency_denom->GetBinContent(binx);
			double an = h_efficiency_numer->GetBinContent(binx);
			double ed = h_efficiency_denom->GetBinError(binx);
			double en = h_efficiency_numer->GetBinError(binx);
			eff[i] = deff[i] = 0.0; 
			if(ad != 0.0 && an != 0.0) { /* more or less normal */ 
				eff[i] = an / ad;
				deff[i] = eff[i] * sqrt(pow(en/an,2.0)+pow(ed/ad,2.0));
			}
			printf("efficiency(bin=%d) = %f +/- %f\n", binx, eff[i], deff[i]);
		}


		for(j=0,biny=start;j<nt;++j,++biny) {
			double acc = 0.0;
			double a = h_efficiency_denom->GetBinContent(biny);
			for(i=0,binx=start;i<nr;++i,++binx) {
				R[i][j] = response->GetBinContent(binx, biny);
				printf("R[%d][%d] = %f\n", i, j, R[i][j]);
				acc += R[i][j];
				dR[i][j] = response->GetBinError(binx, biny);
			}
			for(i=0;i<nr;++i) { R[i][j] /= a; }
			double old_eff = eff[j];
			eff[j] = 0.0;
			for(i=0;i<nr;++i) { eff[j] += R[i][j]; }
			printf("compare(%d->%d) effs = (%f, %f). diff = %f\n", j, biny, old_eff, eff[j], eff[j] / old_eff);
			
		}
#endif

	}

/* the inverse response matrix */
	Rinv = new double * [ nr ];
	for(i=0;i<nr;++i) Rinv[i] = new double [ nt ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) Rinv[i][j] = 0.0; }

	Minv = new double * [ nr ];
	for(i=0;i<nr;++i) Minv[i] = new double [ nt ];
	for(i=0;i<nt;++i) { for(j=0;j<nr;++j) Minv[i][j] = 0.0; }

	y = new double [ nt ]; /* the solution */
	n = new double [ nr ]; /* the data */
	mean = new double [ nt ];
	rms = new double [ nt ];
	bias = new double [ nt ];

/* intermediate variables */
	p = new double [ nt ];
	accr = new double [ nr ];
	acct = new double [ nt ];
	z = new double [ nt ];
	guess = new double [ nt ]; /* the initial guess */
	prior = new double [ nt ]; /* the initial guess */
	dprior = new double [ nt ]; /* the initial guess */
	y_true = new double [ nt ]; 
	closure_ratio = new double [ nt ];
	A = new double [ nr ];
	B = new double [ nr ];
	C = new double * [ nr ];
	J = new double * [ nr ];
	cov = new double * [ nr ];
	icov = new double * [ nr ];

	for(i=0;i<nr;++i) C[i] = new double [ nr ];
	for(i=0;i<nr;++i) J[i] = new double [ nr ];
	for(i=0;i<nr;++i) cov[i] = new double [ nr ];
	for(i=0;i<nr;++i) icov[i] = new double [ nr ];

/*** initialize arrays ***/
	for(i=0;i<nr;++i) { for(j=0;j<nr;++j) J[i][j] = C[i][j] = cov[i][j] = icov[i][j] = 0.0; }
	for(i=0;i<nr;++i) { A[i] = B[i] = 0.0; icov[i][i] = 1.0; } /* default is diag(1,1,1,...,1) for covariance matrix */

	for(i=0;i<nt;++i) {
		y[i] = acct[i] = z[i] = mean[i] = rms[i] = closure_ratio[i] = y_true[i] = bias[i] = 0.0;
	}

	for(j=0;j<nr;++j) n[j] = accr[i] = 0.0; 

	printf("inverting %d X %d matrix...\n", nt, nt);
	bool stat = find_inverse(R, Rinv, nt);
	printf("matrix inversion complete...\n");

	printf("Response Matrix = ");
	for(i=0;i<nt;++i) {
		printf("\nROW(%d) = ", i);
		for(j=0;j<nr;++j) {
			printf("%6.3f ", R[i][j]);
		}
	}
	printf("\n");
//	getchar();

	return stat;
}

#endif

bool Unfold::get_weighted_likelihood_solution(double *y, double *n, bool detail, TH1D **prior) {
	int i, j, k, trial;
	int counter = counter0;
	double *ytemp = new double [ nt ];
	double *ycand = new double [ nr ];
	double *ysave = new double [ nr ];
	double *theta = new double [ nt ];
	double *mu = new double [ nr ];
	double *v1 = new double [ nt ];
	double *v0 = new double [ nt ]; /* jsv. don't need array for this */

/* for normalizing the Poisson weights */
	double *prior_mean = new double [ nt ];
	for(i=0;i<nt;++i) prior_mean[i] = prior[i]->GetMean();
	double *prior_prob = new double [ nr ];
	calculate_response(prior_mean, mu);
	for(j=0;j<nr;++j) prior_prob[j] = TMath::Poisson(n[j], mu[j]); 
// for(i=0;i<nt;++i) printf("prior %d has mean %f and prob = %f\n", i, prior_mean[i], prior_prob[i]);

	for(i=0;i<nt;++i) { v0[i] = v1[i] = ysave[i] = ycand[i] = 0.0; }
	for(i=0;i<nr;++i) { A[i] = B[i] = 0.0; for(j=0;j<nr;++j) C[i][j] = 0.0; } 
	bool converge = false, flag;
	trials = -1; /* nothing good has happened yet */
	for(trial=0;trial<max_trials;++trial) {

		if(progress_report_frequency && ((trial % progress_report_frequency) == 0)) {
			printf("%d pass. %d total\n", trial, max_trials);
		}

	/* draw random solution from the priors PDF */
		for(i=0;i<nt;++i) { theta[i] = prior[i]->GetRandom(); }
		calculate_response(theta, mu);
		double acc = 1.0;
		for(j=0;j<nr;++j) {
			double t = TMath::Poisson(n[j], mu[j]);
// printf("bin %d probs = %f %f coming from THETA=%f MU=%f N=%f\n", j, t, prior_prob[j], theta[j], mu[j], n[j]);
			acc = acc * t / prior_prob[j]; 
		}
// printf("ACC=%f\n", acc);
		for(i=0;i<nt;++i) {
// printf("THETA(%d) = %f\n", i, theta[i]);
			v1[i] += theta[i] * acc;
			v0[i] += acc;
		}
		++trial;

#if 0
/* jsv TODO got to figure this out */
		if(detail) {
			for(i=0;i<nr;++i) {
				double s = mu[i];
				double t = TMath::Log(s);
				A[i] += s; 
				B[i] += t;
				for(j=0;j<nr;++j) C[j][i] += mu[j] * t; 
			}
		}
#endif

		for(j=0;j<nr;++j) {
			ysave[j] = ycand[j]; /* save previous state */
			ycand[j] = v1[j] / v0[j]; /* new state */
// printf("cand(%d) = %f = %f / %f\n", j, ycand[j], v1[j], v0[j]);
		}; 

// getchar();

		converge = true; /* assume convergence */
		for(j=0;j<nr;++j) {
			double y_old = ysave[j];
			double y_new = ycand[j];
			double y_ave = 0.5 * (y_new + y_old);
			if(fabs(y_old - y_new) > epsilon * y_ave) { 
				converge = false; /* no convergence yet */
				counter = counter0; /* reset counter */
				break; /* no need to continue after decision about no convergence */
			}
		}; 

		if(converge) {
			--counter; /* count down */
			if(counter == 0) { /* enough successive trials have converged */ 
				printf("convergence criteria reached with %d trials!\n", trial); 
				for(i=0;i<nt;++i) y[i] = ycand[i];
#if 0
/* jsv TODO got to figure this out */
				if(detail) {
				/* jsv. is assumption about the y-convergence valid for A, B and C? */
					for(j=0;j<nr;++j) { A[j] = A[j] / a; }
					for(j=0;j<nr;++j) { B[j] = B[j] / a; }
					for(i=0;i<nr;++i) { for(j=0;j<nr;++j) { C[i][j] = C[i][j] / a; } }
					for(i=0;i<nr;++i) { for(j=0;j<nr;++j) { J[i][j] = C[i][j] - A[i] * B[j]; } }
				}
#endif

				break;
			}
		}

	}

#if 0
	printf("J = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", J[i][j]);
		}
		printf("\n");
	}

	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			cov[i][j] = 0.0;
			for(int k=0;k<nr;++k) {
				for(int m=0;m<nr;++m) {
					cov[i][j] += icov[k][m] * J[i][k] * J[j][m];
				}
			}
		}
	}

	printf("INPUT COVARIANCE = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", icov[i][j]);
		}
		printf("\n");
	}

	printf("OUTPUT COVARIANCE = \n");
	for(i=0;i<nr;++i) {
		for(j=0;j<nr;++j) {
			printf("%5.3f ", cov[i][j]);
		}
		printf("\n");
	}
#endif

	delete [] theta;
	delete [] mu;
	delete [] v0;
	delete [] v1;
	delete [] ytemp;
	delete [] ycand;
	delete [] ysave;
	delete [] prior_mean;
	delete [] prior_prob;

	trials = trial; /* return the number of trials required for convergence */

	return converge && (trials < max_trials);

}

