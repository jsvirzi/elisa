#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"

#include "stdlib.h"

#include "response_matrix.h"

bool debug = false, verbose = false;

double true_integral(double a, double x_lo, double x_hi) {
	double t_lo = (1.0 + a * x_lo) * TMath::Exp(-a * x_lo);
	double t_hi = (1.0 + a * x_hi) * TMath::Exp(-a * x_hi);
	return (t_lo - t_hi) / (a * a);
}

int main(int argc, char **argv) {
	bool gun = false;
	int nevts = 1000000000, seed = 4357;
	char str[1024];

/* X distribution parameters */
	double x_min = 0.0, x_max = 2.0, x_a = 7.5, x_mean = 0.95, x_width = 0.1;

	std::string ofile("response.root");

	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
		} else if(strcmp("-gun", argv[i]) == 0) { gun = true;
		} else if(strcmp("-o", argv[i]) == 0) { ofile = argv[++i];
		} else if(strcmp("-n", argv[i]) == 0) { nevts = atoi(argv[++i]);
		} else if(strcmp("-a", argv[i]) == 0) { x_a = atof(argv[++i]);
		} else if(strcmp("-seed", argv[i]) == 0) { seed = atoi(argv[++i]);
		} else { printf("unrecognized argument [%s]\n", argv[i]); exit(1); 
		}
	}

	TRandom3 *rndm = new TRandom3(seed);
	seed = rndm->GetSeed();
	gRandom->SetSeed(seed);

	printf("seed = %d\n", seed);
	printf("output file = %s\n", ofile.c_str());

/* X */
	int nbins_x = 10;

	double x_max_generation = 1.5 * x_max;
	TH1D *h_true_x = new TH1D("true_x", "X", nbins_x, x_min, x_max); 
	TH1D *h_meas_x = new TH1D("meas_x", "X", nbins_x, x_min, x_max); 

/* create and initialize response matrix */
	ResponseMatrix *response_matrix = new ResponseMatrix("x");
	response_matrix->set_true(h_true_x);
	response_matrix->set_meas(h_meas_x);
	response_matrix->initialize();

	sprintf(str, "0.95 - 0.45 * TMath::Exp(-2.0 * (x - %f) / %f)", x_min, x_max - x_min);
	TF1 *fxeff = new TF1("fxeff", str, x_min, x_max);

	sprintf(str, "TMath::Gaus(x, %f, %f)", x_mean, x_width);
	TF1 *fxres = new TF1("fxres", str, 0.0, 2.0);
	
	sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
	TF1 *fx0 = new TF1("fx0", str, x_min, x_max_generation); 
#if 0
/* numerical evaluation of the integral */
	double total_integral = fx0->Integral(x_min, x_max_generation);
#else
/* in case we know the exact integral */
	double total_integral = true_integral(x_a, x_min, x_max_generation);
#endif

	printf("total integral of true distribution = %f\n", total_integral);

/* for the purposes of generating the response matrix,
   we go past the boundaries, so as to capture events
   that feed into the measurement region.
   We do not attempt to unfold outside of (x_min, x_max) however */

	for(int i=0;i<(nbins_x+1);++i) {
		double a_min = h_true_x->GetXaxis()->GetBinLowEdge(i+1);
		double a_max = (i < nbins_x) ? h_true_x->GetXaxis()->GetBinUpEdge(i+1) : x_max_generation;
		if(gun) { printf("using particle gun to generate uniform random x on (%f, %f)\n", a_min, a_max);
		} else { printf("using truth distribution to generate random x on (%f, %f)\n", a_min, a_max); }
		sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
		TF1 *fx = new TF1("fx", str, a_min, a_max); 
#if 0
	/* numerical evaluation of the integral */
		double slice_integral = fx->Integral(a_min, a_max);
#else
	/* in case we know the exact integral */
		double slice_integral = true_integral(x_a, a_min, a_max);
#endif
		double weight = slice_integral / total_integral;
		printf("slice integral of true distribution = %f for (%f, %f)\n", slice_integral, a_min, a_max);
		printf("net weight = %g\n", weight);

		response_matrix->set_weight(weight);

		for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
			if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
			double true_x = gun ? rndm->Uniform(a_min, a_max) : fx->GetRandom(); /* random x on truth interval */
			double r_res = fxres->GetRandom(); /* resolution of x */
			double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
		/* generate new measured values of x until we find one inside the volume of measurement */
		/* the efficiency curve will account for missing values of x */
			int i_try = 0, max_tries = 1000;
			while((meas_x < x_min || meas_x >= x_max) && (i_try++ < max_tries)) {
				r_res = fxres->GetRandom(); /* resolution of x */
				meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
			}
			if(i_try == max_tries) {
				printf("timeout error. retry event %d for x in (%f, %f)\n", evt, a_min, a_max);
				continue;
			}

			double eff = fxeff->Eval(true_x); /* efficiency at x */
			double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

		/* efficiency simulation */
			h_true_x->Fill(true_x); /* fill truth unconditionally */
			if(r_eff <= eff) { /* passed efficiency */
				h_meas_x->Fill(meas_x);
				response_matrix->hit(true_x, meas_x); /* fill response matrix */
			} else { /* failed efficiency */
				response_matrix->miss(true_x);
			}

			if(debug) {
				printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
				printf("r_eff=%f. eff=%f\n", r_eff, eff);
				getchar();
			}

			++evt;
		}

		delete fx;
	}

#if 0

	sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
	TF1 *fx = new TF1("fx", str, x_min, x_max_generation); 

	for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
		if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
		double true_x = fx->GetRandom(); /* random x according to distribution */
		double r_res = fxres->GetRandom(); /* resolution of x */
		double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
	/* generate new measured values of x until we find one inside the volume of measurement */
	/* the efficiency curve will account for missing values of x */
		int i_try = 0, max_tries = 1000;
		while((meas_x < x_min || meas_x >= x_max) && (i_try++ < max_tries)) {
			r_res = fxres->GetRandom(); /* resolution of x */
			meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
		}
		if(i_try == max_tries) {
			printf("timeout error. retry event %d for x in (%f, %f)\n", evt, a_min, a_max);
			continue;
		}

		double eff = fxeff->Eval(true_x); /* efficiency at x */
		double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

	/* efficiency simulation */
		h_true_x->Fill(true_x); /* fill truth unconditionally */
		if(r_eff <= eff) { /* passed efficiency */
			h_meas_x->Fill(meas_x);
			response_matrix->hit(true_x, meas_x); /* fill response matrix */
		} else { /* failed efficiency */
			response_matrix->miss(true_x);
		}

		if(debug) {
			printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
			printf("r_eff=%f. eff=%f\n", r_eff, eff);
			getchar();
		}

		++evt;
	}

#endif

	delete fx0;
	delete fxres;
	delete fxeff;

#if 0

	if(gun) {

		for(int i=0;i<(nbins_x+1);++i) {
			double a_min = h_true_x->GetXaxis()->GetBinLowEdge(i+1);
			double a_max = (i < nbins_x) ? h_true_x->GetXaxis()->GetBinUpEdge(i+1) : x_max_generation;
			printf("using particle gun to generate uniform random x on (%f, %f)\n", a_min, a_max);

			for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
				if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
				double true_x = rndm->Uniform(a_min, a_max); /* uniform random x on the truth interval */
				double r_res = fxres->GetRandom(); /* resolution of x */
				double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
			/* generate new measured values of x until we find one inside the volume of measurement */
			/* the efficiency curve will account for missing values of x */
				int i_try = 0, max_tries = 1000;
				while((meas_x < x_min || meas_x >= x_max) && (i_try++ < max_tries)) {
					r_res = fxres->GetRandom(); /* resolution of x */
					meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
				}
				if(i_try == max_tries) {
					printf("timeout error. retry event %d for x in (%f, %f)\n", evt, a_min, a_max);
					continue;
				}

				double eff = fxeff->Eval(true_x); /* efficiency at x */
				double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

			/* efficiency simulation */
				h_true_x->Fill(true_x); /* fill truth unconditionally */
				if(r_eff <= eff) { /* passed efficiency */
					h_meas_x->Fill(meas_x);
					response_matrix->hit(true_x, meas_x); /* fill response matrix */
				} else { /* failed efficiency */
					response_matrix->miss(true_x);
				}

				if(debug) {
					printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
					printf("r_eff=%f. eff=%f\n", r_eff, eff);
					getchar();
				}

				++evt;
			}
		}

	} else {

	/* for the purposes of generating the response matrix,
	   we go past the boundaries, so as to capture events
	   that feed into the measurement region.
	   We do not attempt to unfold outside of (x_min, x_max) however */

		sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
		TF1 *fx = new TF1("fx", str, x_min, x_max_generation); 

		for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
			if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
			double true_x = fx->GetRandom(); /* random x according to distribution */
			double r_res = fxres->GetRandom(); /* resolution of x */
			double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
		/* generate new measured values of x until we find one inside the volume of measurement */
		/* the efficiency curve will account for missing values of x */
			int i_try = 0, max_tries = 1000;
			while((meas_x < x_min || meas_x >= x_max) && (i_try++ < max_tries)) {
				r_res = fxres->GetRandom(); /* resolution of x */
				meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
			}
			if(i_try == max_tries) {
				printf("timeout error. retry event %d for x in (%f, %f)\n", evt, a_min, a_max);
				continue;
			}

			double eff = fxeff->Eval(true_x); /* efficiency at x */
			double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

		/* efficiency simulation */
			h_true_x->Fill(true_x); /* fill truth unconditionally */
			if(r_eff <= eff) { /* passed efficiency */
				h_meas_x->Fill(meas_x);
				response_matrix->hit(true_x, meas_x); /* fill response matrix */
			} else { /* failed efficiency */
				response_matrix->miss(true_x);
			}

			if(debug) {
				printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
				printf("r_eff=%f. eff=%f\n", r_eff, eff);
				getchar();
			}

			++evt;
		}

		delete fx;
	}

#endif

	response_matrix->finalize(ofile);

	return 0;
}

