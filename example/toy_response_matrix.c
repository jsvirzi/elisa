#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"

#include "stdlib.h"

#include "response_matrix.h"

bool debug = false, verbose = false;

int main(int argc, char **argv) {
	int nevts = 1000000000, seed = 4357;
	char str[1024];

/* X distribution parameters */
	double x_min = 0.0, x_max = 2.0, x_a = 7.5, x_mean = 0.95, x_width = 0.1;

	std::string ofile("response.root");

	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
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

	sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
	TF1 *fx = new TF1("fx", str, x_min, x_max); 

	sprintf(str, "TMath::Gaus(x, %f, %f)", x_mean, x_width);
	TF1 *fxres = new TF1("fxres", str, 0.0, 2.0);
	
	sprintf(str, "0.95 - 0.45 * TMath::Exp(-2.0 * (x - %f) / %f)", x_min, x_max - x_min);
	TF1 *fxeff = new TF1("fxeff", str, x_min, x_max);

	TH1D *h_true_x = new TH1D("true_x", "X", nbins_x, x_min, x_max); 
	TH1D *h_meas_x = new TH1D("meas_x", "X", nbins_x, x_min, x_max); 

/* create and initialize response matrix */
	ResponseMatrix *response_matrix_x = new ResponseMatrix("x");
	response_matrix_x->set_true(h_true_x);
	response_matrix_x->set_meas(h_meas_x);
	response_matrix_x->initialize();

	for(int evt=0;evt<nevts;++evt) {
		if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
		double true_x = fx->GetRandom(); /* random x according to distribution */
		double r_res = fxres->GetRandom(); /* resolution of x */
		double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
	/* generate new measured values of x until we find one inside the volume of measurement */
	/* the efficiency curve will account for missing values of x */
		while(meas_x < x_min || meas_x >= x_max) {
			r_res = fxres->GetRandom(); /* resolution of x */
			meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
		}

		double eff = fxeff->Eval(true_x); /* efficiency at x */
		double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

	/* efficiency simulation */
		h_true_x->Fill(true_x); /* fill truth */
		if(r_eff <= eff) { /* passed efficiency */
			h_meas_x->Fill(meas_x);
			response_matrix_x->hit(true_x, meas_x); /* fill response matrix */
		} else { /* failed efficiency */
			response_matrix_x->miss(true_x);
		}

		if(debug) {
			printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
			printf("r_eff=%f. eff=%f\n", r_eff, eff);
			getchar();
		}

	}

	response_matrix_x->finalize(ofile);

	return 0;
}

