#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"

#include "stdlib.h"

bool debug = false, verbose = false;

double true_integral(double a, double x_lo, double x_hi) {
	double t_lo = (1.0 + a * x_lo) * TMath::Exp(-a * x_lo);
	double t_hi = (1.0 + a * x_hi) * TMath::Exp(-a * x_hi);
	return (t_lo - t_hi) / (a * a);
}

int main(int argc, char **argv) {
	int nevts = 1000000, seed = 5743;
	char str[1024];

/* X distribution parameters */
	double x_min = 0.0, x_max = 2.0, x_a = 7.5, x_mean = 0.98, x_width = 0.05;

	std::string ofile, name = "x";

	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
		} else if(strcmp("-o", argv[i]) == 0) { ofile = argv[++i];
		} else if(strcmp("-name", argv[i]) == 0) { name = argv[++i];
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

	TFile *fp = new TFile(ofile.c_str(), "update");

/* X */
	int nbins_x = 10;

	double x_max_generation = 1.5 * x_max;
	TH1D *h_true = new TH1D("true_x", "X", nbins_x, x_min, x_max); 
	TH1D *h_meas = new TH1D("meas_x", "X", nbins_x, x_min, x_max); 

	sprintf(str, "0.95 - 0.45 * TMath::Exp(-2.0 * (x - %f) / %f)", x_min, x_max - x_min);
	TF1 *fxeff = new TF1("fxeff", str, x_min, x_max);

	sprintf(str, "TMath::Gaus(x, %f, %f)", x_mean, x_width);
	TF1 *fxres = new TF1("fxres", str, 0.0, 2.0);
	
	sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
	TF1 *fx = new TF1("fx", str, x_min, x_max_generation); 

	for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
		if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
		double true_x = fx->GetRandom(); /* random x on truth interval */
		double r_res = fxres->GetRandom(); /* resolution of x */
		double meas_x = r_res * true_x; /* reconstructed x is resolution * truth x */
		double eff = fxeff->Eval(true_x); /* efficiency at x */
		double r_eff = rndm->Uniform(); /* random dice to see if we pass efficiency */

	/* efficiency simulation */
		h_true->Fill(true_x); /* fill truth unconditionally */
		if(r_eff <= eff) h_meas->Fill(meas_x); /* passed efficiency */

		if(debug) {
			printf("true_x = %f. meas_x = %f. r_res = %f\n", true_x, meas_x, r_res);
			printf("r_eff=%f. eff=%f\n", r_eff, eff);
			getchar();
		}

		++evt;
	}

	fp->cd();
	sprintf(str, "meas_%s", name.c_str());
	h_meas->Write(str);
	sprintf(str, "true_%s", name.c_str());
	h_true->Write(str);

	delete fx;
	delete fxres;
	delete fxeff;

	return 0;
}

