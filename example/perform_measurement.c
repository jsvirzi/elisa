#include "TF1.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"

#include "stdlib.h"

#include "response_matrix.h"

bool debug = false, verbose = false;

double true_integral(double a, double x_lo, double x_hi) {
	double t_lo = (1.0 + a * x_lo) * TMath::Exp(-a * x_lo);
	double t_hi = (1.0 + a * x_hi) * TMath::Exp(-a * x_hi);
	return (t_lo - t_hi) / (a * a);
}

int main(int argc, char **argv) {
	int i, j, nevts = 1000000, seed = 5743;
	char str[1024];

/* X distribution parameters */
	double x_min = 0.0, x_max = 2.2, x_a = 7.5;

	std::string ofile, rfile = "example/response_matrix.dat", name = "x";

	for(i=1;i<argc;++i) {
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

/* X */
	int nbins_x = 11;

	double x_max_generation = 3.6; /* matches the parameter in response matrix generation */
	TH1D *h_true = new TH1D("true_x", "X", nbins_x, x_min, x_max); 
	TH1D *h_meas = new TH1D("meas_x", "X", nbins_x, x_min, x_max); 

	sprintf(str, "x * TMath::Exp(-%f * x)", x_a);
	TF1 *fx = new TF1("fx", str, x_min, x_max_generation); 

	for(int evt=0;evt<nevts;) { /* evt incremented after successful generation */
		if(evt && ((evt % 100000) == 0)) printf("%d / %d events processed\n", evt, nevts);
		double true_x = fx->GetRandom(); /* random x on truth interval */
		h_true->Fill(true_x); /* fill truth unconditionally */
		++evt;
	}

	ResponseMatrix *response = new ResponseMatrix(rfile.c_str());
	int nt = response->get_nt(), nr = response->get_nr(); 
	double **R = new double * [ nt ]; 
	for(i=0;i<nt;++i) R[i] = new double [ nr ];
	response->get_response_matrix(R);
	response->print();

	double *y0 = new double [ nt ];
	for(i=0;i<nt;++i) y0[i] = h_true->GetBinContent(1 + i); 
	for(i=0;i<nt;++i) printf("TRUE(%d) = %f\n", i, y0[i]); 
	for(j=0;j<nr;++j) {
		double mu = 0.0;
		for(i=0;i<nt;++i) mu += y0[i] * R[i][j]; 
		printf("<MU(%d)> = %f => ", j, mu);
		mu = rndm->Poisson(mu);
		printf("%f\n", mu);
		h_meas->SetBinContent(j+1, mu);
	}

	TFile *fp = new TFile(ofile.c_str(), "update");
	sprintf(str, "meas_%s", name.c_str());
	h_meas->Write(str);
	sprintf(str, "true_%s", name.c_str());
	h_true->Write(str);
	fp->Write();
	fp->Close();
	delete fp;

	delete fx;

	return 0;
}

