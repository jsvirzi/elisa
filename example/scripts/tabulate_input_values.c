#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"

const int InputValue = 1;
const int CentralValue = 0;

void extract_central_value(const char *file, double *y, int nbins) {
	TBranch *b_status = 0, *b_y = 0, *b_n = 0;
	int status, event, nevents;
	double *ytemp = new double [ nbins + 2 ];
	double *ntemp = new double [ nbins + 2 ];
	TChain chain("extraction");
	chain.AddFile(file);
	b_status = b_y = b_n = 0;
	chain.SetBranchAddress("status", &status, &b_status);
	chain.SetBranchAddress("y", ytemp, &b_y);
	chain.SetBranchAddress("n", ntemp, &b_n);
	nevents = chain.GetEntries();
	printf("processing bayes with %d entries\n", nevents);
	for(event=0;event<nevents;++event) {
		b_status->GetEntry(event);
		b_y->GetEntry(event);
		b_n->GetEntry(event);
		if(status != CentralValue) continue;
		else break;
	}
	for(int i=0;i<nbins;++i) y[i] = ytemp[i];
printf("getting up\n");
	delete [] ytemp;
	delete [] ntemp;
printf("going home\n");
}

int tabulate_input_values() {
	char str[1024];
	TFile *fp;

	fp = new TFile("../measurement_1.root", "read");
	TH1D *h_meas_1 = (TH1D *)fp->Get("meas_x");
	TH1D *h_true_1 = (TH1D *)fp->Get("true_x");
	h_meas_1->SetDirectory(0);
	h_true_1->SetDirectory(0);
	fp->Close();
	delete fp;

#if 0
	fp = new TFile("../measurement_2.root", "read");
	TH1D *h_meas_2 = (TH1D *)fp->Get("meas_x");
	TH1D *h_true_2 = (TH1D *)fp->Get("true_x");
	h_meas_2->SetDirectory(0);
	h_true_2->SetDirectory(0);
	fp->Close();
	delete fp;
#endif

	int i, nbins = h_true_1->GetNbinsX();
	printf("%d bins detected\n", nbins);

#if 0
	fp = new TFile("../posterior_1.root", "read");
	double *yelisa = new double [ nbins ];
	for(i=0;i<nbins;++i) {
		sprintf(str, "posterior_bin%d", i);
		TH1D *h = (TH1D *)fp->Get(str);
		yelisa[i] = h->GetMean();
	}
	fp->Close();
	delete fp;
#endif

	TChain *chain;
	TBranch *b_status = 0, *b_y = 0, *b_n = 0;
	int status, event, nevents;

	double *ybayes = new double [ nbins ];
	double *nbayes = new double [ nbins ];
	double *yelisa = new double [ nbins ];
	double *nelisa = new double [ nbins ];

	extract_central_value("../bayes_analysis.root", ybayes, nbins);
	extract_central_value("../elisa_analysis.root", yelisa, nbins);

	printf("processing final output\n");
	for(i=0;i<nbins;++i) {
		double true_1 = h_true_1->GetBinContent(i+1);
		// double true_2 = h_true_2->GetBinContent(i+1);
		printf("%10.1f & %10.1f & %10.1f & %10.1f & %.4f & %10.1f & %10.1f & %.4f\n",
			h_true_1->GetBinLowEdge(i+1), h_true_1->GetBinLowEdge(i+2),
			true_1, h_meas_1->GetBinContent(i+1), yelisa[i] / true_1,
		//	true_2, h_meas_2->GetBinContent(i+1),
			ybayes[i], nbayes[i], ybayes[i] / true_1
			);
	}

	return 0;

	delete [] yelisa;
	delete [] nelisa;
	delete [] ybayes;
	delete [] nbayes;
	delete h_meas_1;
	delete h_true_1;
	// delete h_meas_2;
	// delete h_true_2;

	exit(0);
}
