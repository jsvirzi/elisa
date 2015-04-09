#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"

#define NBINS 12

struct info {
	int bin, npoints;
	double xmin, xmax;
};

struct info info1[NBINS+1] = {
	{ 0, 100, 1764000, 1776000 },
	{ 1, 100, 1430000, 1440000 }, 
	{ 2, 100, 547000, 554000 }, 
	{ 3, 100, 173000, 177000 }, 
	{ 4, 100, 49500, 51400 }, 
	{ 5, 100, 13400, 14600 }, 
	{ 6, 100, 3400, 4050 }, 
	{ 7, 100, 650, 1000 }, 
	{ 8, 100, 170, 360 }, 
	{ 9, 100, 5, 80 }, 
	{ 10, 100, 0, 30 }, 
	{ 11, 100, 0, 15 }, 
	{ 0, 0, 0.0, 0.0 }
};

struct info info2[NBINS+1] = {
	{ 0, 100, 1761000, 1773000 }, 
	{ 1, 100, 1427000, 1439000 }, 
	{ 2, 100, 548000, 555000 }, 
	{ 3, 100, 174000, 178000 }, 
	{ 4, 100, 50000, 52200 }, 
	{ 5, 100, 13000, 14200 }, 
	{ 6, 100, 3150, 3850 }, 
	{ 7, 100, 750, 1150 }, 
	{ 8, 100, 120, 320 }, 
	{ 9, 100, 10, 110 }, 
	{ 10, 100, 0, 30 }, 
	{ 11, 100, 0, 15 },
	{ 0, 0, 0.0, 0.0 }
};

int make_prior(int variant) {
	char str[1024];
	TChain *solution = new TChain("solution");
	sprintf(str, "../elisa_prior_from_%d.root", variant);
	solution->AddFile(str);
	int bin, nbins = NBINS;

	for(bin=1;bin<=nbins;++bin) {
		struct info *p = (variant == 1) ? &info1[bin-1] : &info2[bin-1];
		sprintf(str, "elisa_bin%d", p->bin);
		TH1D themp(str, str, p->npoints, p->xmin, p->xmax);
		sprintf(str, "y[%d] >>elisa_bin%d", p->bin, p->bin);
		solution->Draw(str, "status == 1");
		int i, n = themp.GetNbinsX();
		for(i=1;i<=n;++i) {
			double a = themp.GetBinContent(i);
			a = a * themp.GetBinCenter(i);
			themp.SetBinContent(i, a);
		}
		sprintf(str, "prior%d.root", variant);
		TFile fp(str, "update");
		themp.Write();
		fp.Close();
	}

	return 0;
}

int make_prior() {
	make_prior(1);
	make_prior(2);
}

