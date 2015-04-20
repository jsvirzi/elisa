#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <vector>
#include <string>

#include <dirent.h> 
#include <stdio.h>

/* jsv. read from file */
#define Intermediate 5

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
	std::string idir = "../experiment_2", ofile = "../posterior_2.root", name = "posterior_bin";
	struct dirent *dir;
	DIR *d = opendir(idir.c_str());
	std::vector<std::string> ifiles;
	char str[1024];
	double *y = new double [ NBINS ];
	int i, file_index, status;
	std::vector<std::string>::const_iterator ifile;
	TH1D **h = new TH1D * [ NBINS ];

	while((dir = readdir(d)) != NULL) {
		printf("%s\n", dir->d_name);
		std::string s = dir->d_name;
		if(s.find(".root") == std::string::npos) continue;
		ifiles.push_back(std::string(idir + "/" + s));
	}
	closedir(d);

/* create the histograms */
	for(i=0;i<NBINS;++i) {
		sprintf(str, "%s%d", name.c_str(), i);
		h[i] = new TH1D(str, "posterior", pdf_info[i].nbins, pdf_info[i].x_min, pdf_info[i].x_max);
		h[i]->SetDirectory(0);
	}

/* and fill them one ntuple at a time, instead of chaining them all together */
	for(ifile=ifiles.begin(),file_index=0;ifile!=ifiles.end();++ifile,++file_index) {
		printf("processing file %s\n", ifile->c_str());

	/* open the chain one file at a time and process */
		TChain chain("solution");
		chain.Add(ifile->c_str());
		chain.SetBranchAddress("y", y);
		chain.SetBranchAddress("status", &status);

	/* turn on only the branches we need */
		chain.SetBranchStatus("*", 0);
		chain.SetBranchStatus("y", 1);
		chain.SetBranchStatus("status", 1);

		int ientry, nentries = chain.GetEntries();
		for(ientry=0;ientry<nentries;++ientry) {
			chain.GetEvent(ientry);
			if(status != Intermediate) continue;
			if((ientry % 10000000) == 0) printf("processed %d/%d entries\n", ientry, nentries);
			// for(i=0;i<NBINS;++i) printf("y[%d] = %f\n", i, y[i]);
			// getchar();
			for(i=0;i<NBINS;++i) h[i]->Fill(y[i]);
		}
		printf("finished processing %s\n", ifile->c_str());
		TFile fp(ofile.c_str(), "update");
		for(i=0;i<NBINS;++i) {
			printf("BIN(%d) MEAN=%f RMS=%f\n", i, h[i]->GetMean(), h[i]->GetRMS());
			sprintf(str, "%s%d_file%d", name.c_str(), i, file_index);
			h[i]->Write(str);
		}
		fp.Write();
		fp.Close();
	}

	printf("finished processing all files\n");

	TFile fp(ofile.c_str(), "update");
	for(i=0;i<NBINS;++i) {
		sprintf(str, "%s%d", name.c_str(), i);
		h[i]->Write(str);
	}
	fp.Write();
	fp.Close();

	delete [] y;
	return 0;
}

