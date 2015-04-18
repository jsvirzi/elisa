#include "TH1D.h"
#include "TFile.h"

#include <vector>
#include <string>

#include <dirent.h> 
#include <stdio.h>

#define NBINS 12

int main(int argc, char **argv) {
	std::string idir = "../pdf_1_prior_2", ofile = "../pdf_1_prior_2.root", name = "pdf_bin";
	struct dirent *dir;
	DIR *d = opendir(idir.c_str());
	std::vector<std::string> ifiles;
	while((dir = readdir(d)) != NULL) {
		printf("%s\n", dir->d_name);
		std::string s = dir->d_name;
		if(s.find(".root") == std::string::npos) continue;
		ifiles.push_back(std::string(idir + "/" + s));
		
	}
	closedir(d);

	int bin;
	char str[1024];
	for(bin=0;bin<NBINS;++bin) {
		sprintf(str, "%s%d", name.c_str(), bin);
		std::vector<std::string>::const_iterator ifile = ifiles.begin();
		TFile fp(ifile->c_str(), "read");
		TH1D *h0 = (TH1D *)fp.Get(str);
		for(++ifile;ifile!=ifiles.end();++ifile) {
			TFile fp(ifile->c_str(), "read");
			TH1D *h = (TH1D *)fp.Get(str);
			h0->Add(h);
		}
		TFile fpo(ofile.c_str(), "update");
		h0->Write();
		fpo.Write();
		fpo.Close();
	}
}

