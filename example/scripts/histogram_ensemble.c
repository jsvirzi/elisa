#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

bool debug = false, verbose = false;

int main(int argc, char **argv) {
	bool update = false;
	std::string ifiles, ofile, oname;
	int nbins = 0, bin = -1;
	double xmin = 0.0, xmax = 0.0;
	for(int i=1;i<argc;++i) {
		if(strcmp("-debug", argv[i]) == 0) { debug = true;
		} else if(strcmp("-verbose", argv[i]) == 0) { verbose = true;
		} else if(strcmp("-update", argv[i]) == 0) { update = true;
		} else if(strcmp("-bin", argv[i]) == 0) { bin = atoi(argv[++i]); 
		} else if(strcmp("-i", argv[i]) == 0) { ifiles = argv[++i];
		} else if(strcmp("-o", argv[i]) == 0) { ofile = argv[++i]; oname = argv[++i];
		} else if(strcmp("-hist", argv[i]) == 0) { 
			nbins = atoi(argv[++i]); 
			xmin = atof(argv[++i]);
			xmax = atof(argv[++i]);
		} else { printf("unrecognized argument [%s]\n", argv[i]); exit(1); 
		}
	}

	char str[1024];

	printf("input spec = [%s]\n", ifiles.c_str());

	TChain *chain = new TChain("extraction");
	chain->Add(ifiles.c_str());
	sprintf(str, "bin%d", bin);
	TH1D *h = new TH1D(oname.c_str(), str, nbins, xmin, xmax);
	sprintf(str, "y[%d] >>%s", bin, oname.c_str());
	printf("draw command = [%s]\n", str);
	chain->Draw(str, "status == 2");

	TFile fp(ofile.c_str(), update ? "update" : "recreate");
	h->Write();
	fp.Close();

	delete chain;

	return 0;
}
