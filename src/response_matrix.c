#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

#include "math.h"

#include "response_matrix.h"
#include <stdio.h>

int find_bin(double a, int nbins, double *edges, bool uf = false, bool ov = false) {
	if(a < edges[0]) return uf ? 0 : -1; 
	if(a > edges[nbins]) return ov ? nbins - 1 : nbins;
	for(int bin=0;bin<nbins;++bin) {
		if(edges[bin] <= a && a < edges[bin+1]) {
			return bin;
		}
	}
	return -1;
}

ResponseMatrix::~ResponseMatrix() {
	if(eff) delete [] eff;
	if(eff_numer) delete [] eff_numer;
	if(eff_denom) delete [] eff_denom;
	if(true_bin_edges) delete [] true_bin_edges;
	if(meas_bin_edges) delete [] meas_bin_edges;
	if(R) {
		// for(int i=0;i<nt;++i) delete [] R[i];
		// delete [] R;
	}
}

bool ResponseMatrix::get_response_matrix(double **R) {
	for(int i=0;i<nt;++i) {
		for(int j=0;j<nr;++j) {
			R[i][j] = this->R[i][j];
		}
	}
	return true;
}

ResponseMatrix::ResponseMatrix(int nt, int nr) :
	verbose(false), debug(false),
	eff(0), eff_numer(0), eff_denom(0),
	nt(0), nr(0), true_dimensions(1), meas_dimensions(1),
//	ov(false), uf(false),
	true_bin_edges(0), meas_bin_edges(0), R(0) {
	// true_bin_edges(0), meas_bin_edges(0), R(0), weight(1.0) {
#if 0
	h_x_true(0), h_x_meas(0), 
	h_x_y_true(0), h_x_y_meas(0), 
	h_x_y_z_true(0), h_x_y_z_meas(0),
	h_efficiency(0), h_dim(0), dimensions(0)
#endif
//	this->name = name;
	R = new double * [ nt ];
	for(int i=0;i<nt;++i) R[i] = new double [ nr ];
#if 0
	int len = strlen(name);
	this->name = new char [ len + 1 ];
	sprintf(this->name, name);
	str = new char [ len + 256 ]; /* general purpose string to manipulate name */
#endif
}

#if 0
bool ResponseMatrix::set_output_file(const char *file) {
	ofile = file;
	return true;
}
#endif

bool ResponseMatrix::print() {
	printf("Response Matrix (%d X %d) = ", nt, nr);
	for(int i=0;i<nt;++i) {
		printf("\nROW(%d) = ", i);
		for(int j=0;j<nr;++j) printf("%f ", R[i][j]);
	}
	printf("\n"); 
	return true;
}

bool ResponseMatrix::set_true(int nt, double *bin_edges) {
	this->nt = nt;
	if(true_bin_edges) delete [] true_bin_edges;
	true_bin_edges = new double [ nt + 1 ];
	for(int i=0;i<=nt;++i) true_bin_edges[i] = bin_edges[i];
	initialize_true();
	return true;
}

bool ResponseMatrix::set_true(int nt, double a_min, double a_max) {
	this->nt = nt;
	if(true_bin_edges) delete [] true_bin_edges;
	true_bin_edges = new double [ nt + 1 ];
	double da = (a_max - a_min) / nt;
	for(int i=0;i<=nt;++i) true_bin_edges[i] = a_min + i * da; 
	initialize_true();
	return true;
}

bool ResponseMatrix::set_meas(int nr, double *bin_edges) {
	this->nr = nr;
	if(meas_bin_edges) delete [] meas_bin_edges;
	meas_bin_edges = new double [ nr + 1 ];
	for(int i=0;i<=nr;++i) meas_bin_edges[i] = bin_edges[i];
	initialize_meas();
	return true;
}

bool ResponseMatrix::set_meas(int nr, double a_min, double a_max) {
	this->nr = nr;
	if(meas_bin_edges) delete [] meas_bin_edges;
	meas_bin_edges = new double [ nr + 1 ];
	double da = (a_max - a_min) / nr;
	for(int i=0;i<=nr;++i) meas_bin_edges[i] = a_min + i * da; 
	initialize_meas();
	return true;
}

#if 0

bool ResponseMatrix::set_true(TH1D *h) {
	sprintf(str, "response_distribution_true_%s", name);
	h_x_true = new TH1D(*h);
	h_x_true->SetName(str);
	h_x_true->SetDirectory(0);
	h_x_true->Sumw2();
	return true;
}

bool ResponseMatrix::set_meas(TH1D *h) {
	sprintf(str, "response_distribution_meas_%s", name);
	h_x_meas = new TH1D(*h);
	h_x_meas->SetName(str);
	h_x_meas->SetDirectory(0);
	h_x_meas->Sumw2();
	return true;
}

bool ResponseMatrix::set_true(TH2D *h) {
	sprintf(str, "response_distribution_true_%s", name);
	h_x_y_true = new TH2D(*h);
	h_x_y_true->SetName(str);
	h_x_y_true->SetDirectory(0);
	h_x_y_true->Sumw2();
	return true;
}

bool ResponseMatrix::set_meas(TH2D *h) {
	sprintf(str, "response_distribution_meas_%s", name);
	h_x_y_meas = new TH2D(*h);
	h_x_y_meas->SetName(str);
	h_x_y_meas->SetDirectory(0);
	h_x_y_meas->Sumw2();
	return true;
}

bool ResponseMatrix::set_true(TH3D *h) {
	sprintf(str, "response_distribution_true_%s", name);
	h_x_y_z_true = new TH3D(*h);
	h_x_y_z_true->SetName(str);
	h_x_y_z_true->SetDirectory(0);
	h_x_y_z_true->Sumw2();
	return true;
}

bool ResponseMatrix::set_meas(TH3D *h) {
	sprintf(str, "response_distribution_meas_%s", name);
	h_x_y_z_meas = new TH3D(*h);
	h_x_y_z_meas->SetName(str);
	h_x_y_z_meas->SetDirectory(0);
	h_x_y_z_meas->Sumw2();
	return true;
}

#endif

bool ResponseMatrix::initialize_meas() {
	return true;
}

bool ResponseMatrix::initialize_true() {
#if 0
	int nbins_true = 0, nbins_meas = 0; 
	if(h_x_true && h_x_meas) {
		dimensions = 1;
		nbins_true = h_x_true->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0);
		nbins_meas = h_x_meas->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0);
	} else if(h_x_y_true && h_x_y_meas) {
		nbins_true = 
			(h_x_y_true->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_true->GetNbinsY() + (ov ? 1 : 0) + (uf ? 1 : 0));
		nbins_meas = 
			(h_x_y_meas->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_meas->GetNbinsY() + (ov ? 1 : 0) + (uf ? 1 : 0));
		dimensions = 2;
	} else if(h_x_y_z_true && h_x_y_z_meas) {
		nbins_true = 
			(h_x_y_z_true->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_z_true->GetNbinsY() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_z_true->GetNbinsZ() + (ov ? 1 : 0) + (uf ? 1 : 0));
		nbins_meas = 
			(h_x_y_z_meas->GetNbinsX() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_z_meas->GetNbinsY() + (ov ? 1 : 0) + (uf ? 1 : 0)) * 
			(h_x_y_z_meas->GetNbinsZ() + (ov ? 1 : 0) + (uf ? 1 : 0));
		dimensions = 3;
	}
	sprintf(str, "dimensions_%s", name);
	h_dim = new TH2D(str, "dimensions", dimensions, -0.5, dimensions - 0.5, dimensions, -0.5, dimensions - 0.5);
	sprintf(str, "response_data_%s", name);
	response = new TH2D(str, "true vs meas", 
		nbins_meas, 0.5, nbins_meas + 0.5, 
		nbins_true, 0.5, nbins_true + 0.5);  
	response->SetDirectory(0);
	response->Sumw2();
	sprintf(str, "efficiency_%s", name);
	h_efficiency = new TH1D(str, "efficiency", nbins_true, 0.5, nbins_true + 0.5);
	h_efficiency->SetDirectory(0);
	h_efficiency->Sumw2();
	sprintf(str, "efficiency_denom_%s", name);
	h_efficiency_denom = new TH1D(str, "efficiency denom", nbins_true, 0.5, nbins_true + 0.5);
	h_efficiency_denom->SetDirectory(0);
	h_efficiency_denom->Sumw2();
	sprintf(str, "efficiency_numer_%s", name);
	h_efficiency_numer = new TH1D(str, "efficiency numer", nbins_true, 0.5, nbins_true + 0.5);
	h_efficiency_numer->SetDirectory(0);
	h_efficiency_numer->Sumw2();
#endif
	if(eff) delete [] eff;
	eff = new double [ nt ];
	if(eff_numer) delete [] eff_numer;
	eff_numer = new double [ nt ];
	if(eff_denom) delete [] eff_denom;
	eff_denom = new double [ nt ];
	return true;
}

#if 0
bool ResponseMatrix::set_weight(double weight) {
	this->weight = weight;
	return true;
}
#endif

bool ResponseMatrix::hit(double x_true, double x_meas, double weight) {
	int bin_true = find_bin(x_true, nt, true_bin_edges, uf, ov);
	int bin_meas = find_bin(x_meas, nr, meas_bin_edges, uf, ov);
	if(bin_true < 0 || bin_meas < 0) return false;
	if(bin_true >= nt || bin_meas >= nr) return false;
	R[bin_true][bin_meas] += weight;
	eff_denom[bin_true] += weight;
	eff_numer[bin_true] += weight;
	return true;
}

bool ResponseMatrix::miss(double x_true, double weight) {
	int bin_true = find_bin(x_true, nt, true_bin_edges, uf, ov);
	if(bin_true < 0) return false;
	eff_denom[bin_true] += weight;
	return true;
}

#if 0

bool ResponseMatrix::hit(double x_true, double x_meas) {
	int bin_true = h_x_true->FindBin(x_true);
	int bin_meas = h_x_meas->FindBin(x_meas);
	response->Fill(bin_meas, bin_true, weight);
	h_x_true->Fill(x_true, weight);
	h_x_meas->Fill(x_meas, weight);
	h_efficiency_numer->Fill(bin_true, weight);
	h_efficiency_denom->Fill(bin_true, weight);
	return true;
}

bool ResponseMatrix::hit(double x_true, double y_true, double x_meas, double y_meas) {
	int bin_true = h_x_y_true->FindBin(x_true, y_true);
	int bin_meas = h_x_y_meas->FindBin(x_meas, y_meas);
	response->Fill(bin_meas, bin_true, weight);
	h_efficiency_numer->Fill(bin_true, weight);
	h_efficiency_denom->Fill(bin_true, weight);
	h_x_y_true->Fill(x_true, y_true, weight);
	h_x_y_meas->Fill(x_meas, y_meas, weight);
	return true;
}

bool ResponseMatrix::hit(double x_true, double y_true, double z_true, double x_meas, double y_meas, double z_meas) {
	int bin_true = h_x_y_z_true->FindBin(x_true, y_true, z_true);
	int bin_meas = h_x_y_z_meas->FindBin(x_meas, y_meas, z_meas);
	response->Fill(bin_meas, bin_true, weight);
	h_efficiency_numer->Fill(bin_true, weight);
	h_efficiency_denom->Fill(bin_true, weight);
	h_x_y_z_true->Fill(x_true, y_true, z_true, weight);
	h_x_y_z_meas->Fill(x_meas, y_meas, z_meas, weight);
	return true;
}

bool ResponseMatrix::miss(double x_true) {
	h_x_true->Fill(x_true, weight);
	int bin_true = h_x_true->FindBin(x_true);
	h_efficiency_denom->Fill(bin_true, weight);
	return true;
}

bool ResponseMatrix::miss(double x_true, double y_true) {
	h_x_y_true->Fill(x_true, y_true, weight);
	int bin_true = h_x_y_true->FindBin(x_true, y_true);
	h_efficiency_denom->Fill(bin_true, weight);
	return true;
}

bool ResponseMatrix::miss(double x_true, double y_true, double z_true) {
	h_x_y_z_true->Fill(x_true, y_true, z_true, weight);
	int bin_true = h_x_y_z_true->FindBin(x_true, y_true);
	h_efficiency_denom->Fill(bin_true, weight);
	return true;
}

ResponseMatrix::~ResponseMatrix() {
	if(response) delete response;
	if(h_x_true) delete h_x_true;
	if(h_x_meas) delete h_x_meas;
	if(h_x_y_true) delete h_x_y_true;
	if(h_x_y_meas) delete h_x_y_meas;
	if(h_x_y_z_true) delete h_x_y_z_true;
	if(h_x_y_z_meas) delete h_x_y_z_meas;
	if(h_efficiency) delete h_efficiency;
	if(h_efficiency_numer) delete h_efficiency_numer;
	if(h_efficiency_denom) delete h_efficiency_denom;
	if(h_dim) delete h_dim;
}

#endif

bool ResponseMatrix::write(std::string &ofile) {
	for(int bin=0;bin<nt;++bin) {
		double ad = eff_denom[bin];
		double an = eff_numer[bin];
		if(ad != 0.0) { eff[bin] = an / ad; }
		else { eff[bin] = 0.0; }
	}
	if(ofile.length()) {
		FILE *fp = fopen(ofile.c_str(), "w");
		fwrite(&nt, 1, sizeof(int), fp);
		fwrite(&nr, 1, sizeof(int), fp);
		fwrite(true_bin_edges, nt+1, sizeof(double), fp);
		fwrite(meas_bin_edges, nr+1, sizeof(double), fp);
	/* write out the raw (unnormalized) response matrix */
		for(int i=0;i<nt;++i) fwrite(R[i], nr, sizeof(double), fp);

	/* truth distribution */
		double *x = new double [ nt ];
		for(int i=0;i<nt;++i) {
			double acc = 0.0;
			for(int j=0;j<nr;++j) acc += R[i][j];
			x[i] = acc;
		}
		fwrite(x, nt, sizeof(double), fp);
		delete [] x;

	/* measured distribution */
		double *y = new double [ nr ];
		for(int j=0;j<nr;++j) {
			double acc = 0.0;
			for(int i=0;i<nt;++i) acc += R[i][j];
			y[j] = acc;
		}
		fwrite(y, nr, sizeof(double), fp);
		delete [] y;

	/* write out the final response matrix */
		for(int i=0;i<nt;++i) {
			double acc = 0.0;
			for(int j=0;j<nr;++j) acc += R[i][j];
			for(int j=0;j<nr;++j) R[i][j] = R[i][j] * eff[i] / acc;
			fwrite(R[i], nr, sizeof(double), fp);
			printf("R[%d]: (eff = %f)\n", i, eff[i]);
			for(int j=0;j<nr;++j) printf(" %f", R[i][j]);
		}
		fclose(fp);
#if 0
		fprintf(fp, "<response_matrix>\n");
		fprintf(fp, "\t<true_dimensions>%d</true_dimensions>\n", true_dimensions);
		fprintf(fp, "\t<meas_dimensions>%d</meas_dimensions>\n", meas_dimensions);
		fprintf(fp, "\t<true_bins>%d</true_bins>\n", nt);
		for(int bin=0;bin<nt;++bin) {
			fprintf(fp, "\t\t<true_bin_edge>%f</true_bin_edge>\n", true_bin_edges[i]);
		}
		fprintf(fp, "\t<meas_bins>%d</meas_bins>\n", nt);
		for(int bin=0;bin<nt;++bin) {
			fprintf(fp, "\t\t<meas_bin_edge>%f</meas_bin_edge>\n", meas_bin_edges[i]);
		}

		fprintf(fp, "</response_matrix>\n"
#endif
	}
}

ResponseMatrix::ResponseMatrix(const char *file) {
	int i;
	FILE *fp = fopen(file, "r");
	printf("opening file [%s]\n", file);
	if(fp == NULL) { printf("blah\n"); return; }
	fread(&nt, 1, sizeof(int), fp);
	fread(&nr, 1, sizeof(int), fp);
printf("nt = %d. nr = %d. sizeof(int) = %d\n", nt, nr, sizeof(int));
	R = new double * [ nt ];
	for(i=0;i<nt;++i) R[i] = new double [ nr ];
	true_bin_edges = new double [ nt ];
	meas_bin_edges = new double [ nr ];
	fread(true_bin_edges, nt+1, sizeof(double), fp);
	fread(meas_bin_edges, nr+1, sizeof(double), fp);
	/* read in the raw (unnormalized) response matrix. not useful for the moment */
	for(i=0;i<nt;++i) fread(R[i], nr, sizeof(double), fp);
	double *x = new double [ nt ];
	fread(x, nt, sizeof(double), fp);
	delete [] x;
	double *y = new double [ nr ];
	fread(y, nr, sizeof(double), fp);
	delete [] y;
	for(i=0;i<nt;++i) fread(R[i], nr, sizeof(double), fp);
	fclose(fp);
}

#if 0

bool ResponseMatrix::finalize(std::string &ofile) {
	int nbins_true = h_efficiency->GetNbinsX();
	for(int bin=0;bin<=nbins_true;++bin) {
		double ad = h_efficiency_denom->GetBinContent(bin);
		double an = h_efficiency_numer->GetBinContent(bin);
		double ed = h_efficiency_denom->GetBinError(bin);
		double en = h_efficiency_numer->GetBinError(bin);
		if(ad != 0.0 && an != 0.0) { /* more or less normal */ 
			double eff = an / ad;
			double deff = eff * sqrt(pow(en/an,2.0)+pow(ed/ad,2.0));
			h_efficiency->SetBinContent(bin, eff);
			h_efficiency->SetBinError(bin, deff);
		} else {
			h_efficiency->SetBinContent(bin, 0.0);
			h_efficiency->SetBinError(bin, 0.0);
		}
	}
	if(ofile.length()) {
		fp = new TFile(ofile.c_str(), "update");
		if(response) response->Write();
		if(h_x_true) h_x_true->Write();
		if(h_x_meas) h_x_meas->Write();
		if(h_x_y_true) h_x_y_true->Write();
		if(h_x_y_meas) h_x_y_meas->Write();
		if(h_x_y_z_true) h_x_y_z_true->Write();
		if(h_x_y_z_meas) h_x_y_z_meas->Write();
		if(h_efficiency) h_efficiency->Write();
		if(h_efficiency_numer) h_efficiency_numer->Write();
		if(h_efficiency_denom) h_efficiency_denom->Write();
		if(h_dim) h_dim->Write();

	/* process the response matrix. columns (truth) are normalized to efficiency */
		int nbinsx = response->GetNbinsX(), nbinsy = response->GetNbinsY();
		for(int biny=0;biny<=(nbinsy+1);++biny) {
			if(!uf && (biny == 0)) continue;
			if(!ov && (biny == (nbinsy+1))) continue;
			double acc = 0.0, eff = h_efficiency->GetBinContent(biny);
			for(int binx=0;binx<=(nbinsx+1);++binx) {
				if(!uf && (binx == 0)) continue;
				if(!ov && (binx == (nbinsx+1))) continue;
				acc += response->GetBinContent(binx, biny);
			}
			for(int binx=0;binx<=(nbinsx+1);++binx) {
				if(!uf && (binx == 0)) continue;
				if(!ov && (binx == (nbinsx+1))) continue;
				double a = response->GetBinContent(binx, biny);
				a = a * eff / acc;
				response->SetBinContent(binx, biny, a);
			}
		}
		sprintf(str, "response_%s", name);
		response->Write(str);
		fp->Write();
		delete fp;
	}
	return true;
}

#endif
