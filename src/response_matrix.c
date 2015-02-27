#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

#include "math.h"

#include "response_matrix.h"

ResponseMatrix::ResponseMatrix(const char *name) :
	verbose(false), debug(false),
	ov(false), uf(false),
	h_x_true(0), h_x_meas(0), 
	h_x_y_true(0), h_x_y_meas(0), 
	h_x_y_z_true(0), h_x_y_z_meas(0),
	h_efficiency(0), h_dim(0), response(0), dimensions(0), weight(1.0) {
	int len = strlen(name);
	this->name = new char [ len + 1 ];
	sprintf(this->name, name);
	str = new char [ len + 256 ]; /* general purpose string to manipulate name */
}

bool ResponseMatrix::set_output_file(const char *file) {
	ofile = file;
	return true;
}

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

bool ResponseMatrix::initialize() {
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
	return true;
}

bool ResponseMatrix::set_weight(double weight) {
	this->weight = weight;
	return true;
}

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

