#ifndef RESPONSEMATRIX_H
#define RESPONSEMATRIX_H

#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

class ResponseMatrix {
	public:
	ResponseMatrix(const char *name = "response_matrix");
	~ResponseMatrix();
	bool hit(double x_true, double x_meas);
	bool hit(double x_true, double y_true, double x_meas, double y_meas);
	bool hit(double x_true, double y_true, double z_true, double x_meas, double y_meas, double z_meas);
	bool miss(double x_true);
	bool miss(double x_true, double y_true);
	bool miss(double x_true, double y_true, double z_true);
	bool set_true(TH1D *h);
	bool set_meas(TH1D *h);
	bool set_true(TH2D *h);
	bool set_meas(TH2D *h);
	bool set_true(TH3D *h);
	bool set_meas(TH3D *h);
	bool initialize();
	bool finalize(std::string &ofile);
	bool set_weight(double weight);
	bool set_output_file(const char *file);
	bool enable_overflow(bool flag = true) { ov = flag; return true; };
	bool enable_underflow(bool flag = true) { uf = flag; return true; };
	protected:
	bool verbose, debug, ov, uf;
	TH1D *h_x_true, *h_x_meas;
	TH2D *h_x_y_true, *h_x_y_meas;
	TH3D *h_x_y_z_true, *h_x_y_z_meas;
	TH1D *h_efficiency, *h_efficiency_denom, *h_efficiency_numer;
	TH2D *h_dim;
	TH2D *response;
	int dimensions;
	char *name, *str;
	double weight; // , x_true_lo, x_true_hi, x_meas_lo, x_meas_hi;
	std::string ofile;
	TFile *fp;
};

#endif
