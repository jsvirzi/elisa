#ifndef RESPONSEMATRIX_H
#define RESPONSEMATRIX_H

#include <vector>

#if 0
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#endif

class ResponseMatrix {
	public:
	ResponseMatrix(int nt, int nr);
	ResponseMatrix(const char *file);
	~ResponseMatrix();
	bool hit(double x_true, double x_meas, double weight);
	bool hit(double x_true, double y_true, double x_meas, double y_meas, double weight);
	bool hit(double x_true, double y_true, double z_true, double x_meas, double y_meas, double z_meas, double weight);
	bool miss(double x_true, double weight);
	bool miss(double x_true, double y_true, double weight);
	bool miss(double x_true, double y_true, double z_true, double weight);
	bool set_true(int nt, double *bin_edges);
	bool set_true(int nt, double a_min, double a_max);
	bool set_meas(int nr, double *bin_edges);
	bool set_meas(int nt, double a_min, double a_max);
	bool get_response_matrix(double **R);
	int get_nt() { return nt; }
	int get_nr() { return nr; }
	bool print();
#if 0
	bool set_true(TH1D *h);
	bool set_meas(TH1D *h);
	bool set_true(TH2D *h);
	bool set_meas(TH2D *h);
	bool set_true(TH3D *h);
	bool set_meas(TH3D *h);
#endif
	bool initialize_true(), initialize_meas();
	bool write(std::string &ofile);
	// bool set_weight(double weight);
	// bool set_output_file(const char *file);
	bool enable_overflow(bool flag = true) { ov = flag; return true; };
	bool enable_underflow(bool flag = true) { uf = flag; return true; };
	protected:
	bool verbose, debug;
	bool ov, uf;
	int nt, nr, true_dimensions, meas_dimensions;
#if 0
	TH1D *h_x_true, *h_x_meas;
	TH2D *h_x_y_true, *h_x_y_meas;
	TH3D *h_x_y_z_true, *h_x_y_z_meas;
	TH1D *h_efficiency, *h_efficiency_denom, *h_efficiency_numer;
	TH2D *h_dim;
	TH2D *response;
	int dimensions;
	TFile *fp;
#endif
	double **R, *true_bin_edges, *meas_bin_edges, *eff, *eff_numer, *eff_denom;
	// char *name; // , *str;
	// double weight; // , x_true_lo, x_true_hi, x_meas_lo, x_meas_hi;
	std::string name, ofile;
};

#endif
