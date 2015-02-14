#ifndef UNFOLD_H
#define UNFOLD_H

class Unfold {

	public:

	Unfold(int N = 0);
	~Unfold();

/* the different options for unfolding */
	enum {
		UnfoldElisa,
		UnfoldFullBayesian,
		UnfoldMaximumLikelihood,
	};

	void set_seed(int seed) { this->seed = seed; };
	int get_seed() { return seed; };
	void set_dimension(int N) { this->N = N; };
	bool set_true(double *true_bins, int n = 0);
	bool set_meas(double *meas_bins, int n = 0);
	bool invert_matrix(double **M, int n = 0);
	bool initialize();
	bool fill(double xtrue, double xmeas);
	bool miss(double xtrue);
	bool save(const char *filename, int mode = 0);

	protected:

	int N;
	bool initialized;
	bool use_underflow, use_overflow;
	int seed; /* random number seed */
	double **R, **Rinv; /* the response matrix (and its inverse) */
	double *true_bins; /* the edges of the bins containing truth occupancies */
	double *meas_bins; /* the edges of the bins containing measured value occupancies */ 

};

#endif
