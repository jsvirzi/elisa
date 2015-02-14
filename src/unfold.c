#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

#include <vector>

#include "math.h"

#include "unfold.h"

Unfold::Unfold(int N) : initialized(false), use_underflow(false), use_overflow(false) {
	int i;
	this->N = N;
	if(N) {
		R = new double * [ N + 2 ];
		for(i=0;i<(N+2);++i) R[i] = new double [ N + 2 ];
		Rinv = new double * [ N + 2 ];
		for(i=0;i<(N+2);++i) Rinv[i] = new double [ N + 2 ];
	}
}

Unfold::~Unfold() {
	int i;
	if(R) {
		for(i=0;i<(N+2);++i) delete [] R[i];
		delete [] R;
	}

	if(Rinv) {
		for(i=0;i<(N+2);++i) delete [] Rinv[i];
		delete [] Rinv;
	}
}

bool Unfold::fill(double xtrue, double xmeas) {
	return true;
}

bool Unfold::miss(double xtrue) {
	return true;
}
