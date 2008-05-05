#include "nrutil.h"
using namespace std;

namespace dem {

#define NRANSI
#define MAXM 50

bool zrhqr(long double a[], int m, long double rtr[], long double rti[])
{
	void balanc(long double **a, int n);
	bool hqr(long double **a, int n, long double wr[], long double wi[]);
	int j,k;
	long double **hess,xr,xi;

	if (m > MAXM || a[m] == 0.0) // nrerror("bad args in zrhqr");
	    return false;

	hess=matrix(1,MAXM,1,MAXM);
	for (k=1;k<=m;k++) {
		hess[1][k] = -a[m-k]/a[m];
		for (j=2;j<=m;j++) hess[j][k]=0.0;
		if (k != m) hess[k+1][k]=1.0;
	}
	balanc(hess,m);
	if (!hqr(hess,m,rtr,rti)) {
	    free_matrix(hess,1,MAXM,1,MAXM);
	    return false;
	}

	for (j=2;j<=m;j++) {
		xr=rtr[j];
		xi=rti[j];
		for (k=j-1;k>=1;k--) {
			if (rtr[k] <= xr) break;
			rtr[k+1]=rtr[k];
			rti[k+1]=rti[k];
		}
		rtr[k+1]=xr;
		rti[k+1]=xi;
	}
	free_matrix(hess,1,MAXM,1,MAXM);
	return true;
}

#undef MAXM
#undef NRANSI

} // namespace dem ends
