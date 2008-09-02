#include <stdlib.h>
void TDdXsi_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi);
void TDdXsi_surf_1(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_2(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_3(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_4(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_5(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_6(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_7(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_8(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi,double* z);
void TDdXsi_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dXsi,double* ddXsi) {
	double* z = (double*) malloc(11799*sizeof(double));
	TDdXsi_surf_1(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_2(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_3(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_4(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_5(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_6(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_7(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	TDdXsi_surf_8(params,Xsi,Xa,Ya,Za,Cmat,dXsi,ddXsi,z);
	free(z);
}
