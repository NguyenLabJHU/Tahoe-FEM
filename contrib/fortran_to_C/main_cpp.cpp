/* $Id: main_cpp.cpp,v 1.4 2003-03-15 02:21:25 paklein Exp $ */
#include <iostream.h>
#include <math.h>
#include <iomanip.h>
#include "fortran_names.h"

extern "C" {
void FORTRAN_NAME(double_it)(double*, int*);
void FORTRAN_NAME(sqrt_f)(double*);
void FORTRAN_NAME(sqrt_f_c)(double*);
void FORTRAN_NAME(print_double)(double*);
void FORTRAN_NAME(print_int)(int*);
}

static bool d_check(double a_ref, double a_test) {

	double tol = 1.0e-14;

	if (fabs(a_ref) < tol)
		return fabs(a_ref - a_test) < tol;
	else
		return fabs(a_ref - a_test)/a_ref < tol;
};
static bool i_check(int a_ref, int a_test)
{
	return a_ref == a_test;
};

int main (int, char**)
{
	cout.precision(12);
	cout.setf(ios::showpoint);
	cout.setf(ios::right, ios::adjustfield);
	cout.setf(ios::scientific, ios::floatfield);
	int wd = 15;
	int wi = 5;

	double d1, d2;
	int i1, i2;
	
	d1 = d2 = cos(0.1);
	i1 = i2 = 94117;
	
	cout << "start: d = " << setw(wd) << d1 << ", i = " << setw(wi) << i1 << '\n';
	
	FORTRAN_NAME(double_it)(&d1, &i1);
	d2 *= 2.0;
	i2 *= 2;
	cout << "double an int and a double: ";
	if (d_check(d1,d2) && i_check(i1,i2))
		cout << "PASS" << endl;
	else
		cout << "FAIL" << endl;
	
	FORTRAN_NAME(sqrt_f)(&d1);
	d2 = sqrt(d2);
	cout << "sqrt of double: ";
	if (d_check(d1,d2))
		cout << "PASS" << endl;
	else
		cout << "FAIL" << endl;

	FORTRAN_NAME(sqrt_f_c)(&d1);
	d2 = sqrt(d2);
	cout << "sqrt of double in C called from Fortran: ";
	if (d_check(d1,d2))
		cout << "PASS" << endl;
	else
		cout << "FAIL" << endl;

	cout << "end: d = " << setw(wd) << d1 << ", i = " << setw(wi) << i1 << '\n';

	cout << "print values from Fortran functions:\n";
	FORTRAN_NAME(print_double)(&d1);
	FORTRAN_NAME(print_int)(&i1);
		
	cout << "Done." << endl;
	return 0;
}
