/* $Id: util.c,v 1.2 2004-12-11 09:27:23 paklein Exp $ */
#include <stdio.h>

/* zero the array */
void psp_clear(double* a, int n) {
	int i;
	for (i = 0; i < n; i++)
		a[i] = 0.0;
};

void psp_clear_int(int* a, int n) {
	int i;
	for (i = 0; i < n; i++)
		a[i] = 0;
};

/* write array to file stream */
void psp_dump(FILE* fp, double* a, int n) {
	const int wrap = 5;
	int count = 0;
	int i;
	for (i = 0; i < n; i++) {
		fprintf(fp, "%12.4g", a[i]);
		if ((++count) % wrap == 0) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	if (count != 0)
		fprintf(fp, "\n");
};
