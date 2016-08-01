#ifndef __NRUTIL_
#define __NRUTIL_

#include <stdlib.h>

void nrerror(const char error_text[]);
float *vector(int nl,int nh);
int *ivector(int nl,int nh);
double *dvector(int nl,int nh);
float **matrix(int nrl,int nrh,int ncl,int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
double ***darray3(int n1l, int n1h, int n2l, int n2h, int n3l, int n3h);
double ****darray4(int n1l, int n1h, int n2l, int n2h, int n3l, int n3h, int n4l, int n4h);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void free_vector(float *v,int nl,int nh);
void free_ivector(int *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_darray3(float **m, int n1l, int n1h, int n2l, int n2h, int n3l, int n3h);
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);
float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);

#endif
