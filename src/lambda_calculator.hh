// Copyright 2008 Michiaki Hamada

#ifndef __H_INCLUDE_LAMBDA_CALCULATOR_HH
#define __H_INCLUDE_LAMBDA_CALCULATOR_HH

#include "nrutil.hh"
#include "lubksb.hh"

// These pointers are 1-based!
double calculate_lambda( const double** mat_b, int alpha_size, double* p, double* q );

void makematrix( const double** mat_b, double **a, double lambda);

typedef struct Lambda {
  double min;
    double max;
	  int flag;    // 1 means there is a range, -1 means no solution possible.
} Lambda;

typedef struct Sum {
  double value;
  int flag;   // 1 means no negative bg_freq, -1 means there is negative bg_freq
} Sum;

Lambda Find_JP(const double** mat_b, double la_min, double la_max, double **JP, double *p_in, double *q_in);

Sum Check_root(const double** mat_b, double **a, double lambda, double *p, double *q);

double Check_det(const double** mat_b, double **a, double lambda);

Sum Nail_lambda(const double** mat_b, int flag_sign, double lambda_min, double lambda_max, double *p, double *q, double *la_add);

double Nail_det(const double** mat_b, int flag_sign, double lambda_min, double lambda_max);

 bool Check_range( const double** mat_b );

double *Locate_det_zero( const double** mat_b, int * ); //pointer to root list are returned with how many of them by int 


#endif // __H_INCLUDE_LAMBDA_CALCULATOR_HH
