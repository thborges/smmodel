#ifndef SM_H
#define SM_H

double get_sm_objective(dataset_histogram *hr, optimization_data_s *opt_data, int opt_atu,
	double f, int servers, double multiplier, double *x0_ub, double *x0_lb, double *fmkspan, double *fcomm);



#endif

