
#ifndef ROUND_H

#include "structs.h"

void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, int x[servers][opt_atu], 
	multiway_histogram_estimate *agg_server, char remove_double_proc,
	double f, double g);

void set_cell_place_from_partial_x(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data);

#endif

