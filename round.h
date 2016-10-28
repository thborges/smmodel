
#ifndef ROUND_H

#include <structs.h>

void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, int x[servers][opt_atu]);

#endif

