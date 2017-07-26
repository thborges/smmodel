
#ifndef ROUND_H
#define ROUND_H

#include "deps.h"

void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, int x[servers][opt_atu], 
	multiway_histogram_estimate *agg_server, double f, bool ignore_multiple,
	bool improve_exchange);

void set_cell_place_from_partial_x_nowh(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data);

void set_cell_place_from_partial_x(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data, int where[pairs]);

void remove_double_processed_items(int servers, optimization_data_s *opt_data, int opt_atu, 
	int x[servers][opt_atu], double f);

void schedule_non_assigned_items(int rem_count, int remaining[rem_count], int servers, 
	optimization_data_s *opt_data, double load[servers], int pairs, int where[pairs],
	int x[servers][pairs], double f);

void improve_transformed_solution(int servers, optimization_data_s *opt_data, int pairs, 
	double load[servers], int where[pairs], int x[servers][pairs], double f);

void improve_transformed_solution_fm(int servers, optimization_data_s *opt_data, int pairs, 
	double load[servers], int where[pairs], int x[servers][pairs], double f);

void improve_transformed_solution_exchange_pairs(int servers, optimization_data_s *opt_data, int pairs, 
	double load[servers], int where[pairs], int x[servers][pairs], double f);

#endif

