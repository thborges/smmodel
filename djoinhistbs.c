
#include <stdlib.h>
#include <glibwrap.h>
#include <limits.h>
#include "deps.h"
#include "round.h"
#include "sm.h"

double trade_off;

int choose_best_server_makespan(double *total_io, int servers, multiway_histogram_estimate *estimate) {
	// This heuristic choose the server with minor load to
	// process the cell. It works best when the
	// job/cell being scheduled is the one with major
	// load among all remaining cells to be scheduled

	// Choose the server with minor load
	int choosed = 1;
	int min_cost = estimate[1].to_pnts;
	for(int s = 2; s < servers; s++) {
		if (estimate[s].to_pnts < min_cost) {
			choosed = s;
			min_cost = estimate[s].to_pnts;
		}
	}

	return choosed;
}

int choose_best_server(double *total_io, int servers, multiway_histogram_estimate *estimate,
	multiway_histogram_estimate *agg_server) {

	// Choose the minor communication node
	int choosed = 1;
	double min_cost = total_io[1];
	for(int s = 2; s < servers; s++) {
		if (total_io[s] < min_cost) {
			choosed = s;
			min_cost = total_io[s];
		}
		else if (total_io[s] - min_cost <= 1e-9) { // same cost
			if (estimate[choosed].to_pnts > estimate[s].to_pnts)
				choosed = s;			
		}
	}

	// balance the Network
	/*int max_server = 1;
	int min_server = 1;
	double max_io_pnts = estimate[1].io_pnts;
	double min_io_pnts = estimate[1].io_pnts;
	for(int s = 2; s < servers; s++) {
		if (estimate[s].io_pnts > max_io_pnts) {
			max_server = s;
			max_io_pnts = estimate[s].io_pnts;
		}
		if (estimate[s].io_pnts < min_io_pnts) {
			min_server = s;
			min_io_pnts = estimate[s].io_pnts;
		}
	}
	double diff = max_io_pnts - min_io_pnts;
	if (diff / max_io_pnts > 0.20)
		choosed = min_server;
	*/

	// balance the CPU 
	int min_server = 1;
	double max_to_pnts = /*agg_server[1].to_pnts + */ estimate[1].to_pnts;
	double min_to_pnts = /*agg_server[1].to_pnts + */ estimate[1].to_pnts;
	for(int s = 2; s < servers; s++) {
		if (/*agg_server[s].to_pnts + */ estimate[s].to_pnts > max_to_pnts) {
			max_to_pnts = /*agg_server[s].to_pnts + */ estimate[s].to_pnts;
		}
		if (/*agg_server[s].to_pnts + */ estimate[s].to_pnts < min_to_pnts) {
			min_server = s;
			min_to_pnts = /*agg_server[s].to_pnts + */ estimate[s].to_pnts;
		}
	}
	double diff = max_to_pnts - min_to_pnts;
	if (diff / max_to_pnts > trade_off) {
		choosed = min_server;

		// find minor cost among unbalanced servers
		/*double min_cost = total_io[choosed];
		for(int s = 1; s < servers; s++) {
			if ((max_to_pnts - estimate[s].to_pnts) / max_to_pnts > trade_off) {
				if (min_cost > total_io[s]) {
					min_cost = total_io[s];
					choosed = s;
				}
			}
		}*/
	}

	return choosed;
}

int sort_optdata_decreasing(const void *x, const void *y) {
	optimization_data_s *iy = (optimization_data_s*)y;
	optimization_data_s *ix = (optimization_data_s*)x;
	if (ix->pnts == iy->pnts) return 0;
	if (ix->pnts < iy->pnts) return 1;
	return -1;
}

void bs_optimize_hr(dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
	int pairs, multiway_histogram_estimate *agg_server, double f) {

	//TODO: fix djoinhist to pass the correct value
	servers = servers - 1;

	int x[servers][pairs];
	memset(x, 0, sizeof x);

	int where[pairs];
	memset(where, -1, sizeof where);

	double load[servers];
	memset(load, 0, sizeof load);
	
	int unassigned_count = pairs;
	int unassigned[pairs];
	for(int i = 0; i < pairs; i++)
		unassigned[i] = i;

	schedule_non_assigned_items(unassigned_count, unassigned, servers, opt_data, load, pairs, where, x, f); 
	improve_transformed_solution(servers, opt_data, pairs, load, where, x, f);
	set_cell_place_from_partial_x(hr, servers, pairs, x, opt_data, where);

	double final_mkspan, final_comm;
	double Zheur = get_sm_objective(hr, opt_data, pairs, f, servers, 1, NULL, NULL, &final_mkspan, &final_comm);
	printf("After GR\nZ\tMkspan\tComm\nSM_GR %.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

}

void bs_optimize_hr_old(dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
	int opt_atu, multiway_histogram_estimate *agg_server) {

	char *toff = getenv("BS_TRADEOFF");
	trade_off = toff ? atof(toff) : 0.0;
	if (trade_off <= 0.0)
		trade_off = 0.2;
	//printf("Tradeoff: %f\n", trade_off);

	multiway_histogram_estimate estimate[servers];
	memset(estimate, 0, sizeof(multiway_histogram_estimate)*servers);

	// sort in decreasing order of points
	qsort(opt_data, opt_atu, sizeof(optimization_data_s), sort_optdata_decreasing);

	for(int c = 0; c < opt_atu; c++) {
		int xl = opt_data[c].xl;
		int yl = opt_data[c].yl;

		histogram_cell *resultcell = hr->get_cell(hr, xl, yl);

		// choose the best server for the L cell
		int choosed = choose_best_server(opt_data[c].comm, servers, estimate, agg_server);
		resultcell->place = choosed;
		SET_IN_PLACE(resultcell->copies, choosed);

		for(int i = 0; i < opt_data[c].rcells_size; i++) {
			histogram_cell *rc = opt_data[c].rcells[i].cell;
			if (!IS_IN_PLACE(rc->copies, choosed)) {
				SET_IN_PLACE(rc->copies, choosed);
				estimate[choosed].io_pnts += rc->points;
			}
		}

		histogram_cell *lc = opt_data[c].lcell;
		if (!IS_IN_PLACE(lc->copies, choosed)) {
			SET_IN_PLACE(lc->copies, choosed);
			estimate[choosed].io_pnts += lc->points;
		}

		estimate[choosed].to_pnts += opt_data[c].pnts;
	}

	double final_mkspan, final_comm;
	double Zheur = get_sm_objective(hr, opt_data, opt_atu, servers, servers, 1, NULL, NULL, &final_mkspan, &final_comm);
	printf("After BS\nZ\tMkspan\tComm\n%.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

}

