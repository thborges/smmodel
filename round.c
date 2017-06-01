
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <multiway-join.h>
#include "structs.h"
#include <float.h>
//#include "system.h"

typedef struct {
	double f;
	double g;
	int servers;
	optimization_data_s* opt_data;
} thunk_sort_lagrange;

decl_qsort_p_cmp(sort_remaining_decreasing_larger_cost, x, y, thunk) {
	int iy = *(int*)y;
	int ix = *(int*)x;
	thunk_sort_lagrange *thunkd = (thunk_sort_lagrange*)thunk;

	double comm_ix2 = MAX(thunkd->opt_data[ix].comm[1], thunkd->opt_data[ix].comm[2]);
	double comm_ix = MIN(thunkd->opt_data[ix].comm[1], thunkd->opt_data[ix].comm[2]);
	for(int s = 1; s <= thunkd->servers; s++) {
		if (comm_ix < thunkd->opt_data[ix].comm[s]) {
			comm_ix2 = comm_ix;
			comm_ix = thunkd->opt_data[ix].comm[s];
		}
	}
	comm_ix = comm_ix2 - comm_ix;

	double comm_iy2 = MAX(thunkd->opt_data[iy].comm[1], thunkd->opt_data[iy].comm[2]);
	double comm_iy = MIN(thunkd->opt_data[iy].comm[1], thunkd->opt_data[iy].comm[2]);
	for(int s = 1; s <= thunkd->servers; s++) {
		if (comm_iy < thunkd->opt_data[iy].comm[s]) {
			comm_iy2 = comm_iy;
			comm_iy = thunkd->opt_data[iy].comm[s];
		}
	}
	comm_iy = comm_iy2 - comm_iy;

	double cost_x = thunkd->f * thunkd->opt_data[ix].pnts + 
					thunkd->g * comm_ix;
	double cost_y = thunkd->f * thunkd->opt_data[iy].pnts +
					thunkd->g * comm_iy;

	if (cost_x == cost_y) return 0;
	if (cost_x < cost_y) return 1;
	return -1;
}

decl_qsort_p_cmp(sort_remaining_decreasing_low_comm, x, y, thunk) {
	int iy = *(int*)y;
	int ix = *(int*)x;
	thunk_sort_lagrange *thunkd = (thunk_sort_lagrange*)thunk;
	double low_comm_ix = DBL_MAX;
	for(int s = 1; s <= thunkd->servers; s++) {
		if (low_comm_ix > thunkd->opt_data[ix].comm[s])
			low_comm_ix = thunkd->opt_data[ix].comm[s];
	}
	double low_comm_iy = DBL_MAX;
	for(int s = 1; s <= thunkd->servers; s++) {
		if (low_comm_iy > thunkd->opt_data[iy].comm[s])
			low_comm_iy = thunkd->opt_data[iy].comm[s];
	}

	if (low_comm_ix == low_comm_iy) return 0;
	if (low_comm_ix < low_comm_iy) return 1;
	return -1;
}

void remove_double_processed_items(int servers, optimization_data_s *opt_data, int opt_atu, 
	int x[servers][opt_atu], double f, double g) {

	// remove additional runs of multiple processed items
	int iters = 0;
	int dbl_proc_cnt;
	int dbl_processed[opt_atu];
	do {
		dbl_proc_cnt = 0;
		double server_mkspan[servers];
		memset(server_mkspan, 0, sizeof server_mkspan);

		for(int cell = 0; cell < opt_atu; cell++) {
			int integrally = 0;

			for(int s = 0; s < servers; s++) {
				if (x[s][cell] > 0) {
					integrally++;
				}
			}

			if (integrally > 1) {
				dbl_processed[dbl_proc_cnt] = cell;
				dbl_proc_cnt++;

				for(int s = 0; s < servers; s++) {
					if (x[s][cell] > 0)
						server_mkspan[s] += opt_data[cell].pnts;
				}
			}
		}

		if (dbl_proc_cnt > 0) {
			int max_s = 0;
			double max_dbl_processed = server_mkspan[0];
			for(int s = 1; s < servers; s++) {
				if (max_dbl_processed < server_mkspan[s]) {
					max_s = s;
					max_dbl_processed = server_mkspan[s];
				}
			}
			// remove all multiple processed items from max_s
			for(int i = 0; i < dbl_proc_cnt; i++) {
				x[max_s][dbl_processed[i]] = 0;
			}
		}

		iters++;
	} while (dbl_proc_cnt > 0);

	// assert method works
	for(int cell = 0; cell < opt_atu; cell++) {
		int integrally = 0;
		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0)
				integrally++;
		}
		assert(integrally <= 1);
	}
	//printf("Remove double iters %d\n", iters);
}

typedef struct {
	double mkspan;
	int id;
} aux_server_mkspan;

int double_increasing_compare(const void* p1, const void* p2)
{ 
   double i1 = *(double*) p1;
   double i2 = *(double*) p2;
   if (i1 > i2) return 1;
   else if (i1 == i2) return 0;
   else return -1;
}

void print_makespan_gauge(int server, double read, double total) {
	const static char *progress_gauge_equal = "========================================";
	const static char *progress_gauge_empty = "                                        ";
	const static int size = 40;
	const static double resolution = 100.0/40.0;
	int percent = ((read*100) / total);
	int pres = percent / resolution;
	fprintf(stderr, "%2d %15.2f [%.*s%.*s] %3d%%\n", server, read, pres, progress_gauge_equal,
		size - pres, progress_gauge_empty, percent);
}

void schedule_non_assigned_items(int servers, optimization_data_s *opt_data, int pairs, 
	int x[servers][pairs], double f, double g) {

	aux_server_mkspan server_mkspan[servers];
	for(int s=0; s< servers; s++) {
		server_mkspan[s].id = s;
		server_mkspan[s].mkspan = 0;
	}

	// find non assigned items
	int rem_count = 0;
	int remaining[pairs];
	for(int p=0; p < pairs; p++) {
		int count = 0;
		for(int s=0; s < servers; s++) {
			if (x[s][p] == 1) {
				server_mkspan[s].mkspan += opt_data[p].pnts;
				count++;
				break;
			}
		}
		if (count == 0) {
			remaining[rem_count] = p;
			rem_count++;
			//printf("Item %d non scheduled.\n", p);
		}
	}

	// sort by decreasing larger cost
	thunk_sort_lagrange sthunk;
	sthunk.f = f;
	sthunk.g = g;
	sthunk.servers = servers;
	sthunk.opt_data = opt_data;
	qsort_p(remaining, rem_count, sizeof(int), sort_remaining_decreasing_larger_cost, &sthunk);

	// sort by increasing makespan
	qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
	double x0 = server_mkspan[servers-1].mkspan;
	//for(int s = 0; s < servers; s++)
	//	print_makespan_gauge(server_mkspan[s].id, server_mkspan[s].mkspan, x0);	

	// put each non scheduled item in the server where it
	// least increase the cost
	for(int i = 0; i < rem_count; i++) {
		int p = remaining[i];
		int minor = -1;
		int minor_aux = -1;
		double minor_cost = DBL_MAX;
		double smaller_mkspan_inc = DBL_MAX;
		for(int s = 0; s < servers; s++) {
			int orig_s = server_mkspan[s].id;
			double cost_increase = g * opt_data[p].comm[orig_s+1];
			double mkspan_increase = opt_data[p].pnts - (x0 - server_mkspan[s].mkspan);
			if (mkspan_increase > 0)
				cost_increase += f * mkspan_increase;
			if (cost_increase < minor_cost ||
			   (cost_increase == minor_cost && mkspan_increase < smaller_mkspan_inc)) {
				minor_cost = cost_increase;
				minor = orig_s;
				minor_aux = s;
				smaller_mkspan_inc = mkspan_increase;
			}
		}

		assert(minor != -1);

		//printf("Item %d scheduled to %d.\n", p, minor);
		x[minor][p] = 1;
		server_mkspan[minor_aux].mkspan += opt_data[p].pnts;
		if (smaller_mkspan_inc > 0) {
			qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
			x0 = server_mkspan[servers-1].mkspan;
		}
	}
}

void improve_transformed_solution(int servers, optimization_data_s *opt_data, int pairs, 
	int x[servers][pairs], double f, double g) {

	// calc server makespan
	double server_mkspan[servers];
	memset(server_mkspan, 0, sizeof server_mkspan);

	int rem_count = 0;
	int where[pairs];
	for(int p=0; p < pairs; p++) {
		for(int s=0; s < servers; s++) {
			if (x[s][p] == 1) {
				where[p] = s;
				server_mkspan[s] += opt_data[p].pnts;
				break;
			}
		}
	}

	double old_x0;
	double x0, x1;

	// identify the larger server
	int x0_s = server_mkspan[0] > server_mkspan[1] ? 0 : 1;
	x1 = MIN(server_mkspan[0], server_mkspan[1]);
	x0 = MAX(server_mkspan[0], server_mkspan[1]);
	for(int s = 2; s < servers; s++) {
		if (x0 < server_mkspan[s]) {
			x1 = x0;
			x0 = server_mkspan[s];
			x0_s = s;
		}
	}
	old_x0 = x0;
	
	for(int p = 0; p < pairs; p++) {
		int minor = -1;
		double minor_cost = 0;
		for(int s = 0; s < servers; s++) {
			if (s == where[p])
				continue;
			
			double new_x0 = x0;
			if (s == x0_s) // removing p from x0_s will reduce makespan
				new_x0 -= MIN(x0-x1, opt_data[p].pnts);
			new_x0 = MAX(new_x0, server_mkspan[s] + opt_data[p].pnts);

			double cost_increase = g * (opt_data[p].comm[s+1] - opt_data[p].comm[where[p]+1]);
			cost_increase += f * (new_x0 - x0);

			if (cost_increase < minor_cost) {
				minor_cost = cost_increase;
				minor = s;
			}
		}

		if (minor != -1) {
			x[minor][p] = 1;
			x[where[p]][p] = 0;
			server_mkspan[minor] += opt_data[p].pnts;
			server_mkspan[where[p]] -= opt_data[p].pnts;
			where[p] = minor;

			x0_s = server_mkspan[0] > server_mkspan[1] ? 0 : 1;
			x1 = MIN(server_mkspan[0], server_mkspan[1]);
			x0 = MAX(server_mkspan[0], server_mkspan[1]);
			for(int s = 2; s < servers; s++) {
				if (x0 < server_mkspan[s]) {
					x1 = x0;
					x0 = server_mkspan[s];
					x0_s = s;
				}
			}
		}
	}
	//printf("Old x0 %f, new x0 %f\n", old_x0, x0);
}

void improve_transformed_solution_old(int servers, optimization_data_s *opt_data, int pairs, 
	int x[servers][pairs], double f, double g) {

	// calc server makespan
	aux_server_mkspan server_mkspan[servers];
	for(int s=0; s< servers; s++) {
		server_mkspan[s].id = s;
		server_mkspan[s].mkspan = 0;
	}

	int rem_count = 0;
	int remaining[pairs];
	for(int p=0; p < pairs; p++) {
		for(int s=0; s < servers; s++) {
			if (x[s][p] == 1) {
				server_mkspan[s].mkspan += opt_data[p].pnts;
				break;
			}
		}
	}

	double old_x0;
	double x0;

	// identify the large server
	qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
	x0 = server_mkspan[servers-1].mkspan;

	old_x0 = x0;
	double makespan_difference;
	int pair_moved[pairs]; // to prevent infinite loop
	memset(pair_moved, 0, sizeof pair_moved);
	int improved;
	do {
		if (servers < 2)
			break;

		improved = 0;
		int cs = servers-1;
		int cs_orig = server_mkspan[cs].id;
		makespan_difference = server_mkspan[cs].mkspan - server_mkspan[cs-1].mkspan;

		// try to remove some item, without increasing any cost
		for(int p=0; p < pairs; p++) {
			if (x[cs_orig][p] == 0 || pair_moved[p] >= 100) // prevent deadlock
				continue;

			int large_mkspan_excess_server = -1;
			double large_mkspan_excess = 0;
			double comm_atu = opt_data[p].comm[cs_orig+1];
			for(int ds = 0; ds < servers; ds++) {
				if (ds == cs)
					continue;

				int ds_orig = server_mkspan[ds].id;

				double mkspan_excess = x0 - server_mkspan[ds].mkspan;
				if (comm_atu >= opt_data[p].comm[ds_orig+1] && 
				   (mkspan_excess >= opt_data[p].pnts)) {
					if (mkspan_excess > large_mkspan_excess) {
						large_mkspan_excess = mkspan_excess;
						large_mkspan_excess_server = ds;
					}
				}
			}

			if (large_mkspan_excess_server != -1) {
				int ds = large_mkspan_excess_server;
				int ds_orig = server_mkspan[ds].id;
				x[ds_orig][p] = 1;
				x[cs_orig][p] = 0;
				server_mkspan[ds].mkspan += opt_data[p].pnts;
				server_mkspan[cs].mkspan -= opt_data[p].pnts;
				pair_moved[p]++;
	
				improved++;

				// identify the large server again
				qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
				x0 = server_mkspan[servers-1].mkspan;
				cs_orig = server_mkspan[cs].id;
				makespan_difference = server_mkspan[cs].mkspan - server_mkspan[cs-1].mkspan;
				//printf("Pair %d from %d to %d\n", p, cs_orig, ds_orig);
				if (makespan_difference < 0) // server is not the larger anymore
					break;
			}

		}

	} while (improved > 0);
	//printf("Makespan reduced from %f to %f: %f%%, ", old_x0, x0, (old_x0-x0)/old_x0);

	// try to reduce communication cost
	// identify the large server
	double comm_reduced = 0;
	qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
	x0 = server_mkspan[servers-1].mkspan;
	old_x0 = x0;
	for(int cs = servers-1; cs >= 0; cs--) {
		int cs_orig = server_mkspan[cs].id;
		for(int p=0; p < pairs; p++) {
			if (x[cs_orig][p] == 0)
				continue;

			double comm_atu = opt_data[p].comm[cs_orig+1];
			for(int ds = 0; ds < servers; ds++) {
				if (ds == cs)
					continue;

				int ds_orig = server_mkspan[ds].id;

				if (comm_atu > opt_data[p].comm[ds_orig+1] && 
				   (x0 - server_mkspan[ds].mkspan >= opt_data[p].pnts)) {
					x[ds_orig][p] = 1;
					x[cs_orig][p] = 0;
					server_mkspan[ds].mkspan += opt_data[p].pnts;
					server_mkspan[cs].mkspan -= opt_data[p].pnts;
					comm_reduced += (comm_atu - opt_data[p].comm[ds_orig+1]);
					//printf("Pair %d from %d to %d\n", p, cs_orig, ds_orig);
		
					x0 = server_mkspan[0].mkspan;
					for(int s = 1; s < servers; s++) {
						if (x0 < server_mkspan[s].mkspan)
							x0 = server_mkspan[s].mkspan;
					}
					break;
				}
			}
		}
	}
	qsort(server_mkspan, servers, sizeof(aux_server_mkspan), &double_increasing_compare);
	x0 = server_mkspan[servers-1].mkspan;


	//printf(" x0 again %f%c, comm. cost reduced %f\n", x0, old_x0 > x0 ? 'x' : ' ', comm_reduced);

	/*printf("-----\n");
	for(int s = 0; s < servers; s++)
		print_makespan_gauge(server_mkspan[s].id, server_mkspan[s].mkspan, x0);	
	exit(1); */
}


void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, int x[servers][opt_atu],
	multiway_histogram_estimate *agg_server, char remove_double_proc,
	double f, double g) {

	int qtd = 10;
	int qtdatu = 0;
	double remaining_mkspan = 0.0;
	int *remaining = g_new(int, qtd);

	double server_mksp[servers+1];
	server_mksp[0] = 0;
	// start with prior cumulated makespan
	for(int s = 0; s < servers; s++) {
		server_mksp[s+1] = 0; //agg_server[s+1].to_pnts;
	}

	if (remove_double_proc) {
		remove_double_processed_items(servers, opt_data, opt_atu, x, f, g);
		schedule_non_assigned_items(servers, opt_data, opt_atu, x, f, g); 
		improve_transformed_solution(servers, opt_data, opt_atu, x, f, g);
	}

	for(int cell = 0; cell < opt_atu; cell++) {
		int integrally = 0;
		int used_server = 0;

		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0) {
				integrally++;
				used_server = s+1;
			}
		}

		histogram_cell *rcell = hr->get_cell(hr, opt_data[cell].xl, opt_data[cell].yl);
		if (integrally == 1) {
			//printf("map\t%d\t%d\t1\n", cell, used_server);

			rcell->place = used_server;
			rcell->copies = 0;
			SET_IN_PLACE(rcell->copies, used_server);

			server_mksp[used_server] += opt_data[cell].pnts;

			for(int c = 0; c < opt_data[cell].rcells_size; c++) {
				histogram_cell *rc = opt_data[cell].rcells[c].cell;
				if (!IS_IN_PLACE(rc->copies, used_server)) {
					SET_IN_PLACE(rc->copies, used_server);
				}
			}
					
			histogram_cell *lc = opt_data[cell].lcell;
			if (!IS_IN_PLACE(lc->copies, used_server)) {
				SET_IN_PLACE(lc->copies, used_server);
			}

		}
		else { // !integrally
			printf("Pair %d in: ", cell);
			for(int s = 0; s < servers; s++) {
				if (x[s][cell] > 0)
					printf("%d ", s);
			}
			printf("\n");
			
			remaining_mkspan += opt_data[cell].pnts;

			remaining[qtdatu] = cell;
			qtdatu++;
			if (qtdatu == qtd) {
				qtd *= 2;
				remaining = g_renew(int, remaining, qtd);
			}
		}
	}

	assert(qtdatu == 0);

	thunk_sort_lagrange sthunk;
	sthunk.servers = servers;
	sthunk.opt_data = opt_data;
	qsort_p(remaining, qtdatu, sizeof(int), sort_remaining_decreasing_low_comm, &sthunk);

	/*for(int i = 0; i < qtdatu; i++) {
		int cell = remaining[i];
		printf("%d: mkspan: %.0f, comm: ", cell, opt_data[cell].pnts);
		for(int s = 1; s <= servers; s++)
			printf("%.0f, ", opt_data[cell].comm[s]);
		printf("\n");
	}*/

	// GREEDY ROUND HEURISTIC
	for(int i = 0; i < qtdatu; i++) {
		
		// find server with minor makespan 
		int minor = 1;
		int minor_pnts = server_mksp[minor];
		for(int s = 2; s <= servers; s++) {
			if (server_mksp[s] < minor_pnts) {
				minor_pnts = server_mksp[s];
				minor = s;
			}
		}

		int cell = remaining[i];

		/*int minor = 1;
		double minor_pnts = opt_data[cell].comm[minor];
		for(int s = 1; s <= servers; s++) {
			if (minor_pnts > opt_data[cell].comm[s]) {
				minor_pnts = opt_data[cell].comm[s];
				minor = s;
			}
		}*/

		histogram_cell *rcell = hr->get_cell(hr, opt_data[cell].xl, opt_data[cell].yl);
		rcell->place = minor;
		rcell->copies = 0;
		SET_IN_PLACE(rcell->copies, minor);

		server_mksp[minor] += opt_data[cell].pnts;

		//printf("set %d to server: %d\n", cell, minor);

		for(int c = 0; c < opt_data[cell].rcells_size; c++) {
			histogram_cell *rc = opt_data[cell].rcells[c].cell;
			if (!IS_IN_PLACE(rc->copies, minor)) {
				SET_IN_PLACE(rc->copies, minor);
			}
		}

		histogram_cell *lc = opt_data[cell].lcell;
		if (!IS_IN_PLACE(lc->copies, minor)) {
			SET_IN_PLACE(lc->copies, minor);
		}
	}

	free(remaining);
}


void set_cell_place_from_partial_x(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data) {

	for(int cell = 0; cell < pairs; cell++) {
		int integrally = 0;
		int used_server = 0;

		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0) {
				integrally++;
				used_server = s+1;
			}
		}

		histogram_cell *rcell = hr->get_cell(hr, opt_data[cell].xl, opt_data[cell].yl);
		if (integrally == 1) {
			//printf("map\t%d\t%d\t1\n", cell, used_server);
			rcell->place = used_server;
		} else {
			rcell->place = 0;
		}
	}
}


