
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

decl_qsort_p_cmp(sort_remaining_decreasing_regret_cost, x, y, thunk) {
	int iy = *(int*)y;
	int ix = *(int*)x;
	thunk_sort_lagrange *thunkd = (thunk_sort_lagrange*)thunk;

	double comm_ix_max = thunkd->opt_data[ix].comm[1];
	double comm_ix_min = thunkd->opt_data[ix].comm[1];
	for(int s = 2; s <= thunkd->servers; s++) {
		if (comm_ix_max < thunkd->opt_data[ix].comm[s])
			comm_ix_max = thunkd->opt_data[ix].comm[s];
		if (comm_ix_min > thunkd->opt_data[ix].comm[s])
			comm_ix_min = thunkd->opt_data[ix].comm[s];
		
	}
	double regret_ix = comm_ix_max - comm_ix_min;

	double comm_iy_max = MAX(thunkd->opt_data[iy].comm[1], thunkd->opt_data[iy].comm[2]);
	double comm_iy_min = MIN(thunkd->opt_data[iy].comm[1], thunkd->opt_data[iy].comm[2]);
	for(int s = 1; s <= thunkd->servers; s++) {
		if (comm_iy_max < thunkd->opt_data[iy].comm[s])
			comm_iy_max = thunkd->opt_data[iy].comm[s];
		if (comm_iy_min > thunkd->opt_data[iy].comm[s])
			comm_iy_min = thunkd->opt_data[iy].comm[s];
	}
	double regret_iy = comm_iy_max - comm_iy_min;

	double cost_x = thunkd->f * thunkd->opt_data[ix].pnts + 
					thunkd->g * regret_ix;
	double cost_y = thunkd->f * thunkd->opt_data[iy].pnts +
					thunkd->g * regret_iy;

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
	/*for(int cell = 0; cell < opt_atu; cell++) {
		int integrally = 0;
		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0)
				integrally++;
		}
		assert(integrally <= 1);
	}*/
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
	qsort_p(remaining, rem_count, sizeof(int), sort_remaining_decreasing_regret_cost, &sthunk);

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

void three_most_loaded_server(int servers, double server_load[servers], double values[3], int indxs[3]) {
	// identify the most loaded server
	assert(servers >= 3);
	int large = 0;
	if (server_load[1] > server_load[large]) large = 1;
	if (server_load[2] > server_load[large]) large = 2;
	int smaller = 0;
	if (server_load[1] < server_load[smaller]) smaller = 1;
	if (server_load[2] < server_load[smaller]) smaller = 2;
	int middle = 3 - large - smaller;
	values[0] = server_load[large];
	values[1] = server_load[middle];
	values[2] = server_load[smaller];
	indxs[0] = large;
	indxs[1] = middle;
	indxs[2] = smaller;
	for(int s = 3; s < servers; s++) {
		if (server_load[s] > values[0]) {
			indxs[2] = indxs[1];
			indxs[1] = indxs[0];
			indxs[0] = s;
			values[2] = values[1];
			values[1] = values[0];
			values[0] = server_load[s];
		} else if (server_load[s] > values[1]) {
			indxs[2] = indxs[1];
			indxs[1] = s;
			values[2] = values[1];
			values[1] = server_load[s];
		} else if (server_load[s] > values[2]) {
			indxs[2] = s;
			values[2] = server_load[s];
		}
	}
}

void improve_transformed_solution_exchange_pairs(int servers, optimization_data_s *opt_data, int pairs, 
	int x[servers][pairs], double f, double g) {

	// calc server makespan
	double server_load[servers];
	memset(server_load, 0, sizeof server_load);

	int rem_count = 0;
	int where[pairs];
	for(int p=0; p < pairs; p++) {
		for(int s=0; s < servers; s++) {
			if (x[s][p] == 1) {
				where[p] = s;
				server_load[s] += opt_data[p].pnts;
				break;
			}
		}
	}

	double old_x0;
	double mostloaded[3];
	int mostloadedidx[3];
	three_most_loaded_server(servers, server_load, mostloaded, mostloadedidx);
	old_x0 = mostloaded[0];
	double x0 = mostloaded[0];
	
	int exchanged_pairs_cnt = 0;
	double total_cost_reduction = 0;
	for(int j1 = 0; j1 < pairs; j1++) {
		int best_exchange = -1;
		double best_exchange_cost = 0;
		int pj1 = where[j1];

		for(int j2 = 0; j2 < pairs; j2++) {
			if (j1 == j2) continue;
			if (where[j1] == where[j2]) continue;

			int pj2 = where[j2];

			// increase makespan?
			double aux_mostloaded[3];
			aux_mostloaded[0] = mostloaded[0];
			aux_mostloaded[1] = mostloaded[1];
			aux_mostloaded[2] = mostloaded[2];
			int server_inc = opt_data[j1].pnts > opt_data[j2].pnts ? pj2 : pj1;
			int server_dec = server_inc == pj1 ? pj2 : pj1;
			double pdiff = fabs(opt_data[j1].pnts - opt_data[j2].pnts);
			for(int i = 0; i < 3; i++) {
				if (server_inc == mostloadedidx[i])
					aux_mostloaded[i] += pdiff;
			}
			for(int i = 0; i < 3; i++) {
				if (server_dec == mostloadedidx[i])
					aux_mostloaded[i] -= pdiff;
			}
			double new_x0 = aux_mostloaded[0];
			if (new_x0 < aux_mostloaded[1]) new_x0 = aux_mostloaded[1];
			if (new_x0 < aux_mostloaded[2]) new_x0 = aux_mostloaded[2];

			// the server for which the load increases becomes the new x0?
			double new_server_inc = server_load[server_inc] + pdiff;
			if (new_x0 < new_server_inc) new_x0 = new_server_inc;

			// cost increase/decrease
			double cost_inc = - opt_data[j1].comm[pj1+1]
							  - opt_data[j2].comm[pj2+1]
							  + opt_data[j1].comm[pj2+1]
							  + opt_data[j2].comm[pj1+1];

			double worthiness = f * (new_x0 - x0) + g * cost_inc;

			if (worthiness < best_exchange_cost) {
				best_exchange_cost = worthiness;
				best_exchange = j2;
			}
		}

		if (best_exchange != -1) {
			exchanged_pairs_cnt++;
			total_cost_reduction += best_exchange_cost;

			int j2 = best_exchange;
			int pj2 = where[j2];
			assert(x[pj1][j1]);
			assert(x[pj2][j2]);
			x[pj1][j1] = 0;
			x[pj2][j2] = 0;
			x[pj2][j1] = 1;
			x[pj1][j2] = 1;
			where[j1] = pj2;
			where[j2] = pj1;
			server_load[pj1] += - opt_data[j1].pnts + opt_data[j2].pnts;
			server_load[pj2] += - opt_data[j2].pnts + opt_data[j1].pnts;

			three_most_loaded_server(servers, server_load, mostloaded, mostloadedidx);
			x0 = mostloaded[0];
		}
	}
	//printf("Exchanged %d, z decreased %f, Old x0 %f, new x0 %f\n", exchanged_pairs_cnt, total_cost_reduction, old_x0, x0);
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
		improve_transformed_solution_exchange_pairs(servers, opt_data, opt_atu, x, f, g);
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

	//assert(qtdatu == 0);

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

