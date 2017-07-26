
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "deps.h"

typedef struct {
	double f;
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

	double cost_x = thunkd->f * thunkd->opt_data[ix].pnts + regret_ix;
	double cost_y = thunkd->f * thunkd->opt_data[iy].pnts + regret_iy;

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

void remove_double_processed_items2(int servers, optimization_data_s *opt_data, int opt_atu, 
	int x[servers][opt_atu], double f) {

	// remove additional runs of multiple processed items
	double adc_load[servers];
	memset(adc_load, 0, sizeof adc_load);

	for(int cell = 0; cell < opt_atu; cell++) {
		int integrally = 0;

		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0) {
				integrally++;
			}
		}

		if (integrally > 1) {
			int stay_on = -1;
			double min_cost = DBL_MAX;
			for(int s = 0; s < servers; s++) {
				if (x[s][cell] > 0) {
					x[s][cell] = 0;
					double scost = f * opt_data[cell].pnts + opt_data[cell].comm[s+1];
					if (scost < min_cost) {
						stay_on = s;
						min_cost = scost;
					}
				}
			}
			x[stay_on][cell] = 1;
		}
	}
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

void remove_double_processed_items(int servers, optimization_data_s *opt_data, int opt_atu, 
	int x[servers][opt_atu], double f) {

	// remove additional runs of multiple processed items
	int iters = 0;
	int dbl_proc_cnt;
	int dbl_processed[opt_atu];
	do {
		dbl_proc_cnt = 0;
		double adc_load[servers];
		memset(adc_load, 0, sizeof adc_load);

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
						adc_load[s] += opt_data[cell].pnts;
				}
			}
		}

		if (dbl_proc_cnt > 0) {
			int max_s = 0;
			double max_dbl_processed = adc_load[0];
			for(int s = 1; s < servers; s++) {
				if (max_dbl_processed < adc_load[s]) {
					max_s = s;
					max_dbl_processed = adc_load[s];
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

void schedule_non_assigned_items(int unassigned_count, int unassigned[unassigned_count], int servers, 
	optimization_data_s *opt_data, double load[servers], int pairs, int where[pairs], 
	int x[servers][pairs], double f) {

	// sort by decreasing larger cost
	thunk_sort_lagrange sthunk;
	sthunk.f = f;
	sthunk.servers = servers;
	sthunk.opt_data = opt_data;
	qsort_p(unassigned, unassigned_count, sizeof(int), sort_remaining_decreasing_regret_cost, &sthunk);

	double x0 = maxdbl(load, 0, servers);
	//for(int s = 0; s < servers; s++)
	//	print_makespan_gauge(s, load[s], x0);	

	// put each non scheduled item in the server where it
	// least increase the cost
	for(int i = 0; i < unassigned_count; i++) {
		int p = unassigned[i];
		int minor = -1;
		double minor_cost = DBL_MAX;
		double smaller_mkspan_inc = DBL_MAX;
		for(int s = 0; s < servers; s++) {
			double cost_increase = opt_data[p].comm[s+1];
			double mkspan_increase = opt_data[p].pnts - (x0 - load[s]);
			if (mkspan_increase > 0)
				cost_increase += f * mkspan_increase;
			if (cost_increase < minor_cost ||
			   (cost_increase == minor_cost && mkspan_increase < smaller_mkspan_inc)) {
				minor_cost = cost_increase;
				minor = s;
				smaller_mkspan_inc = mkspan_increase;
			}
		}

		assert(minor != -1);

		//printf("Item %d scheduled to %d.\n", p, minor);
		x[minor][p] = 1;
		where[p] = minor;
		load[minor] += opt_data[p].pnts;
		if (load[minor] > x0)
			x0 = load[minor];
	}
	/*for(int p = 0; p < pairs; p++) {
		assert(where[p] != -1);
		assert(x[where[p]][p] == 1);
	}*/
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
	double load[servers], int where[pairs], int x[servers][pairs], double f) {

	double old_x0;
	double mostloaded[3];
	int mostloadedidx[3];
	three_most_loaded_server(servers, load, mostloaded, mostloadedidx);
	old_x0 = mostloaded[0];
	double x0 = mostloaded[0];
	
	int exchanged_pairs_cnt = 0;
	double total_cost_reduction = 0;
	for(int j1 = 0; j1 < pairs; j1++) {
		int best_exchange = -1;
		double best_exchange_cost = 0;
		int pj1 = where[j1];

		// don't need to check j1,j2 and j2,j1
		for(int j2 = j1+1; j2 < pairs; j2++) {
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
			double new_server_inc = load[server_inc] + pdiff;
			if (new_x0 < new_server_inc) new_x0 = new_server_inc;

			// cost increase/decrease
			double cost_inc = - opt_data[j1].comm[pj1+1]
							  - opt_data[j2].comm[pj2+1]
							  + opt_data[j1].comm[pj2+1]
							  + opt_data[j2].comm[pj1+1];

			double worthiness = f * (new_x0 - x0) + cost_inc;

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
			load[pj1] += - opt_data[j1].pnts + opt_data[j2].pnts;
			load[pj2] += - opt_data[j2].pnts + opt_data[j1].pnts;

			three_most_loaded_server(servers, load, mostloaded, mostloadedidx);
			x0 = mostloaded[0];
		}
	}
	//printf("Exchanged %d, z decreased %f, Old x0 %f, new x0 %f\n", exchanged_pairs_cnt, total_cost_reduction, old_x0, x0);
}

void improve_transformed_solution_fm(int servers, optimization_data_s *opt_data, int pairs, 
	double load[servers], int where[pairs], int x[servers][pairs], double f) {

	typedef struct {
		int key;
		short usedc[64]; //TODO: This limit to max 64 servers
		UT_hash_handle hh;
	} used_rcells;
	used_rcells *rcells = NULL;

	// compute the comm cost
	// lcell accounts only in lcell->place by design
	for(int p = 0; p < pairs; p++) {
		for(int c = 0; c < opt_data[p].rcells_size; c++) {
			right_opt_data *rcell = &opt_data[p].rcells[c];
			int rid = rcell->xr + (rcell->yr<<16);
			used_rcells *r;
			HASH_FIND_INT(rcells, &rid, r);
			if (r == NULL) {
				r = g_new0(used_rcells, 1);
				r->key = rid;
				HASH_ADD_INT(rcells, key, r);
			}
			r->usedc[where[p]]++;
		}
	}

	double old_x0;
	double x0, x1;

	// identify the larger server load
	int x0_s = load[0] > load[1] ? 0 : 1;
	x1 = MIN(load[0], load[1]);
	x0 = MAX(load[0], load[1]);
	for(int s = 2; s < servers; s++) {
		if (x0 < load[s]) {
			x1 = x0;
			x0 = load[s];
			x0_s = s;
		} else if (x1 < load[s]) {
			x1 = load[s];
		}
	}
	old_x0 = x0;
	
	for(int p = 0; p < pairs; p++) {
		int minor = -1;
		double minor_cost = 0;
		for(int s = 0; s < servers; s++) {
			if (s == where[p]) // move p to where it alread is
				continue;
			
			double new_x0 = x0;
			if (s == x0_s) // removing p from x0_s will reduce makespan
				new_x0 -= MIN(x0-x1, opt_data[p].pnts);
			new_x0 = MAX(new_x0, load[s] + opt_data[p].pnts);

			double cost_increase = 0;
			for(int c = 0; c < opt_data[p].rcells_size; c++) {
				right_opt_data *rcell = &opt_data[p].rcells[c];
				int rid = rcell->xr + (rcell->yr<<16);
				used_rcells *r;
				HASH_FIND_INT(rcells, &rid, r);
				assert(r != NULL);
				// remove c from where[p] will reduce costs?
				if (r->usedc[where[p]] == 1) // only this pair put r there?
					cost_increase -= ((rcell->cell->place-1) == where[p] ? 0.0 : rcell->cell->points);
				// adding c to s will increase cost?
				if (r->usedc[s] == 0) // increase if r are not in s
					cost_increase += ((rcell->cell->place-1) == s ? 0.0 : rcell->cell->points);
			}

			// lcell
			histogram_cell *lcell = opt_data[p].lcell;
			cost_increase -= (lcell->place-1) == where[p] ? 0.0 : lcell->points;
			cost_increase += (lcell->place-1) == s ? 0.0 : lcell->points;
		
			cost_increase += f * (new_x0 - x0);
			if (cost_increase < minor_cost) {
				minor_cost = cost_increase;
				minor = s;
			}
		}

		if (minor != -1) {
			printf("Moving %d from %d to %d reduce cost: %f\n", p, where[p], minor, minor_cost);
			int s = minor;
			double new_x0 = x0;
			if (s == x0_s) // removing p from x0_s will reduce makespan
				new_x0 -= MIN(x0-x1, opt_data[p].pnts);
			new_x0 = MAX(new_x0, load[s] + opt_data[p].pnts);
			printf("\t makespan from %f to %f\n", x0, new_x0);

			histogram_cell *lcell = opt_data[p].lcell;
			printf("\t l cell actual %f new %f\n", 
				(lcell->place-1) == where[p] ? 0.0 : lcell->points,
				(lcell->place-1) == s ? 0.0 : lcell->points);
	
			for(int c = 0; c < opt_data[p].rcells_size; c++) {
				right_opt_data *rcell = &opt_data[p].rcells[c];
				int rid = rcell->xr + (rcell->yr<<16);
				used_rcells *r;
				HASH_FIND_INT(rcells, &rid, r);
				assert(r != NULL);
				// remove c from where[p] will reduce costs?
				double actual = 0;
				if (r->usedc[where[p]] == 1) // only this pair put r there?
					actual = ((rcell->cell->place-1) == where[p] ? 0.0 : rcell->cell->points);
				// adding c to s will increase cost?
				double next = 0;
				if (r->usedc[s] == 0) // increase if r are not in s
					next = ((rcell->cell->place-1) == s ? 0.0 : rcell->cell->points);
				printf("\t\trid: %d from %f to %f\n", rid, actual, next);
			}
	
			x[minor][p] = 1;
			x[where[p]][p] = 0;
			load[minor] += opt_data[p].pnts;
			load[where[p]] -= opt_data[p].pnts;

			for(int c = 0; c < opt_data[p].rcells_size; c++) {
				right_opt_data *rcell = &opt_data[p].rcells[c];
				int rid = rcell->xr + (rcell->yr<<16);
				used_rcells *r;
				HASH_FIND_INT(rcells, &rid, r);
				assert(r != NULL);
				r->usedc[where[p]]--;
				assert(r->usedc[where[p]] >= 0);
				r->usedc[minor]++;
			}

			where[p] = minor;

			x0_s = load[0] > load[1] ? 0 : 1;
			x1 = MIN(load[0], load[1]);
			x0 = MAX(load[0], load[1]);
			for(int s = 2; s < servers; s++) {
				if (x0 < load[s]) {
					x1 = x0;
					x0 = load[s];
					x0_s = s;
				} else if (x1 < load[s]) {
					x1 = load[s];
				}
			}
		}
	}

	used_rcells *current, *tmp;
	HASH_ITER(hh, rcells, current, tmp) {
	    HASH_DEL(rcells, current);
    	g_free(current);
	}

	//printf("Old x0 %f, new x0 %f\n", old_x0, x0);
	/*for(int p = 0; p < pairs; p++) {
		assert(where[p] != -1);
		assert(x[where[p]][p] == 1);
	}*/
}

void improve_transformed_solution(int servers, optimization_data_s *opt_data, int pairs, 
	double load[servers], int where[pairs], int x[servers][pairs], double f) {

	double old_x0;
	double x0, x1;

	// identify the larger server
	int x0_s = load[0] > load[1] ? 0 : 1;
	x1 = MIN(load[0], load[1]);
	x0 = MAX(load[0], load[1]);
	for(int s = 2; s < servers; s++) {
		if (x0 < load[s]) {
			x1 = x0;
			x0 = load[s];
			x0_s = s;
		} else if (x1 < load[s]) {
			x1 = load[s];
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
			new_x0 = MAX(new_x0, load[s] + opt_data[p].pnts);

			double cost_increase = (opt_data[p].comm[s+1] - opt_data[p].comm[where[p]+1]);
			cost_increase += f * (new_x0 - x0);

			if (cost_increase < minor_cost) {
				minor_cost = cost_increase;
				minor = s;
			}
		}

		if (minor != -1) {
			x[minor][p] = 1;
			x[where[p]][p] = 0;
			load[minor] += opt_data[p].pnts;
			load[where[p]] -= opt_data[p].pnts;
			where[p] = minor;

			x0_s = load[0] > load[1] ? 0 : 1;
			x1 = MIN(load[0], load[1]);
			x0 = MAX(load[0], load[1]);
			for(int s = 2; s < servers; s++) {
				if (x0 < load[s]) {
					x1 = x0;
					x0 = load[s];
					x0_s = s;
				} else if (x1 < load[s]) {
					x1 = load[s];
				}
			}
		}
	}

	//printf("Old x0 %f, new x0 %f\n", old_x0, x0);
	/*for(int p = 0; p < pairs; p++) {
		assert(where[p] != -1);
		assert(x[where[p]][p] == 1);
	}*/
}

void get_where_for_pair(int pairs, int servers, int x[servers][pairs], int where[pairs]) {
	for(int cell = 0; cell < pairs; cell++) {
		int count = 0;
		int used_server = 0;

		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0) {
				count++;
				used_server = s;
			}
		}

		if (count == 1) {
			where[cell] = used_server;
		}
		else {
			printf("where[%d] = %d\n", cell, count);
			where[cell] = -1;
		}
	}
}

void set_cell_place_from_partial_x(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data, int where[pairs]) {

	for(int cell = 0; cell < pairs; cell++) {
		histogram_cell *rcell = hr->get_cell(hr, opt_data[cell].xl, opt_data[cell].yl);
		if (where[cell] > -1) {
			//printf("map\t%d\t%d\t1\n", cell, used_server);
			int used_server = where[cell]+1;
			rcell->place = used_server;
			rcell->copies = 0;
			SET_IN_PLACE(rcell->copies, used_server);

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
		} else {
			printf("where %d set to %d.\n", cell, where[cell]);
			rcell->place = 0;
		}
	}
}

void set_cell_place_from_partial_x_nowh(dataset_histogram *hr, int servers, int pairs, 
	int x[servers][pairs], optimization_data_s *opt_data) {
	int where[pairs];
	get_where_for_pair(pairs, servers, x, where);
	set_cell_place_from_partial_x(hr, servers, pairs, x, opt_data, where);
}

void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int pairs, int x[servers][pairs],
	multiway_histogram_estimate *agg_server, double f, bool ignore_multiple,
	bool improve_exchange) {

	double load[servers];
	memset(load, 0, sizeof(double)*servers);

	int where[pairs];
	memset(where, -1, sizeof where);

	if (ignore_multiple)
		remove_double_processed_items2(servers, opt_data, pairs, x, f);

	// find non assigned items
	int unassigned_count = 0;
	int unassigned[pairs];
	for(int p=0; p < pairs; p++) {
		int s1;
		int count = 0;
		for(int s=0; s < servers; s++) {
			if (x[s][p] == 1) {
				x[s][p] = 0; // disable for the case of multiple processed
				s1 = s;
				count++;
			}
		}

		if (count == 1) {
			load[s1] += opt_data[p].pnts;
			x[s1][p] = 1; // reenable
			where[p] = s1;
		}
		else {
			unassigned[unassigned_count] = p;
			unassigned_count++;
			//printf("Item %d non scheduled.\n", p);
		}
	}

	if (unassigned_count > 0)
		schedule_non_assigned_items(unassigned_count, unassigned, servers, opt_data, load, pairs, where, x, f);
	if (improve_exchange)
		improve_transformed_solution_exchange_pairs(servers, opt_data, pairs, load, where, x, f);
	else
		improve_transformed_solution(servers, opt_data, pairs, load, where, x, f);

	/*for(int p = 0; p < pairs; p++) {
		assert(where[p] != -1);
		assert(x[where[p]][p] == 1);
	}*/

	set_cell_place_from_partial_x(hr, servers, pairs, x, opt_data, where);
}

