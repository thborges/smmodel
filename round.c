
#include <float.h>
#include <structs.h>
#include <glibwrap.h>
#include <stdlib.h>

typedef struct {
	int servers;
	optimization_data_s* opt_data;
} thunk_sort_lagrange;

// declare and call qsort_r on OSX and Linux, in a portable way
#ifdef __MACH__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(void *thunk, const void *x, const void *y)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, thunk, compar)
#elif __linux__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(const void *x, const void *y, void *thunk)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, compar, thunk)
#endif

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

void lp_optimize_hr_round_decreasing_low_comm(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, int x[servers][opt_atu]) {

	int qtd = 10;
	int qtdatu = 0;
	int *remaining = g_new(int, qtd);

	for(int cell = 0; cell < opt_atu; cell++) {
		int integrally = 0;
		int used_server = 0;

		for(int s = 0; s < servers; s++) {
			if (x[s][cell] > 0) {
				integrally++;
				used_server = s+1;
			}
		}

		if (integrally == 1) {
			//printf("map\t%d\t%d\t1\n", cell, used_server);

			histogram_cell *rcell = GET_HISTOGRAM_CELL(hr, opt_data[cell].xl, opt_data[cell].yl);
			rcell->place = used_server;
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

		}
		else { // !integrally || integrally > 1
			remaining[qtdatu] = cell;
			qtdatu++;
			if (qtdatu == qtd) {
				qtd *= 2;
				remaining = g_renew(int, remaining, qtd);
			}
		}
	}

	printf("Remaining items: %d\n", qtdatu);
	thunk_sort_lagrange sthunk;
	sthunk.servers = servers;
	sthunk.opt_data = opt_data;
	qsort_p(remaining, qtdatu, sizeof(int), sort_remaining_decreasing_low_comm, &sthunk);

	printf("id\tmkspan\tcomm\t");
	for(int s = 1; s <= servers; s++)
		printf("%d\t", s);
	printf("\n");
	
	for(int i = 0; i < qtdatu; i++) {
		int cell = remaining[i];
		printf("%d\t%.0f\t\t", cell, opt_data[cell].pnts);
		for(int s = 1; s <= servers; s++)
			printf("%.0f\t", opt_data[cell].comm[s]);
		printf("\n");
	}

	printf("Round result:\nid\tserver\n");
	// GREEDY ROUND HEURISTIC
	for(int i = 0; i < qtdatu; i++) {
		
		// find server with minor points
		/*int minor = 1;
		int minor_pnts = estimate[minor].to_pnts;
		for(int s = 2; s <= servers; s++) {
			if (estimate[s].to_pnts < minor_pnts) {
				minor_pnts = estimate[s].to_pnts;
				minor = s;
			}
		}*/

		int cell = remaining[i];

		int minor = 1;
		double minor_pnts = opt_data[cell].comm[minor];
		for(int s = 1; s <= servers; s++) {
			if (minor_pnts > opt_data[cell].comm[s]) {
				minor_pnts = opt_data[cell].comm[s];
				minor = s;
			}
		}

		histogram_cell *rcell = GET_HISTOGRAM_CELL(hr, opt_data[cell].xl, opt_data[cell].yl);
		rcell->place = minor;
		SET_IN_PLACE(rcell->copies, minor);

		printf("%d\t%d\n", cell, minor);

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
