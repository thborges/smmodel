/*
 * main.c
 *
 *  Created on: 10/23/2016
 *      Author: Thiago Borges de Oliveira 
 */

#include <stdio.h>
#include <stdlib.h>
#include <glibwrap.h>
#include <assert.h>
#include <float.h>
#include "structs.h"
#include "cpphash.h"
#include "utils.h"

bool verbose;

void lp_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f, double dualvalues[opt_atu]);

void bs_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f);

void lagrange_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f, double dualvalues[opt_atu]);

void lpi_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	double f, bool only_root_node);

void totalize_estimate(dataset_histogram *hr, multiway_histogram_estimate *estimate, 
	int servers, optimization_data_s *opt_data, int pairs);
void multiway_totalize_estimate(multiway_histogram_estimate *estimate, int servers,
	double *totalpnts, double *totalcomm, 
	double *mkspan, double *max_comm,
	double *stdev_mkspan, double *stdev_comm,
	double *mkspan_gap);
void histogram_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers);
void reset_opt_data_copies(optimization_data_s *opt_data, int opt_atu);

void histogram_set_data_grid(dataset *ds, dataset_histogram_persist *hp);
histogram_cell *histogram_get_cell(dataset_histogram *dh, unsigned short x, unsigned short y);

double runtime_diff_ms(struct timespec *start, struct timespec *end) {
	return ( end->tv_sec - start->tv_sec ) * 1000.0 + (double)( end->tv_nsec - start->tv_nsec ) / 1E6;
}

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		printf("%s input.bin\n", argv[0]);
		exit(1);
	}

	// check if MW_VERBOSE is set
	char *auxf = getenv("MW_VERBOSE");
	if (auxf)
		verbose = 1;

	// check if MW_F or MW_F0... is set
	auxf = getenv("MW_F");
	double f = auxf ? atof(auxf) : 1;
	printf("Tradeoff f: %f\n",  f);

	// open instance
	FILE *finst = fopen(argv[1], "r");
	if (!finst) {
		printf("Error opening %s\n", argv[1]);
		exit(1);
	}

	int opt_atu, servers;
	fread(&opt_atu, sizeof(int), 1, finst);
	fread(&servers, sizeof(int), 1, finst);

	dataset ds;
	dataset_histogram *hr = &ds.metadata.hist;
	fread(hr, sizeof(dataset_histogram), 1, finst);

	int size;
	fread(&size, sizeof(int), 1, finst);
	char hdata[size];
	fread(hdata, size, 1, finst);

	dataset_histogram_persist *hp = (dataset_histogram_persist*)hdata;
	// code from dataset_set_histogram
	ds.metadata.hist.htype = hp->htype;
	ds.metadata.hist.get_cell = histogram_get_cell;
	histogram_set_data_grid(&ds, hp);

	optimization_data_s opt_data[opt_atu];	
	CppMap *map = map_create();

	for(int i = 0; i < opt_atu; i++) {
		fread(&opt_data[i], sizeof(optimization_data_s), 1, finst);
		
		opt_data[i].lcell = (histogram_cell*)malloc(sizeof(histogram_cell));
		fread(opt_data[i].lcell, sizeof(histogram_cell), 1, finst);

		opt_data[i].comm = (double*)malloc(sizeof(double)*(servers+1));
		fread(opt_data[i].comm, sizeof(double), servers+1, finst);

		opt_data[i].rcells = (right_opt_data*)malloc(sizeof(right_opt_data)*
			opt_data[i].rcells_size);
		for(int j = 0; j < opt_data[i].rcells_size; j++) {
			fread(&opt_data[i].rcells[j], sizeof(right_opt_data), 1, finst);

			histogram_cell hc_aux;
			fread(&hc_aux, sizeof(histogram_cell), 1, finst);

			int id = GET_ROD_ID(&opt_data[i].rcells[j]);
			void *hc = map_get(map, id);
			if (!hc) {
				hc = g_memdup(&hc_aux, sizeof(histogram_cell));
				map_put(map, id, hc);
			}
			opt_data[i].rcells[j].cell = hc;
		}
	}	

	map_destroy(map);

	multiway_histogram_estimate agg_server[servers+1];
	fread(agg_server, sizeof(multiway_histogram_estimate), servers+1, finst);

	// print optimization data
	/*printf("pair\tpoints\tlcomm\trcells");
	for(int s = 1; s <= servers; s++)
		printf("\t%d", s);
	printf("\n");

	for(int i = 0; i < opt_atu; i++) {
		printf("%d\t%f\t%f", i, opt_data[i].pnts, opt_data[i].lcell->points);
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\trc.%d.%d", opt_data[i].rcells[c].xr, opt_data[i].rcells[c].yr);
		}
		printf("\n\t\tplace");
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\t%d", opt_data[i].rcells[c].cell->place);
		}
		printf("\n\t\tpoints");
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\t%f", opt_data[i].rcells[c].cell->points);
		}
		printf("\n");
	}*/
	
	print_instance_and_solution_fo_file(opt_data, opt_atu, servers, NULL);

	reset_opt_data_copies(opt_data, opt_atu);


	#ifdef LAGRANGE_FINDF
	// call lagrangian optimization
	// code to find f
	double smaller_diff = DBL_MAX;
	double most_similar_f;
	double current_f = f;

	for(int i = 0; i < 10; i++) {
		double aux_makespan, aux_comm;
		printf("\e[1;34mFor f=%.10f\e[0m\n", current_f);

		// calculate an upper bound using greedy algorithm
		double dualvalues[opt_atu];
		lp_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, dualvalues);
		reset_opt_data_copies(opt_data, opt_atu);
		lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, current_f, dualvalues); //, &aux_makespan, &aux_comm);

		double aux = fabs(current_f * aux_makespan - aux_comm);
		if (aux < smaller_diff) {
			smaller_diff = aux;
			most_similar_f = current_f;
		}

		double new_f = round(aux_comm/aux_makespan*100000000.0)/100000000.0;
		if (fabs(new_f-current_f) < 0.0001)
			break;
		else
			current_f = new_f;
	}
	printf("f should be: %.10f\n", most_similar_f);
	#endif

	struct timespec cs;
	clock_gettime(CLOCK_REALTIME, &cs);

	#ifdef LAGRANGE
	// transform the instance to PA
/*	double wjs[] = {4, 7, 8, 6, 4, 6, 7, 9, 8, 2, 3, 1, 9, 8, 4, 7, 1, 7, 1, 1, 2, 8, 9, 8, 3, 1, 3, 2, 8, 6, 8, 3, 7, 9, 1, 2, 4, 7, 2, 9, 8, 9, 2, 8, 5, 8, 8, 5, 6, 1, 4, 7, 2, 8, 7, 4, 3, 1, 7, 9, 6, 9, 3, 1, 7, 6, 7, 2, 4};
//	double wjs[] = {4,2,3,1,4,1,2,3,5,2,3,3,1,5,4,2,1,2,1,2,2,4,3,5,3,1,3,2,3,1,5,3,3,4,5,3,4,2,2,3,5,4,2,5,5,3,4,5,1,1,4,2,2,1,2,4,3,1,2,2,1,1,5,3,2,1,2,5,4};

	double multiplier = 1000.0;
	for(int j=0; j < opt_atu; j++) {
		opt_data[j].pnts = round(opt_data[j].pnts / multiplier);
		//opt_data[j].pnts = wjs[j];
		//opt_data[j].pnts = random()%5+1;
		//printf("%.0f,", opt_data[j].pnts);
		for(int i=0; i < servers; i++) {
			opt_data[j].comm[i+1] = round(opt_data[j].comm[i+1] / multiplier);
			if (opt_data[j].comm[i+1] < 10000)
				opt_data[j].comm[i+1] *= 2;
		}
	}*/

	// calculate an upper bound using greedy algorithm
	double aux_makespan, aux_comm;

	double *dualvalues = NULL;
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f);
//	double dualvalues[opt_atu];
//	lp_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, dualvalues);

	reset_opt_data_copies(opt_data, opt_atu);
	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, dualvalues);//, &aux_makespan, &aux_comm);

//	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, true);
	#endif

	#ifdef MIP
	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, false);
	#endif

	struct timespec cf;
	clock_gettime(CLOCK_REALTIME, &cf);
	double runtime = runtime_diff_ms(&cs, &cf);
	printf("runtime: %.2f\n", runtime);

	// print the summary 
	multiway_histogram_estimate estimate[servers];
	memset(estimate, 0, sizeof estimate);

	totalize_estimate(hr, estimate, servers, opt_data, opt_atu);
	if (verbose)
		histogram_print_estimate("server", estimate, servers);

	printf("\nValues for obj, with costs observing FM rules\n");
	double totalcomm, mksp, stdevmksp;
	multiway_totalize_estimate(estimate, servers-1, NULL, &totalcomm, 
		&mksp, NULL, &stdevmksp, NULL, NULL);
	printf("FM %s\t%d\t%d\t%12.0f\t%12.0f\t%12.0f\t%10.2f\n", "", opt_atu, servers, 
		mksp, totalcomm, stdevmksp, runtime);

	fclose(finst);
}

int compare_right_opt_data(const void* x, const void*y) {
	right_opt_data *xc = *(right_opt_data**)x;
	right_opt_data *yc = *(right_opt_data**)y;
	if (xc->xr > yc->xr) return 1;
	if (xc->xr < yc->xr) return -1;
	if (xc->xr == yc->xr) {
		if (xc->yr < yc->yr) return 1;
		else if (xc->yr > yc->yr) return -1;
	}
	return 0; 	
}

void totalize_estimate(dataset_histogram *hr, multiway_histogram_estimate *estimate, 
	int servers, optimization_data_s *opt_data, int pairs) {

	// clear the estimate structure to refill with the optimized schedule
	memset(estimate, 0, sizeof(multiway_histogram_estimate)*(servers+1));

	// print optimization data
	/*printf("pair\tpoints\tlcomm\trcells");
	for(int s = 1; s <= servers; s++)
		printf("\t%d", s);
	printf("\n");

	for(int i = 0; i < pairs; i++) {
		printf("%d\t%f\t%f", i, opt_data[i].pnts, opt_data[i].lcell->points);
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\trc.%d.%d", opt_data[i].rcells[c].xr, opt_data[i].rcells[c].yr);
		}
		printf("\n\t\tplace");
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\t%d", opt_data[i].rcells[c].cell->place);
		}
		printf("\n\t\tpoints");
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			printf("\t%f", opt_data[i].rcells[c].cell->points);
		}
		printf("\n");
	}*/

	int rcells_atu = 0;
	int rcells_count = pairs; // at least!
	right_opt_data **rcells = g_new(right_opt_data*, rcells_count);

	for(int i = 0; i < pairs; i++) {
		// to_pnts
		histogram_cell *resultcell = hr->get_cell(hr, opt_data[i].xl, opt_data[i].yl);
		estimate[resultcell->place].to_pnts += opt_data[i].pnts;

		// lcell io_pnts
		if (opt_data[i].lcell->place != resultcell->place)
			estimate[resultcell->place].io_pnts += opt_data[i].lcell->points;

		// rcell's
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			rcells[rcells_atu] = &opt_data[i].rcells[c];

			rcells_atu++;
			if (rcells_atu >= rcells_count) {
				rcells_count *= 2;
				rcells = g_renew(right_opt_data*, rcells, rcells_count);
			}
		}
	}

	qsort(rcells, rcells_atu, sizeof(right_opt_data*), compare_right_opt_data);
	
	histogram_cell *last_cell = NULL;
	for(int i=0; i < rcells_atu; i++) {
		if (last_cell != rcells[i]->cell) {
			for(int s = 1; s <= servers; s++) {
				if (rcells[i]->cell->place != s && IS_IN_PLACE(rcells[i]->cell->copies, s))
					estimate[s].io_pnts += rcells[i]->cell->points;
			}
			last_cell = rcells[i]->cell;
		}
	}

	g_free(rcells);
};


void histogram_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers) {

	int to_pnts = 0, io_pnts = 0;
	//printf("%6s   Makespan       Comm\n", name);
	//printf("------ ---------- ----------\n");
    int max_to_pnts = (int)estimate[1].to_pnts;
    int max_io_pnts = (int)estimate[1].io_pnts;
    int min_to_pnts = (int)estimate[1].to_pnts;
    int min_io_pnts = (int)estimate[1].io_pnts;
	for(int s = 1; s <= servers; s++) {
		printf("%6d %10d %10d\n", s, (int)estimate[s].to_pnts, (int)estimate[s].io_pnts);		
		to_pnts += (int)estimate[s].to_pnts;
		io_pnts += (int)estimate[s].io_pnts;
        if (max_to_pnts < (int)estimate[s].to_pnts)
            max_to_pnts = (int)estimate[s].to_pnts;
        if (max_io_pnts < (int)estimate[s].io_pnts)
            max_io_pnts = (int)estimate[s].io_pnts;
        if (min_to_pnts > (int)estimate[s].to_pnts)
            min_to_pnts = (int)estimate[s].to_pnts;
        if (min_io_pnts > (int)estimate[s].io_pnts)
            min_io_pnts = (int)estimate[s].io_pnts;
	}
	//printf("------ ---------- ----------\n");
	//printf("total  %10d %10d\n", to_pnts, io_pnts);
	//printf("max    %10d %10d\n", max_to_pnts, max_io_pnts);
	//printf("min    %10d %10d\n", min_to_pnts, min_io_pnts);
}

void reset_opt_data_copies(optimization_data_s *opt_data, int opt_atu) {
	// reset previous copies
	for(int cell = 0; cell < opt_atu; cell++) {
		histogram_cell *lc = opt_data[cell].lcell;
		lc->copies = 0;
		SET_IN_PLACE(lc->copies, lc->place);

		for(int c = 0; c < opt_data[cell].rcells_size; c++) {
			histogram_cell *rc = opt_data[cell].rcells[c].cell;
			rc->copies = 0;
			SET_IN_PLACE(rc->copies, rc->place);
		}
	}
}

void histogram_alloc(dataset_histogram *dh, int xqtd, int yqtd) {
	assert(xqtd > 0 && yqtd > 0 && "X and Y must be greater than zero.");
	grid_histogram_data *ghd = g_new0(grid_histogram_data, 1);
	dh->extra_data = ghd;
	ghd->xqtd = xqtd;
	ghd->yqtd = yqtd;
	ghd->xtics = g_new(double, xqtd+1);
	ghd->ytics = g_new(double, yqtd+1);
	ghd->hcells = g_new0(histogram_cell, xqtd*yqtd);

	for(int x = 0; x < xqtd; x++) {
		for(int y = 0; y < yqtd; y++) {
			histogram_cell *cell = GET_GRID_HISTOGRAM_CELL(ghd, x, y);
			cell->x = x;
			cell->y = y;
		}
	}
}

void histogram_set_data_grid(dataset *ds, dataset_histogram_persist *hp) {

	ds->metadata.mbr = hp->mbr;

	grid_histogram_data *ghd_orig = (grid_histogram_data*)&hp->extra_data[0];

	histogram_alloc(&ds->metadata.hist, ghd_orig->xqtd, ghd_orig->yqtd);

	grid_histogram_data *ghd = (grid_histogram_data*)ds->metadata.hist.extra_data;
	ghd->xsize = ghd_orig->xsize;
	ghd->ysize = ghd_orig->ysize;

	int pos = sizeof(grid_histogram_data);

	double *xtics = (double*)&hp->extra_data[pos];
	memcpy(ghd->xtics, xtics, sizeof(double)*(ghd_orig->xqtd+1));
	pos += sizeof(double) * (ghd_orig->xqtd+1);

	double *ytics = (double*)&hp->extra_data[pos];
	memcpy(ghd->ytics, ytics, sizeof(double)*(ghd_orig->yqtd+1));
	pos += sizeof(double) * (ghd_orig->yqtd+1);

	histogram_cell *hcells = (histogram_cell*)&hp->extra_data[pos];
	memcpy(ghd->hcells, hcells, sizeof(histogram_cell) * ghd_orig->yqtd * ghd_orig->xqtd);
}

histogram_cell *histogram_get_cell(dataset_histogram *dh, unsigned short x, unsigned short y) {
	grid_histogram_data *ghd = (grid_histogram_data*)dh->extra_data;
	return GET_GRID_HISTOGRAM_CELL(ghd, x, y);
}

double get_multiway_estimate_to_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].to_pnts;
}

double get_multiway_estimate_io_pnts(const void *data, const int n) {
	return ((multiway_histogram_estimate*)data)[n].io_pnts;
}

void multiway_totalize_estimate(multiway_histogram_estimate *estimate, int servers,
	double *totalpnts, double *totalcomm, 
	double *mkspan, double *max_comm,
	double *stdev_mkspan, double *stdev_comm,
	double *mkspan_gap) {

	double to_pnts = 0, io_pnts = 0;
    double max_to_pnts = estimate[1].to_pnts;
    double max_io_pnts = estimate[1].io_pnts;
	for(int s = 1; s <= servers; s++) {
		to_pnts += estimate[s].to_pnts;
		io_pnts += estimate[s].io_pnts;
        if (max_to_pnts < estimate[s].to_pnts)
            max_to_pnts = estimate[s].to_pnts;
        if (max_io_pnts < estimate[s].io_pnts)
            max_io_pnts = estimate[s].io_pnts;
	}

	if (totalpnts)
		*totalpnts = to_pnts;
	if (totalcomm)
		*totalcomm = io_pnts;
	if (mkspan)
		*mkspan = max_to_pnts;
	if (max_comm)
		*max_comm = max_io_pnts;
	if (stdev_mkspan)
		*stdev_mkspan = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_to_pnts);
	if (stdev_comm)
		*stdev_comm = stdevd_ex(estimate, 1, servers+1, get_multiway_estimate_io_pnts);
	if (mkspan_gap) {
		double min_mksp = (to_pnts / servers);
		*mkspan_gap = (max_to_pnts - min_mksp) / min_mksp * 100;
	}
}

