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
#include "structs.h"
#include "cpphash.h"

void histogram_alloc(dataset_histogram *dh, int xqtd, int yqtd);
void bs_optimize_hr(dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
	int opt_atu, multiway_histogram_estimate *agg_server);
void lagrange_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int pairs, multiway_histogram_estimate *agg_server,
	double dualvalues[pairs]);
void lpi_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	bool only_root_node);
void totalize_estimate(dataset_histogram *hr, multiway_histogram_estimate *estimate, 
	int servers, optimization_data_s *opt_data, int pairs);
void histogram_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers, 
	int *mkspan, int *totalcomm);
void reset_opt_data_copies(optimization_data_s *opt_data, int opt_atu);

void histogram_set_data_grid(dataset *ds, dataset_histogram_persist *hp);
histogram_cell *histogram_get_cell(dataset_histogram *dh, unsigned short x, unsigned short y);

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		printf("%s input.bin\n", argv[0]);
		exit(1);
	}

	FILE *f = fopen(argv[1], "r");
	if (!f) {
		printf("Error opening %s\n", argv[1]);
		exit(1);
	}

	int opt_atu, servers;
	fread(&opt_atu, sizeof(int), 1, f);
	fread(&servers, sizeof(int), 1, f);

	dataset ds;
	dataset_histogram *hr = &ds.metadata.hist;
	fread(hr, sizeof(dataset_histogram), 1, f);

	int size;
	fread(&size, sizeof(int), 1, f);
	char hdata[size];
	fread(hdata, size, 1, f);

	dataset_histogram_persist *hp = (dataset_histogram_persist*)hdata;
	// code from dataset_set_histogram
	ds.metadata.hist.htype = hp->htype;
	ds.metadata.hist.get_cell = histogram_get_cell;
	histogram_set_data_grid(&ds, hp);

	optimization_data_s opt_data[opt_atu];	
	CppMap *map = map_create();

	for(int i = 0; i < opt_atu; i++) {
		fread(&opt_data[i], sizeof(optimization_data_s), 1, f);
		
		opt_data[i].lcell = (histogram_cell*)malloc(sizeof(histogram_cell));
		fread(opt_data[i].lcell, sizeof(histogram_cell), 1, f);

		opt_data[i].comm = (double*)malloc(sizeof(double)*(servers+1));
		fread(opt_data[i].comm, sizeof(double), servers+1, f);

		opt_data[i].rcells = (right_opt_data*)malloc(sizeof(right_opt_data)*
			opt_data[i].rcells_size);
		for(int j = 0; j < opt_data[i].rcells_size; j++) {
			fread(&opt_data[i].rcells[j], sizeof(right_opt_data), 1, f);

			histogram_cell hc_aux;
			fread(&hc_aux, sizeof(histogram_cell), 1, f);

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
	fread(agg_server, sizeof(multiway_histogram_estimate), servers+1, f);

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

	reset_opt_data_copies(opt_data, opt_atu);


	#ifdef LAGRANGE
	// calculate an upper bound using greedy algorithm
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server);
	// call lagrangian optimization
	reset_opt_data_copies(opt_data, opt_atu);
	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, NULL);
	#endif

	#ifdef MIP
	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, false);
	#endif

	// print the summary 
	//printf("\n\nSummary:\n");
	multiway_histogram_estimate estimate[servers];
	memset(estimate, 0, sizeof estimate);

	totalize_estimate(hr, estimate, servers, opt_data, opt_atu);
	int makespan, totalcomm;
	histogram_print_estimate("server", estimate, servers, &makespan, &totalcomm);

    char *aux = getenv("LP_COMM");
    double g = aux ? atof(aux) : 1.0;
    aux = getenv("LP_MKSP");
    double fm = aux ? atof(aux) : servers;

	printf("\nValues for obj, with costs observing FM rules\n");
	printf("Obj: %.0f, Makespan: %d, Comm: %d\n\n",
		fm * makespan + g * totalcomm, makespan, totalcomm);

	fclose(f);
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


void histogram_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers, 
	int *mkspan, int *totalcomm) {

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

	if (mkspan)
		*mkspan = max_to_pnts;
	if (totalcomm)
		*totalcomm = io_pnts;
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

