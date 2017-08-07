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
#include "uthash.h"
#include "utils.h"

bool verbose;

void lp_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f, double dualvalues[opt_atu]);

void bs_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f);

void bs_optimize_hr_old(dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
	int opt_atu, multiway_histogram_estimate *agg_server, double trade_off);

void bs_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f);

void lagrange_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_sever,
	double f, double dualvalues[opt_atu], double *m1, double *c1, double *m2, double *c2);

void lpi_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	double f, bool only_root_node, double *mkspan, double *totalcomm);

void lpi_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	double f, char *lname, char *rname, bool only_root_node);

void totalize_estimate(dataset_histogram *hr, multiway_histogram_estimate *estimate, 
	int servers, optimization_data_s *opt_data, int pairs);
void multiway_totalize_estimate(multiway_histogram_estimate *estimate, int servers,
	double *totalpnts, double *totalcomm, 
	double *mkspan, double *max_comm,
	double *stdev_mkspan, double *stdev_comm,
	double *mkspan_gap);

void histogram_alloc(dataset_histogram *dh, int xqtd, int yqtd);
void histogram_print_estimate(char *name, multiway_histogram_estimate *estimate, int servers);
void reset_opt_data_copies(optimization_data_s *opt_data, int opt_atu);

void histogram_set_data_grid(dataset *ds, dataset_histogram_persist *hp);
histogram_cell *histogram_get_cell(dataset_histogram *dh, unsigned short x, unsigned short y);

double runtime_diff_ms(struct timespec *start, struct timespec *end) {
	return ( end->tv_sec - start->tv_sec ) * 1000.0 + (double)( end->tv_nsec - start->tv_nsec ) / 1E6;
}

//const double scale = 100000;
const double scale = 1.0;

void exactpa(double f1, double yf1, double x01, double c1,
		double f2, double yf2, double x02, double c2,
		int level,
		dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
		int opt_atu, multiway_histogram_estimate *agg_server) {
	const char *empty = "                    ";

	double f3 = (c2-c1)/(x01-x02);
	if (f3 <= 0) {
		printf("%.*sPA: cannot continue due f3=%f\n", level, empty, f3);
		return;
	}

	double mkspan, totalcomm;
	reset_opt_data_copies(opt_data, opt_atu);
	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f3, false, &mkspan, &totalcomm);
	double x03 = mkspan/scale;
	double c3 = totalcomm/scale;
	printf("%.*sPA: f3=%f f1: %fx + %f\n", level, empty, f3, x03, c3);

	//find zf for p4
	double zp4 = f3*x01+c1;

	//find zf for p3
	double zp3 = f3*x03+c3;

	// y=ax+b
	double a = (yf1-yf2)/(f1-f2);
	double b = -a * f1 + yf1;
	double zp5 = f3*a + b;

	printf("%.*sPA: zp4: %f, zp3: %f, zp5: %f\n", level, empty, zp4, zp3, zp5);
	if (round(zp4*10) == round(zp3*10)) {
		printf("%.*sPA: case b: found breakpoint f=%f\n", level, empty, f3);
		printf("%.*sGEOG: Segment[(%f,%f),(%f,%f)]\n", level, empty, f1, yf1, f3, zp3);
		printf("%.*sGEOG: Segment[(%f,%f),(%f,%f)]\n", level, empty, f3, zp3, f2, yf2);
	}
	else if (round(zp3*10) == round(zp5*10)) {
		printf("%.*sPA: case c: no breakpoint\n", level, empty);
		printf("%.*sGEOG: Segment[(%f,%f),(%f,%f)]\n", level, empty, f1, yf1, f2, yf2);
	}
	else {
		printf("%.*sPA: case a: recursive call to refine uncertainty area\n", level, empty);
		exactpa(f1, yf1, x01, c1, f3, zp3, x03, c3, level+1, hr, servers, opt_data, opt_atu, agg_server);
		exactpa(f3, zp3, x03, c3, f2, yf2, x02, c2, level+1, hr, servers, opt_data, opt_atu, agg_server);
	}
}

void pa(double f1, double yf1u, double x01u, double c1u, double yf1l,
		double f2, double yf2u, double x02u, double c2u, double yf2l,
		int level,
		dataset_histogram *hr, int servers,	optimization_data_s *opt_data, 
		int opt_atu, multiway_histogram_estimate *agg_server) {
	const char *empty = "                    ";

	if (level > 4) {
		return;
	}

	double f3 = (c2u-c1u)/(x01u-x02u);
	if (f3 <= 0 || isinf(f3) || isnan(f3)) {
		printf("%.*sPA: cannot continue due f3=%f\n", level, empty, f3);
		printf("%.*sPA: Segment[(%f,%f),(%f,%f)]\n", level, empty, f1, yf1u, f2, yf2u);
		return;
	}
	if (!(f3 > f1 && f3 < f2) || !(f3 < f1 && f3 > f2)) {
		printf("%.*sPA: f3 %f not between f1 %f and f2 %f. Will continue with f3=%f\n", 
			level, empty, f3, f1, f2, (f1+f2)/2.0);
		f3 = (f1+f2)/2.0;
	}

	printf("%.*sPA: f3=%f\n", level, empty, f3);

	double lw_mkspan, lw_totalcomm, up_mkspan, up_totalcomm;
	reset_opt_data_copies(opt_data, opt_atu);
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f3);
	reset_opt_data_copies(opt_data, opt_atu);
	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f3, NULL, 
		&lw_mkspan, &lw_totalcomm, &up_mkspan, &up_totalcomm);
	double x03l = lw_mkspan/scale;
	double c3l = lw_totalcomm/scale;
	double x03u = up_mkspan/scale;
	double c3u = up_totalcomm/scale;
	printf("%.*sPA: f3=%f f3: %fx + %f\n", level, empty, f3, x03u, c3u);

	//lower bound
	printf("%.*sPA: Segment[(%f,%f),(%f,%f)]\n", level, empty, f1, yf1l, f3, f3*x03l+c3l);
	printf("%.*sPA: Segment[(%f,%f),(%f,%f)]\n", level, empty, f2, yf2l, f3, f3*x03l+c3l);
	//upper bound
	printf("%.*sPA: Segment[(%f,%f),(%f,%f)]\n", level, empty, f1, f1*x01u+c1u, f3, f3*x03u+c3u);
	printf("%.*sPA: Segment[(%f,%f),(%f,%f)]\n", level, empty, f2, f2*x02u+c2u, f3, f3*x03u+c3u);

	pa(f1, yf1u, x01u, c1u, yf1l, f3, f3*x03u+c3u, x03u, c3u, f3*x03l+c3l, level+1, 
		hr, servers, opt_data, opt_atu, agg_server);
	pa(f3, f3*x03u+c3u, x03u, c3u, f3*x03l+c3l, f2, yf2u, x02u, c2u, yf2l, level+1, 
		hr, servers, opt_data, opt_atu, agg_server);

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

	char *qname = getenv("MW_QNAME");
	if (!qname)
		qname = "none";

	// load instance. This is the M_pa or Q2_3_16M
	const int m = 16;
	const int n = 69;
	int servers = m;
	int opt_atu = n;
	double w[n] = {4,7,8,6,4,6,7,9,8,2,3,1,9,8,4,7,1,7,1,1,2,8,9,8,3,1,3,2,8,6,8,3,7,9,1,2,4,7,2,9,8,9,2,8,5,8,8,5,6,1,4,7,2,8,7,4,3,1,7,9,6,9,3,1,7,6,7,2,4};
	double cost[n][m+1] = {
{0,8,8,8,8,8,6,2,8,8,8,8,8,8,8,8,8},
{0,102,102,102,102,102,102,102,102,102,102,102,102,102,58,44,102},
{0,152,152,152,152,152,152,152,152,152,152,152,152,126,82,94,152},
{0,14,14,14,14,14,14,8,8,14,14,14,14,14,14,14,12},
{0,24,24,24,24,24,22,16,18,24,24,24,24,24,24,16,24},
{0,42,42,42,42,24,20,40,42,42,42,42,42,42,42,42,42},
{0,154,188,188,188,188,188,188,188,188,188,188,188,188,188,130,94},
{0,118,118,118,118,118,118,118,118,118,118,118,118,118,118,60,58},
{0,182,182,182,168,182,182,182,182,182,182,182,182,182,156,98,122},
{0,144,144,136,144,144,144,144,144,144,144,144,144,104,74,118,144},
{0,186,180,186,186,186,186,186,186,186,186,186,142,96,140,186,186},
{0,94,94,94,94,94,94,94,94,94,94,94,50,44,94,94,94},
{0,14,14,14,14,14,14,14,0,14,14,14,14,14,14,14,14},
{0,18,18,18,18,18,18,14,4,18,18,18,18,18,18,18,18},
{0,52,52,52,52,52,6,48,52,52,52,52,52,52,52,52,52},
{0,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,0},
{0,306,306,306,306,306,306,306,306,306,292,306,306,306,306,194,126},
{0,288,288,288,288,288,288,288,288,268,288,288,288,288,134,174,288},
{0,372,372,372,372,372,372,372,328,372,372,372,372,196,218,372,372},
{0,214,214,214,214,214,214,174,214,214,214,214,214,40,214,214,214},
{0,292,292,292,292,292,288,292,292,292,292,292,178,118,292,292,292},
{0,16,16,16,16,16,16,10,6,16,16,16,16,16,16,16,16},
{0,48,48,48,48,48,10,42,48,48,48,48,46,48,48,48,48},
{0,352,348,352,352,352,352,352,352,352,352,352,352,352,352,200,154},
{0,334,334,334,334,334,334,334,334,334,334,334,334,334,152,184,334},
{0,428,428,428,428,428,428,428,428,428,428,428,428,216,244,428,394},
{0,226,226,226,226,226,226,226,226,226,226,226,226,14,226,210,226},
{0,410,410,410,410,410,410,410,410,410,410,410,218,200,404,410,410},
{0,294,294,294,294,294,294,294,294,294,294,230,102,256,294,294,294},
{0,32,32,32,32,24,14,28,32,32,32,32,32,32,32,32,32},
{0,52,52,52,36,28,40,52,52,52,52,52,52,52,52,52,52},
{0,160,160,160,160,160,160,158,160,160,160,160,160,160,160,78,82},
{0,352,352,352,352,352,288,352,352,352,352,352,352,352,286,208,274},
{0,346,346,346,346,276,346,346,346,346,346,346,346,272,206,282,346},
{0,348,348,348,314,348,348,348,348,348,348,348,268,188,274,348,348},
{0,196,196,166,196,196,196,196,196,196,196,196,114,110,196,196,196},
{0,458,458,358,458,458,458,458,458,458,458,364,278,370,458,458,458},
{0,244,286,286,286,286,286,286,286,286,264,162,188,286,286,286,286},
{0,98,98,98,98,98,52,98,98,98,76,68,98,98,98,98,98},
{0,78,78,78,78,78,78,78,78,64,44,48,78,78,78,78,78},
{0,174,174,174,174,174,174,174,174,80,174,174,174,174,174,94,174},
{0,304,304,304,304,304,304,304,240,210,304,304,304,304,236,222,304},
{0,278,278,278,278,278,278,202,214,278,278,278,278,202,212,278,278},
{0,308,308,308,308,308,230,232,308,308,308,308,226,234,308,308,308},
{0,162,162,162,162,162,84,162,162,162,162,158,80,162,162,162,162},
{0,346,346,346,346,256,268,346,346,346,340,252,266,346,346,346,346},
{0,246,246,246,226,156,246,246,224,246,222,152,246,246,246,246,246},
{0,58,58,58,38,58,58,58,58,58,36,42,58,58,58,58,58},
{0,68,68,56,50,68,68,68,68,56,46,68,68,68,68,68,68},
{0,94,168,168,168,168,168,168,168,76,168,168,168,168,168,168,168},
{0,220,296,296,296,296,296,296,232,202,296,296,296,296,296,296,232},
{0,276,276,276,276,276,276,202,214,276,276,276,276,276,276,202,214},
{0,166,166,166,166,166,88,166,166,166,166,166,166,166,78,166,166},
{0,342,342,342,342,254,266,342,342,342,342,342,342,254,256,342,342},
{0,220,220,220,202,132,220,220,220,220,220,220,196,132,220,220,220},
{0,60,60,48,40,60,60,60,60,60,60,54,36,58,60,60,60},
{0,162,162,162,88,162,162,162,162,162,162,162,162,162,74,160,162},
{0,320,320,250,246,320,320,320,320,320,320,320,320,232,232,320,320},
{0,200,182,130,200,200,200,200,200,200,200,200,178,112,200,200,200},
{0,144,144,144,72,144,144,144,74,144,144,144,144,144,144,144,144},
{0,286,286,216,212,286,286,214,216,286,286,286,286,286,286,286,286},
{0,170,152,100,170,170,160,100,170,170,170,170,170,170,170,170,170},
{0,256,256,256,256,256,256,256,186,194,256,256,256,194,194,256,256},
{0,142,142,142,142,142,142,142,72,142,142,142,132,82,142,142,142},
{0,268,268,268,268,268,268,196,198,268,268,266,206,206,268,268,268},
{0,146,146,146,146,146,136,74,146,146,146,144,84,146,146,146,146},
{0,66,66,66,66,66,66,66,66,64,66,66,62,4,66,66,66},
{0,130,130,130,130,130,130,130,130,130,130,128,66,70,130,130,130},
{0,68,68,68,68,68,68,68,68,68,68,64,6,68,68,68,68},
	};

	dataset_histogram hra;
	histogram_alloc(&hra, n, 1);
	//histogram_set_functions_grid(&hra);
	hra.htype = HIST_GRID;
	hra.get_cell = histogram_get_cell;

	histogram_cell lcell;
	histogram_cell rcell;
	right_opt_data ropt;
	ropt.xr = 0;
	ropt.yr = 0;
	ropt.cell = &rcell;

	optimization_data_s opt_data[n];


	for(int a = 0; a < n; a++) {
		opt_data[a].xl = a;
		opt_data[a].yl = 0;
		opt_data[a].lcell = &lcell;
		opt_data[a].pnts = w[a];
		opt_data[a].comm = cost[a];
		opt_data[a].rcells_size = 1;
		opt_data[a].rcells = &ropt;
	}
	dataset_histogram *hr = &hra;

	// open instance
/*	FILE *finst = fopen(argv[1], "r");
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

	map_destroy(map); */

	multiway_histogram_estimate agg_server[servers+1];
	//fread(agg_server, sizeof(multiway_histogram_estimate), servers+1, finst);

	//TODO: disabled
	memset(agg_server, 0, sizeof agg_server);

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

	#ifdef EXACTPA
	{
		// find bounds of f
		double f1 = 0;
		double mkspan, totalcomm;
		reset_opt_data_copies(opt_data, opt_atu);
		lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f1, false, &mkspan, &totalcomm);
		double x01 = mkspan/scale;
		double c1 = totalcomm/scale;
		printf("PA: f1=%f f1: %fx + %f\n", f1, x01, c1);

		double f2=100;
		reset_opt_data_copies(opt_data, opt_atu);
		lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f2, false, &mkspan, &totalcomm);
		double x02 = mkspan/scale;
		double c2 = totalcomm/scale;
		printf("PA: f2=%f f2: %fx + %f\n", f2, x02, c2);

		exactpa(f1, f1*x01+c1, x01, c1, f2, f2*x02+c2, x02, c2,
			1, hr, servers, opt_data, opt_atu, agg_server);
	}
	#endif

	#ifdef PA
	{
		// find bounds of f
		double f1 = 0;
		double mkspan, totalcomm;
		reset_opt_data_copies(opt_data, opt_atu);
		lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f1, false, &mkspan, &totalcomm);
		double x01 = mkspan/scale;
		double c1 = totalcomm/scale;
		printf("PA: f1=%f f1: %fx + %f\n", f1, x01, c1);

		double f2=100;
		reset_opt_data_copies(opt_data, opt_atu);
		bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f2);
		double lw_mkspan, lw_totalcomm;
		double up_mkspan, up_totalcomm;
		reset_opt_data_copies(opt_data, opt_atu);
		lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f2, NULL, 
			&lw_mkspan, &lw_totalcomm, &up_mkspan, &up_totalcomm);
		printf("2: %f %f %f %f\n", lw_mkspan, lw_totalcomm, up_mkspan, up_totalcomm);
		double x02l = lw_mkspan/scale;
		double c2l = lw_totalcomm/scale;
		double x02u = up_mkspan/scale;
		double c2u = up_totalcomm/scale;
		printf("PA: f2=%f f2: %fx + %f\n", f2, x02u, c2u);

		double yf1 = f1*x01+c1;
		double yf2 = f2*x02u+c2u;
		pa(f1, yf1, x01, c1, yf1, f2, yf2, x02u, c2u, f2*x02l+c2l,
			1, hr, servers, opt_data, opt_atu, agg_server);
 

/*		//find zf for p4
		double zp4 = f3*x01+c1;

		//find zf for p3
		double zp3 = f3*x03+c3;

		//equation from (f1,yf1)
		double yf1 = f1*x01+c1;
		double yf2 = f2*x02+c2;
		// y=ax+b
		double a = (yf1-yf2)/(f1-f2);
		double b = -a * f1 + yf1;
		double zp5 = f3*a + b;

		printf("PA: zp4: %f, zp3: %f, zp5: %f\n", zp4, zp3, zp5);
		if (round(zp4*10) == round(zp3*10))
			printf("PA: case b: found breakpoint f=%f\n", f3);
		else if (round(zp3*10) == round(zp5*10))
			printf("PA: case c: no breakpoint\n");
		else {
			printf("PA: case a: recursive call to refine uncertainty area\n");
			pa(f1, yf1, x01, c1, f3, zp3, x03, c3, 1, hr, servers, opt_data, opt_atu, agg_server);
			pa(f3, zp3, x03, c3, f2, yf2, x02, c2, 1, hr, servers, opt_data, opt_atu, agg_server);
		}*/
		

	}
	#endif

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
		lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, current_f, dualvalues, &aux_makespan, &aux_comm);

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
	double totalcomm, mksp, stdevmksp;
	multiway_histogram_estimate estimate[servers+1];

	#ifdef LAGRANGE
	// transform the instance to PA
/*	double wjs[] = {4, 7, 8, 6, 4, 6, 7, 9, 8, 2, 3, 1, 9, 8, 4, 7, 1, 7, 1, 1, 2, 8, 9, 8, 3, 1, 3, 2, 8, 6, 8, 3, 7, 9, 1, 2, 4, 7, 2, 9, 8, 9, 2, 8, 5, 8, 8, 5, 6, 1, 4, 7, 2, 8, 7, 4, 3, 1, 7, 9, 6, 9, 3, 1, 7, 6, 7, 2, 4};
//	double wjs[] = {4,2,3,1,4,1,2,3,5,2,3,3,1,5,4,2,1,2,1,2,2,4,3,5,3,1,3,2,3,1,5,3,3,4,5,3,4,2,2,3,5,4,2,5,5,3,4,5,1,1,4,2,2,1,2,4,3,1,2,2,1,1,5,3,2,1,2,5,4};

	double scale = 1000.0;
	for(int j=0; j < opt_atu; j++) {
		opt_data[j].pnts = round(opt_data[j].pnts / scale);
		//opt_data[j].pnts = wjs[j];
		//opt_data[j].pnts = random()%5+1;
		//printf("%.0f,", opt_data[j].pnts);
		for(int i=0; i < servers; i++) {
			opt_data[j].comm[i+1] = round(opt_data[j].comm[i+1] / scale);
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
	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, dualvalues, &aux_makespan, &aux_comm);

//	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, true, NULL, NULL);
	#endif

	#ifdef MIP
	lpi_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, false, NULL, NULL);
	#endif

	#ifdef MIPFM
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f);
//	reset_opt_data_copies(opt_data, opt_atu);
//	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, NULL, &aux_makespan, &aux_comm);

/*	totalize_estimate(hr, estimate, servers, opt_data, opt_atu);
	multiway_totalize_estimate(estimate, servers, NULL, &totalcomm, 
		&mksp, NULL, &stdevmksp, NULL, NULL);
	printf("----------===================*********************\n");
	printf("FM  pairs servers Z mksp comm stdevmk\n");
	printf("FM %s\t%d\t%d\t%.2f\t%12.0f\t%12.0f\t%12.0f\t%10.2f\n", "", opt_atu, servers, 
		mksp*f+totalcomm, mksp, totalcomm, stdevmksp, 0.0);

	lpi_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, qname, "", true);
*/
	#endif

	#ifdef GR
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f);
	//bs_optimize_hr_old(hr, servers+1, opt_data, opt_atu, agg_server, 0.001);
	#endif

	#ifdef LP
	double dualvalues[opt_atu];
	lp_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, dualvalues);
	#endif

	#ifdef LPI
	bs_optimize_hr(hr, servers+1, opt_data, opt_atu, agg_server, f);
	lagrange_sm_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, NULL, NULL, NULL);
	lpi_optimize_hr(hr, servers, opt_data, opt_atu, agg_server, f, qname, "", false);
	#endif

	struct timespec cf;
	clock_gettime(CLOCK_REALTIME, &cf);
	double runtime = runtime_diff_ms(&cs, &cf);
	printf("runtime: %.2f\n", runtime);

	// print the summary 
	totalize_estimate(hr, estimate, servers, opt_data, opt_atu);
	if (verbose)
		histogram_print_estimate("server", estimate, servers);

	printf("\nValues for obj, with costs observing FM rules\n");
	multiway_totalize_estimate(estimate, servers, NULL, &totalcomm, 
		&mksp, NULL, &stdevmksp, NULL, NULL);
	printf("FM  pairs servers Z mksp comm stdevmk runtime\n");
	printf("FM %s\t%d\t%d\t%.2f\t%12.2f\t%12.2f\t%12.2f\t%10.2f\n", "", opt_atu, servers, 
		mksp*f+totalcomm,mksp, totalcomm, stdevmksp, runtime);

	//fclose(finst);
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

	typedef struct {
		int key;
		uint64_t place;
		UT_hash_handle hh;
	} used_rcells;
	used_rcells *rcells = NULL;

	for(int i = 0; i < pairs; i++) {
		// to_pnts
		histogram_cell *resultcell = hr->get_cell(hr, opt_data[i].xl, opt_data[i].yl);
		estimate[resultcell->place].to_pnts += opt_data[i].pnts;

		// lcell io_pnts
		if (opt_data[i].lcell->place != resultcell->place)
			estimate[resultcell->place].io_pnts += opt_data[i].lcell->points;

		// rcell's
		for(int c = 0; c < opt_data[i].rcells_size; c++) {
			if (opt_data[i].rcells[c].cell->place != resultcell->place) {
				right_opt_data *rcell = &opt_data[i].rcells[c];

				int rid = rcell->xr + (rcell->yr<<16);
				used_rcells *r;
				HASH_FIND_INT(rcells, &rid, r);
				if (r == NULL) {
					r = g_new(used_rcells, 1);
					r->key = rid;
					r->place = 0;
					SET_IN_PLACE64(r->place, resultcell->place);
					HASH_ADD_INT(rcells, key, r);
					estimate[resultcell->place].io_pnts += rcell->cell->points;
				} else if (!IS_IN_PLACE(r->place, resultcell->place)) {
					SET_IN_PLACE64(r->place, resultcell->place);
					estimate[resultcell->place].io_pnts += rcell->cell->points;
				}
			}
		}
	}

	used_rcells *current, *tmp;
	HASH_ITER(hh, rcells, current, tmp) {
	    HASH_DEL(rcells, current);
    	g_free(current);
	}
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

