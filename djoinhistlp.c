
#include <stdlib.h>
#include <glibwrap.h>
#include <limits.h>
#include <glpk.h>
#include "deps.h"
#include "round.h"
#include "sm.h"

#if !ENABLE_CPLEX
#include <glpk_cplex_wrap.h>
#else
#include <ilcplex/cplex.h>
#endif

#define LP_WITH_NAMED_VARS

const char *fcellsname = "cells.txt";
const char *fpointsname = "points.txt";
const char *fnetworname = "networ.txt";

void lp_generate_data(optimization_data_s *opt_data, int opt_atu, int servers, multiway_histogram_estimate *agg_server) {
	
	FILE *fcells = fopen(fcellsname, "w");
	FILE *fpoints = fopen(fpointsname, "w");
	FILE *fnetwor = fopen(fnetworname, "w");

	fprintf(fcells, "data;\n\nset SERVERS := ");
	fprintf(fpoints, "data;\nparam cells_points :=\n");
	fprintf(fnetwor, "data; param cells_networ: ");
	for(int s = 1; s <= servers; s++) {
		fprintf(fnetwor, " %d", s);
		fprintf(fcells, " %d", s);
	}

	fprintf(fcells, ";\n\nparam priour_mkspn :=\n");
	for(int s = 1; s <= servers; s++) {
		fprintf(fcells, "%d\t%f\n", s, agg_server[s].to_pnts);
	}

	fprintf(fnetwor, " := \n");
	fprintf(fcells, ";\n\nset CELLS := \n");
	for(int c = 0; c < opt_atu; c++) {
		fprintf(fpoints, "%d\t%d\n", c, (int)opt_data[c].pnts);
		fprintf(fcells, "%d\n", c);
		fprintf(fnetwor, "%d", c);
		for(int s = 1; s <= servers; s++) {
			fprintf(fnetwor, "\t%d", (int)opt_data[c].comm[s]);
		}
		fprintf(fnetwor, "\n");
	}

	fprintf(fcells, ";\n\nend;\n");
	fprintf(fpoints, ";\n\nend;\n");
	fprintf(fnetwor, ";\n\nend;\n");

	fclose(fpoints);
	fclose(fnetwor);
	fclose(fcells);
}

void lp_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	double f, double dualvalues[opt_atu]) {

	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status = 0;

	env = CPXopenCPLEX(&status);
	lp = CPXcreateprob(env, &status, "dgeo");

	CPXchgobjsen(env, lp, CPX_MIN);

	// create model cols and rows
	int cols = opt_atu * servers + 1;
	double *obj = g_new(double, cols);
	double *ub = g_new(double, cols);

	#ifdef LP_WITH_NAMED_VARS
	char **cnames = g_new(char*, cols);
	#else
	char **cnames = NULL;
	#endif

	for(int i = 0; i < opt_atu; i++) {
		for(int s = 1; s <= servers; s++) {
			int catu = i*servers + s-1;
			obj[catu] = opt_data[i].comm[s];
			ub[catu] = 1;

			#ifdef LP_WITH_NAMED_VARS
			char *name = malloc(sizeof(char)*20);
			sprintf(name, "c%d.%d", i, s);
			cnames[catu] = name;
			#endif
		}
	}

	//makespan
	obj[cols-1] = f; 
	ub[cols-1] = CPX_INFBOUND;

	#ifdef LP_WITH_NAMED_VARS
	cnames[cols-1] = "mksp";
	#endif

	status = CPXnewcols(env, lp, cols, obj, NULL, ub, NULL, cnames);
	
	// rows - constraints: only one server
	int oos_rows = opt_atu;
	int oos_rowsnz = opt_atu * servers;
	int *oos_matbeg = g_new(int, oos_rows);
	double *oos_rhs = g_new(double, oos_rows);
	char *oos_sense = g_new(char, oos_rows);
	int *oos_matind = g_new(int, oos_rowsnz);
	assert(oos_matind != NULL);
	double *oos_matval = g_new(double, oos_rowsnz);
	assert(oos_matval != NULL);

	for(int i = 0; i < opt_atu; i++) {
		int begin = i*servers;
		oos_matbeg[i] = begin;
		oos_sense[i] = 'E';
		oos_rhs[i] = 1.0;
		for(int s = 1; s <= servers; s++) {
			oos_matind[begin + s-1] = begin + s-1;
			oos_matval[begin + s-1] = 1.0;
		}
	}
	status = CPXaddrows(env, lp, 0, oos_rows, oos_rowsnz, oos_rhs, 
		oos_sense, oos_matbeg, oos_matind, oos_matval, NULL, NULL);

	// rows - constraints: makespan
	int mks_rows = servers;
	int mks_rowsnz = (opt_atu+1) * servers;

	int *mks_matbeg = g_new(int, mks_rows);
	double *mks_rhs = g_new(double, mks_rows);
	char *mks_sense = g_new(char, mks_rows);
	int *mks_matind = g_new(int, mks_rowsnz);
	assert(mks_matind != NULL);
	double *mks_matval = g_new(double, mks_rowsnz);
	assert(mks_matval != NULL);

	for(int s = 0; s < servers; s++) {
		int begin = s * (opt_atu+1);
		mks_matbeg[s] = begin;
		mks_sense[s] = 'G';
		//enable to input makespan of previous step
		//mks_rhs[s] = agg_server[s+1].to_pnts;
		mks_rhs[s] = 0;
	
		for(int i = 0; i < opt_atu; i++) {
			mks_matind[begin + i] = i*servers + s;
			mks_matval[begin + i] = - opt_data[i].pnts;
		}

		mks_matind[begin + opt_atu] = cols-1;
		mks_matval[begin + opt_atu] = 1.0;
	}
	status = CPXaddrows(env, lp, 0, mks_rows, mks_rowsnz, mks_rhs, 
		mks_sense, mks_matbeg, mks_matind, mks_matval, NULL, NULL);

	// screen output
	// status = CPXsetintparam(env, 1035, CPX_ON);

	// threads
	status = CPXsetintparam(env, 1067, 1);

	// optimize
	status = CPXlpopt(env, lp);
	
	// get opt data
	//int cur_numrows = CPXgetnumrows(env, lp);
	int cur_numcols = CPXgetnumcols(env, lp);
	double *x = g_new(double, cur_numcols);
	int solstat;
	double objval;

	//char filename[100];
	//sprintf(filename, "debug.lp", opt_atu);
	//CPXwriteprob(env, lp, filename, NULL);

	status = CPXsolution(env, lp, &solstat, &objval, x, NULL, NULL, NULL);

	int xaux[servers][opt_atu];
	for(int cell = 0; cell < opt_atu; cell++) {
		for(int s = 0; s < servers; s++) {
			int catu = cell*servers + s;
			xaux[s][cell] = x[catu] > 0.0 ? 1 : 0;
		}
	}

	double Zheur, final_mkspan, final_comm;
	set_cell_place_from_partial_x_nowh(hr, servers, opt_atu, xaux, opt_data);
    Zheur = get_sm_objective(hr, opt_data, opt_atu, f, servers, 1, NULL, NULL, &final_mkspan, &final_comm);
	printf("Before LP rounding\nZ\tMkspan\tComm\n%.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

	lp_optimize_hr_round_decreasing_low_comm(hr, servers, opt_data, opt_atu, xaux, agg_server, f, true, false);

	Zheur = get_sm_objective(hr, opt_data, opt_atu, f, servers, 1, NULL, NULL, &final_mkspan, &final_comm);
	printf("After LP rounding\nZ\tMkspan\tComm\nSM_LP %.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

	// free problem
	#ifdef LP_WITH_NAMED_VARS
	for(int catu = 0; catu < cols-1; catu++)
		free(cnames[catu]);
	#endif

	// provide dual values for constraints that forces only one server to process an item "= 1"
	if (dualvalues)
		status = CPXgetpi(env, lp, dualvalues, 0, opt_atu);

	status = CPXfreeprob(env, &lp);
	status = CPXcloseCPLEX(&env);
}

