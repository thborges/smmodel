
#ifdef __MACH__
//#include <config.h>
#endif

#include <stdlib.h>
#include <glibwrap.h>
//#include <multiway-join.h>
//#include <dataset.h>
//#include <utils.h>
#include "structs.h"
#include <limits.h>
#include <glpk.h>
#include <round.h>
#include "sm.h"

#define ENABLE_CPLEX 1

#if !ENABLE_CPLEX
#include <glpk_cplex_wrap.h>
#else
#include <ilcplex/cplex.h>
#endif

#define LP_WITH_NAMED_VARS


void lpi_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	bool only_root_node) {

	char *aux = getenv("LP_COMM");
	double g = aux ? atof(aux) : 1.0;
	aux = getenv("LP_MKSP");
	double f = aux ? atof(aux) : servers;
	printf("Tradeoff Comm: %f, Mkspan %f\n", g, f);
	printf("\n============ Building CPLEX model =============\n");

	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status = 0;

	env = CPXopenCPLEX(&status);
	lp = CPXcreateprob(env, &status, "dgeo");

	CPXchgobjsen(env, lp, CPX_MIN);

	// create model cols and rows
	int cols = opt_atu * servers + 1;
	double obj[cols];
	double ub[cols];
	char ctype[cols];

	#ifdef LP_WITH_NAMED_VARS
	char *cnames[cols];
	#else
	char **cnames = NULL;
	#endif

	for(int i = 0; i < opt_atu; i++) {
		for(int s = 1; s <= servers; s++) {
			int catu = i*servers + s-1;
			obj[catu] = opt_data[i].comm[s];
			ub[catu] = g;
			ctype[catu] = 'I'; // integer

			#ifdef LP_WITH_NAMED_VARS
			char *name = malloc(sizeof(char)*20);
			sprintf(name, "map(%d,%d)", i, s);
			cnames[catu] = name;
			#endif
		}
	}

	//makespan
	obj[cols-1] = f; 
	ub[cols-1] = CPX_INFBOUND;
	ctype[cols-1] = 'C'; // continuous

	#ifdef LP_WITH_NAMED_VARS
	cnames[cols-1] = "mksp";
	#endif

	status = CPXnewcols(env, lp, cols, obj, NULL, ub, ctype, cnames);
	
	// rows - constraints: only one server
	int oos_rows = opt_atu;
	int oos_rowsnz = opt_atu * servers;
	int oos_matbeg[oos_rows];
	double oos_rhs[oos_rows];
	char oos_sense[oos_rows];
	int oos_matind[oos_rowsnz];
	double oos_matval[oos_rowsnz];

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
	int mks_matbeg[mks_rows];
	double mks_rhs[mks_rows];
	char mks_sense[mks_rows];
	int mks_matind[mks_rowsnz];
	double mks_matval[mks_rowsnz];

	for(int s = 0; s < servers; s++) {
		int begin = s * (opt_atu+1);
		mks_matbeg[s] = begin;
		mks_sense[s] = 'G';
		mks_rhs[s] = 0; //agg_server[s+1].to_pnts;
	
		for(int i = 0; i < opt_atu; i++) {
			mks_matind[begin + i] = i*servers + s;
			mks_matval[begin + i] = - opt_data[i].pnts;
		}

		mks_matind[begin + opt_atu] = cols-1;
		mks_matval[begin + opt_atu] = 1.0;
	}
	status = CPXaddrows(env, lp, 0, mks_rows, mks_rowsnz, mks_rhs, 
		mks_sense, mks_matbeg, mks_matind, mks_matval, NULL, NULL);


	// provide an initial solution
	int mcnt = 1;
	int nzcnt = opt_atu * servers;
	int beg[mcnt]; beg[0] = 0;
	int varindices[nzcnt];
	double values[nzcnt];
	int actual = 0;
	for(int i = 0; i < opt_atu; i++) {
		for(int s = 1; s <= servers; s++) {
			int index = i*servers + s-1;
			varindices[actual] = index;

			histogram_cell *rcell = hr->get_cell(hr, opt_data[i].xl, opt_data[i].yl);
			values[actual] = rcell->place == s ? 1.0 : 0.0;
			actual++;
		}
	}
	status = CPXaddmipstarts(env, lp, mcnt, nzcnt, beg, varindices, values, NULL, NULL);

	// only root node
	if (only_root_node)
		status = CPXsetintparam(env, CPXPARAM_MIP_Limits_Nodes, 0);

	// screen output
	status = CPXsetintparam(env, 1035, CPX_ON);

	// threads
	// status = CPXsetintparam(env, 1067, 1);

	// stop at mip gap
	status = CPXsetdblparam (env, CPXPARAM_MIP_Tolerances_MIPGap, (double)0.005);
	status = CPXsetdblparam (env, CPX_PARAM_WORKMEM, 60*1024.0); // at most 1G RAM
	status = CPXsetintparam (env, CPX_PARAM_NODEFILEIND, 2); // write node files to disk, uncompressed
	status = CPXsetintparam (env, CPX_PARAM_VARSEL, 2); // use strong branching
	status = CPXsetintparam (env, CPX_PARAM_MIPEMPHASIS, 3); // Emphasize best bound

	// cuts
	status = CPXsetintparam (env, CPX_PARAM_CLIQUES, 3);
	status = CPXsetintparam (env, CPX_PARAM_COVERS, 3);
	status = CPXsetintparam (env, CPX_PARAM_DISJCUTS, 3);
	status = CPXsetintparam (env, CPX_PARAM_FLOWCOVERS, 2);
	status = CPXsetintparam (env, CPX_PARAM_FLOWPATHS, 2);
	status = CPXsetintparam (env, CPX_PARAM_FRACCUTS, 2);
	status = CPXsetintparam (env, CPX_PARAM_GUBCOVERS, 2);
	status = CPXsetintparam (env, CPX_PARAM_IMPLBD, 2);
	status = CPXsetintparam (env, CPX_PARAM_MIRCUTS, 2);
	status = CPXsetintparam (env, CPX_PARAM_MCFCUTS, 2);
	status = CPXsetintparam (env, CPX_PARAM_ZEROHALFCUTS, 2);

	// optimize
	status = CPXmipopt(env, lp);

	// get best objective bound
	double best_objval;
	status = CPXgetbestobjval(env, lp, &best_objval);
	printf("\n\nBest bound objective value: %f\n", best_objval);

	// print solutions in pool objective value

	printf("%3s %-10s %15s %15s %10s\n", "Sol", "Name", "Best", "Objective", "Gap%");
	int numsolns = CPXgetsolnpoolnumsolns (env, lp);
	for(int i = 0; i < numsolns; i++) {
		double objval;
		char name[10];
		int surplus;
		CPXgetsolnpoolsolnname(env, lp, name, 10, &surplus, i);
		status = CPXgetsolnpoolobjval (env, lp, i, &objval);
		if (i == numsolns-1)
			sprintf(name, "%s", "m1");
		printf("%3d %-10s %15.2f %15.2f %10.2f\n", i, name, best_objval, objval, fabs(objval-best_objval)/objval*100);
	}
	printf("\n");

	// char filename[100];
	// sprintf(filename, "debug.lp", opt_atu);
	// CPXwriteprob(env, lp, filename, NULL);

	// fill map
	int cur_numcols = CPXgetnumcols(env, lp);
	double x[cur_numcols];
	int solstat;
	double objval;
	status = CPXsolution(env, lp, &solstat, &objval, x, NULL, NULL, NULL);
	printf("Status = %d, Solution status = %d, Objective value = %f\n", status,
		solstat, objval);

	int map[opt_atu][servers+1];
	memset(map, 0, sizeof map);

	for(int cell = 0; cell < opt_atu; cell++) {
		for(int server = 1; server <= servers; server++) {
			int index = cell*servers + server-1;
			map[cell][server] = x[index] > 0 ? 1 : 0;
		}
	}

	for(int cell = 0; cell < opt_atu; cell++) {
		int used_server = 0;
		for(int server = 1; server <= servers; server++) {
			if (map[cell][server] == 1) {
				used_server = server;
				break;
			}
		}
		assert(used_server > 0);

		histogram_cell *rcell = hr->get_cell(hr, opt_data[cell].xl, opt_data[cell].yl);
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

	double final_mkspan, final_comm;
	double Zheur = get_sm_objective(hr, opt_data, opt_atu, f, g, servers, 1, NULL, NULL, &final_mkspan, &final_comm);
	printf("LPI Z\tMkspan\tComm\n%.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

	// free problem
	#ifdef LP_WITH_NAMED_VARS
	for(int catu = 0; catu < cols-1; catu++)
		free(cnames[catu]);
	#endif

	status = CPXfreeprob(env, &lp);
	status = CPXcloseCPLEX(&env);
}

