
#include <stdlib.h>
#include <glibwrap.h>
#include <limits.h>
#include <glpk.h>
#include "deps.h"

#define ENABLE_CPLEX 1

#if !ENABLE_CPLEX
#include <glpk_cplex_wrap.h>
#else
#include <ilcplex/cplex.h>
#endif

void write_solution_found(int servers, int opt_atu, 
	int map[opt_atu][servers+1], char *lname, char *rname);

typedef struct {
	short xr;
	short yr;
	histogram_cell *cell;
} rcell_aux;

int compare_rcell_id_int(const void* x, const void*y) {
	rcell_aux *xc = (rcell_aux*)x;
	rcell_aux *yc = (rcell_aux*)y;
	if (xc->xr > yc->xr) return 1;
	if (xc->xr < yc->xr) return -1;
	if (xc->xr == yc->xr) {
		if (xc->yr < yc->yr) return 1;
		else if (xc->yr > yc->yr) return -1;
	}
	return 0; 	
}

int filenameidx = 0;
char *fpattern = "loadbalancef_test_%d.dat";

void lp_generate_data_int(optimization_data_s *opt_data, int opt_atu, int servers, multiway_histogram_estimate *agg_server, double f) {
	
	char filename[100];
	sprintf(filename, fpattern, filenameidx++);
	FILE *fcells = fopen(filename, "w");

	fprintf(fcells, "data;\n\nset SERVERS := ");
	for(int s = 1; s <= servers; s++) {
		fprintf(fcells, " %d", s);
	}

	fprintf(fcells, ";\n\nparam priour_mkspn :=\n");
	for(int s = 1; s <= servers; s++) {
		fprintf(fcells, "%d\t%.10f\n", s, agg_server[s].to_pnts);
	}

	int rcells_atu = 0;
	int rcells_count = 10;
	rcell_aux *rcells = g_new(rcell_aux, rcells_count);
	fprintf(fcells, ";\n\nset LCELLS :=\n");
	for(int p = 0; p < opt_atu; p++) {
		fprintf(fcells, " %d", opt_data[p].xl + (opt_data[p].yl<<16));
		for(int r = 0; r < opt_data[p].rcells_size; r++) {
			rcells[rcells_atu].xr = opt_data[p].rcells[r].xr;
			rcells[rcells_atu].yr = opt_data[p].rcells[r].yr;
			rcells[rcells_atu].cell = opt_data[p].rcells[r].cell;

			rcells_atu++;
			if (rcells_atu >= rcells_count) {
				rcells_count *= 2;
				rcells = g_renew(rcell_aux, rcells, rcells_count);
			}
		}
	}
	fprintf(fcells, ";\n\nset RCELLS :=\n");
	qsort(rcells, rcells_atu, sizeof(rcell_aux), compare_rcell_id_int);
	rcell_aux catu = rcells[0];
	for(int r = 1; r < rcells_atu; r++) {
		if (catu.xr != rcells[r].xr || catu.yr != rcells[r].yr) {
			fprintf(fcells, " %d", catu.xr + (catu.yr<<16));
			catu = rcells[r];
		}
	}
	fprintf(fcells, " %d;\n\n", catu.xr + (catu.yr<<16));

	fprintf(fcells, "param : PAIRS : pair_left pair_weight :=\n");
	for(int p = 0; p < opt_atu; p++) {
		fprintf(fcells, "%d\t%d\t%.10f\n", p, opt_data[p].xl + (opt_data[p].yl<<16), opt_data[p].pnts);
	}
	fprintf(fcells, ";\n\n");

	fprintf(fcells, "param pair_right default 0 (tr)\n");
	for(int p = 0; p < opt_atu; p++) {
		fprintf(fcells, ":\t%d :=\n", p);
		for(int r = 0; r < opt_data[p].rcells_size; r++) {
			fprintf(fcells, "%d\t1\n", opt_data[p].rcells[r].xr + (opt_data[p].rcells[r].yr<<16));
		}
		fprintf(fcells, "\n");
	}
	fprintf(fcells, ";\n\n");


	fprintf(fcells, "param lcell_network(tr) :\n");
	for(int s = 1; s <= servers; s++) {
		fprintf(fcells, "\t%d", s);
	}
	fprintf(fcells, " :=\n");

	for(int p = 0; p < opt_atu; p++) {
		fprintf(fcells, "%d", opt_data[p].xl + (opt_data[p].yl << 16));
		for(int s = 1; s <= servers; s++) {
			histogram_cell *lc = opt_data[p].lcell;
			fprintf(fcells, "\t%.10f", (lc->place != s) ? opt_data[p].lcell->points : 0.0); 
		}
		fprintf(fcells, "\n");
	}

	fprintf(fcells, ";\n\nparam rcell_network(tr) :\n");
	for(int s = 1; s <= servers; s++) {
		fprintf(fcells, "\t%d", s);
	}
	fprintf(fcells, " :=\n");

	catu = rcells[0];
	for(int r = 1; r < rcells_atu; r++) {
		if (catu.xr != rcells[r].xr || catu.yr != rcells[r].yr) {
			fprintf(fcells, "%d", catu.xr + (catu.yr << 16));
			for(int s = 1; s <= servers; s++) {
				fprintf(fcells, "\t%.10f", (catu.cell->place != s) ? catu.cell->points : 0.0);
			}
			fprintf(fcells, "\n");
			catu = rcells[r];
		}
	}
	fprintf(fcells, "%d", catu.xr + (catu.yr << 16));
	for(int s = 1; s <= servers; s++)
		fprintf(fcells, "\t%.10f", (catu.cell->place != s) ? catu.cell->points : 0.0);
	fprintf(fcells, ";\n");

	fprintf(fcells, "\n param f := %.10f;\n", f);

	fprintf(fcells, "end;\n");

	g_free(rcells);
	fclose(fcells);
}

void lpi_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int opt_atu, multiway_histogram_estimate *agg_server,
	double f, char *lname, char *rname, bool only_root_node) {
   
	int map[opt_atu][servers+1];
	memset(map, 0, sizeof map);

	glp_tran *tran = glp_mpl_alloc_wksp();
    glp_prob *prob = glp_create_prob();
    glp_iocp iocp;
    glp_bfcp bfcp;
    const char *out_dpy = NULL;
    const char *in_file = "loadbalancef.mod";

    glp_get_bfcp(prob, &bfcp);
    glp_init_iocp(&iocp);
	iocp.presolve = GLP_ON;
	iocp.gmi_cuts = GLP_ON;
	iocp.mir_cuts = GLP_ON;
	iocp.cov_cuts = GLP_ON;
	iocp.clq_cuts = GLP_ON;
	iocp.fp_heur = GLP_ON;
	iocp.bt_tech = GLP_BT_BPH;
	iocp.br_tech = GLP_BR_PCH;

	glp_scale_prob(prob, GLP_SF_AUTO);

	lp_generate_data_int(opt_data, opt_atu, servers, agg_server, f);

    glp_mpl_read_model(tran, in_file, 1);
	char filename[100];
	sprintf(filename, fpattern, filenameidx-1);
    glp_mpl_read_data(tran, filename);
    glp_mpl_generate(tran, out_dpy);
    glp_mpl_build_prob(tran, prob);

	#if ENABLE_CPLEX
	// write to cplex format
	glp_write_lp(prob, NULL, "cplex.lp");
	
	// call cplex solver
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status = 0;

	env = CPXopenCPLEX(&status);
	lp = CPXcreateprob(env, &status, "dgeo");

	#define CPXPARAM_ScreenOutput 1035
	status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
	status = CPXsetintparam (env, CPX_PARAM_MIPDISPLAY, 3);

	//status = CPXsetintparam (env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);
	status = CPXsetdblparam (env, CPXPARAM_MIP_Tolerances_MIPGap, (double)1e-4);

	status = CPXsetdblparam (env, CPX_PARAM_WORKMEM, 1024.0); // at most 1G RAM
	printf("Set workmem status %d\n", status);
	status = CPXsetintparam (env, CPX_PARAM_NODEFILEIND, 2); // write node files to disk, uncompressed
	status = CPXsetintparam (env, CPX_PARAM_VARSEL, 3); // use strong branching

	printf("Loading cplex.lp into CPLEX solver.\n");
	status = CPXreadcopyprob (env, lp, "cplex.lp", NULL);


	// provide an initial solution
	printf("Configuring initial solution\n");
	int mcnt = 1;
	int nzcnt = opt_atu * servers;
	int beg[mcnt]; beg[0] = 0;
	int varindices[nzcnt];
	double values[nzcnt];
	int actual = 0;

	int totalcols = CPXgetnumcols(env, lp);
	for(int c = 0; c < totalcols; c++) {
		char *colnames[1];
		char name[20];
		int aux;
		status = CPXgetcolname(env, lp, colnames, name, 20, &aux, c, c);
		const char *search = "map";
		if (status == 0 && 
			strncmp(name, search, strlen(search)) == 0) {
			int i, s;
			sscanf(&name[4], "%d,%d", &i, &s);
			//printf("Found %s i=%d and s=%d\n", name, i, s);
			varindices[actual] = c;
			histogram_cell *rcell = hr->get_cell(hr, opt_data[i].xl, opt_data[i].yl);
			values[actual] = rcell->place == s ? 1.0 : 0.0;
			actual++;
		}
	}
	printf("Setting initial solution m1\n");
	status = CPXaddmipstarts(env, lp, mcnt, nzcnt, beg, varindices, values, NULL, NULL);

	// read initial solution from file if it exists
	char solfilename[200];
	sprintf(solfilename, "not_sol/%s_%s_%d_%d.txt", lname, rname, servers, opt_atu);
	FILE *fsol = fopen(solfilename, "r");
	if (fsol) {
		
		actual = 0;
		for(int i = 0; i < opt_atu; i++) {
			for(int s = 1; s <= servers; s++) {
				int ii, ss, v;
				fscanf(fsol, "%d:%d=%d", &ii, &ss, &v);
				char name[100];
				sprintf(name, "map(%d,%d)", ii, ss);
				int index;
				status = CPXgetcolindex(env, lp, name, &index);
				assert(status == 0);
				varindices[actual] = index;

				values[actual] = v;
				actual++;
			}
		}
		status = CPXaddmipstarts(env, lp, mcnt, nzcnt, beg, varindices, values, NULL, NULL);
		fclose(fsol);
	}

	// only root node
	if (only_root_node)
		status = CPXsetintparam(env, CPXPARAM_MIP_Limits_Nodes, 0);

	// stop at mip gap
	status = CPXsetdblparam (env, CPXPARAM_MIP_Tolerances_MIPGap, (double)0.0005);
	status = CPXsetdblparam (env, CPX_PARAM_WORKMEM, 60*1024.0); // at most 1G RAM
	status = CPXsetintparam (env, CPX_PARAM_NODEFILEIND, 2); // write node files to disk, uncompressed
	status = CPXsetintparam (env, CPX_PARAM_VARSEL, 2); // use strong branching
	status = CPXsetintparam (env, CPX_PARAM_MIPEMPHASIS, 3); // Emphasize best bound
	//status = CPXsetintparam (env, CPX_PARAM_MIPEMPHASIS, 1); // Emphasize feasibility
	//status = CPXsetintparam (env, CPX_PARAM_MIPEMPHASIS, 4); // Emphasize hidden feasibility

	status = CPXsetdblparam (env, CPX_PARAM_CUTSFACTOR, (double)1.0); // disable all cuts
	printf("Param CutsFactor set status %d\n", status);

	if (only_root_node) {
		// disable probing on variables: time consuming at start
		status = CPXsetintparam (env, CPX_PARAM_PROBE, -1);
	}

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

	// fill map
	int cur_numcols = CPXgetnumcols(env, lp);
	double x[cur_numcols];
	int solstat;
	double objval;
	status = CPXsolution(env, lp, &solstat, &objval, x, NULL, NULL, NULL);

	int mksp_index;
	status = CPXgetcolindex(env, lp, "x0", &mksp_index);
	assert(status == 0);

	printf("Status = %d, Solution status = %d, Objective value = %f Makespan x0 = %f\n", status,
		solstat, objval, x[mksp_index]);

	for(int cell = 0; cell < opt_atu; cell++) {
		for(int server = 1; server <= servers; server++) {
			char name[100];
			sprintf(name, "map(%d,%d)", cell, server);
			int index;
			status = CPXgetcolindex(env, lp, name, &index);
			assert(status == 0);
			map[cell][server] = x[index];
		}
	}

	status = CPXfreeprob(env, &lp);
	status = CPXcloseCPLEX(&env);

	#else

    glp_set_bfcp(prob, &bfcp);
    glp_intopt(prob, &iocp);

    glp_mpl_postsolve(tran, prob, GLP_SOL);

	int ncols = glp_get_num_cols(prob);
	for(int col = 1; col <= ncols; col++) {
		int server, cell;
		const char *colname = glp_get_col_name(prob, col);
		int i = sscanf(colname, "map[%d,%d]", &cell, &server);
		if (i == 2) {
			double v = glp_mip_col_val(prob, col);
			if (v == 1.0) {
				map[cell][server] = 1;
			}
			else
				printf("Error reading integer value from MIP column %s: %f\n", colname, v);
		}
	}
	#endif

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

	// write a file with the solution found
	write_solution_found(servers, opt_atu, map, lname, rname);

    glp_delete_prob(prob);
    glp_mpl_free_wksp(tran);
    //gmp_free_mem();
    glp_free_env();
}

void write_solution_found(int servers, int opt_atu, 
	int map[opt_atu][servers+1], char *lname, char *rname) {
	char filename[100];
	sprintf(filename, "sol/%s_%s_%d_%d.txt", lname, rname, servers, opt_atu);
	FILE *f = fopen(filename, "wb");
	for(int i=0; i < opt_atu; i++) {
		for(int s=1; s <= servers; s++) {
			fprintf(f, "%d:%d=%d\n", i, s, map[i][s]);
		}
	}
	fclose(f);
}

