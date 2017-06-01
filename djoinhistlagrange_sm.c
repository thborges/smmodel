
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <glibwrap.h>
#include <round.h>
#include <structs.h>

typedef struct {
	double g;
	double f;

	double *u;
} lagrange_model_data;

long minknap(int n, int *p, int *w, int *x, int c);

double lagrange_sm_get_opt_value(lagrange_model_data *md,
	int pairs, int servers,
	int x_ijk[servers][pairs], 
	optimization_data_s *opt_data,
	double *x0) {

	// find makespan and communication cost
	int count = 0;
	double commcost = 0;
	double makespan[servers];
	memset(makespan, 0, sizeof makespan);
	for(int ij = 0; ij < pairs; ij++) {
		for(int s = 0; s < servers; s++) {
			if (x_ijk[s][ij] > 0) {
				makespan[s] += opt_data[ij].pnts;
				commcost += md->g * opt_data[ij].comm[s+1] + md->u[ij];
				count++;
			}
		}
	}
	*x0 = 0;
	for(int s = 0; s < servers; s++) {
		if (*x0 < makespan[s])
			*x0 = makespan[s];
	}
	if (*x0 <= 0)
		printf("nenhum item colocado na knapsack: %f, %d\n", *x0, count);

	// find the constant term of the objective function
	double constantterm = 0;
	for(int ij = 0; ij < pairs; ij++)
		constantterm += md->u[ij]; //HERE

	return md->f * (*x0) + commcost - constantterm;
}

double get_sm_objective(dataset_histogram *hr, optimization_data_s *opt_data, int opt_atu,
	double f, double g, int servers, double multiplier, double *x0_ub, double *x0_lb, double *fmkspan, double *fcomm) {
	double x[servers+1];
	memset(x, 0, sizeof x);
	double xm[servers+1];
	memset(xm, 0, sizeof xm);
	double netcost = 0.0;
	double cpucost = 0.0;

	for(int ij = 0; ij < opt_atu; ij++) {
		int xl = opt_data[ij].xl;
		int yl = opt_data[ij].yl;
		histogram_cell *cell = hr->get_cell(hr, xl, yl);
		netcost += opt_data[ij].comm[cell->place];
		x[cell->place] += opt_data[ij].pnts;
		xm[cell->place] += ceil(opt_data[ij].pnts * multiplier);
		cpucost += ceil(opt_data[ij].pnts * multiplier);
	}

	// find x0: makespan
	if (x0_ub) {
		*x0_ub = 0.0;
		for(int s = 1; s <= servers; s++) {
			if (xm[s] > *x0_ub)
				*x0_ub = xm[s];
		}
	}
	double real_mkspan = 0.0;
	for(int s = 1; s <= servers; s++) {
		if (x[s] > real_mkspan)
			real_mkspan = x[s];
	}

	if (x0_lb)
		*x0_lb = cpucost / (double)servers;

	if (fcomm)
		*fcomm = netcost;
	if (fmkspan)
		*fmkspan = real_mkspan;

	return f * real_mkspan + g * netcost;
}

int int_decreasing_compare(const void* p1, const void* p2)
{ 
   int i1 = *(int*) p1;
   int i2 = *(int*) p2;
   if (i1 < i2) return 1;
   else if (i1 == i2) return 0;
   else return -1;
   /* or simply: return i1 - i2; */
}

double find_knapsack_bestf_capacity(optimization_data_s *opt_data, int pairs, int servers, double multiplier) {
	int points[pairs];
	for(int i = 0; i < pairs; i++)
		points[i] = ceil(opt_data[i].pnts * multiplier);
	qsort(points, pairs, sizeof(int), int_decreasing_compare);	

	double serverload[servers];
	memset(serverload, 0, sizeof serverload);
	for(int i = 0; i < pairs; i++) {
		int min = 0;
		double minvalue = serverload[min];
		for(int s = 1; s < servers; s++) {
			if (minvalue > serverload[s]) {
				min = s;
				minvalue = serverload[s];
			}
		}
		serverload[min] += points[i];
	}
	double max_load = 0;
	for(int s = 0; s < servers; s++) {
		if (max_load < serverload[s])
			max_load = serverload[s];
	}
	return max_load;
}

void lagrange_sm_optimize_hr(dataset_histogram *hr, int servers,
	optimization_data_s *opt_data, int pairs, multiway_histogram_estimate *agg_server,
	double dualvalues[pairs], double f, double *return_mkspan, double *return_comm) {

	lagrange_model_data md;
	md.u = g_new(double, pairs);

	//printf("Pairs: %d\n", pairs);

	char *aux = getenv("LP_COMM");
	md.g = aux ? atof(aux) : 1.0;
	aux = getenv("LP_MKSP");
	md.f = aux ? atof(aux) : f;
	//printf("Tradeoff g: %f, f: %f\n", md.g, md.f);

	// find a multiplier to reduce knapsack capacity
	double multiplier = 1;
	double knapsack_bestf = find_knapsack_bestf_capacity(opt_data, pairs, servers, 1);
	double auxmax = knapsack_bestf;
	while (auxmax > 1e6) {
		multiplier /= 10;
		auxmax /= 10;
	}
	multiplier = 1;
	//printf("Knapsack bestf %f: multiplier %e\n", knapsack_bestf, multiplier);

	double heur_x0;
	double x0_lb;
	double Zheur = get_sm_objective(hr, opt_data, pairs, md.f, md.g, servers, multiplier, &heur_x0, &x0_lb, NULL, NULL);
	//printf("Heuristic Z*: %'.2f x0: %f x0_lb: %f\n", Zheur, heur_x0, x0_lb);

	// 0-1 knapsack variables
	double best_u[pairs];
	double best_subgrad[pairs];
	int *weight = g_new(int, pairs);
	int (*profit)[pairs] = malloc(servers * sizeof *profit);
	int (*x_ijk)[pairs] = malloc(servers * sizeof *x_ijk);
	int (*best_x_ijk)[pairs] = malloc(servers * sizeof *best_x_ijk);
	memset(best_x_ijk, 0, sizeof(int)*servers*pairs);

	for(int ij = 0; ij < pairs; ij++) {
		// Initial values for md.u[]
		double min = DBL_MAX;
		for(int s = 0; s < servers; s++) {
			if (min > opt_data[ij].comm[s+1])
				min = opt_data[ij].comm[s+1];
		}
		//md.u[ij] = -((md.f/servers) * opt_data[ij].pnts + min);
		//md.u[ij] = - min;
		//md.u[ij] = ((md.f/servers) * opt_data[ij].pnts + min);
		//md.u[ij] = min;
		//md.u[ij] = - (x0_lb + min);
		md.u[ij] = 0;
		/*if (dualvalues) {
			md.u[ij] = - dualvalues[ij];
			printf("dual %i %f\n", ij, dualvalues[ij]);
		}*/
		
		//printf("(%f\t%f)\t", opt_data[ij].pnts, min);


		weight[ij] = ceil(opt_data[ij].pnts * multiplier);

		for(int s = 0; s < servers; s++) {
			// Initial profit values
			profit[s][ij] = - (/*md.f * opt_data[ij].pnts +*/ (md.g * opt_data[ij].comm[s+1]) + md.u[ij]);
		}

		// set the know solution, for the case no better solution is found
		histogram_cell *cell = hr->get_cell(hr, opt_data[ij].xl, opt_data[ij].yl);
		best_x_ijk[cell->place-1][ij] = 1;
	}

	const int stop = 10000;
	double lambda = 2.0;
	double lambda_reduce = 2.0;
	double tk;
	double Zdk;
	int notimproved = 0;
	double best_zd = -DBL_MAX;
	int zd_processed = 0;
	int best_processed = 0;

	int knapsack_lb = x0_lb;
	knapsack_bestf = find_knapsack_bestf_capacity(opt_data, pairs, servers, multiplier);
	double knapsack_gap = (knapsack_bestf - knapsack_lb) / knapsack_lb;
	//printf("Knapsack lb, bestf, gap: %f, %f, %f\n", x0_lb, knapsack_bestf, knapsack_gap);
	if (knapsack_gap < 0.01) {
		lambda_reduce = 1.5;
		knapsack_bestf *= 1.005;
		//printf("Knapsack increased to %f, due knapsack_gap < 0.01\n", knapsack_bestf);
	}

	double knapsack_atu = knapsack_bestf;

	double sol_gap = 0;
	int stable_sol_gap = 0;

	int solutions = 0;
	int decrease_lambda = 0;
	int k = 0;
	while (k < stop) {
		// run the knapsack algorithm for each server 
		for(int s = 0; s < servers; s++) {
			int z = minknap(pairs, profit[s], weight, x_ijk[s], knapsack_atu);
			//if (z == 0)
			//	minknap(pairs, profit[s], weight, x_ijk[s], 2*knapsack_atu);
		}

		double x0;
		Zdk = lagrange_sm_get_opt_value(&md, pairs, servers, x_ijk, opt_data, &x0);

		if (Zdk >= Zheur) {
			if (knapsack_atu > knapsack_lb) {
				knapsack_atu = MAX(knapsack_lb, knapsack_atu*.99);
				//printf("Knapsack reduced to %f. Zd %f, Zheur %f\n", knapsack_atu, Zdk, Zheur);
				if (knapsack_atu > knapsack_lb)
					continue;
			}
		}
		
		if (stable_sol_gap >= 0) {
			double gap = fabs(Zdk - Zheur) / Zheur;
			if (fabs(sol_gap - gap) < 0.01)
				stable_sol_gap++;
			else
				stable_sol_gap = 0;
			sol_gap = gap;
		}
		if (stable_sol_gap > 30) {
			//printf("Stable gap detected. Solution gap %f!\n", sol_gap);
			if (sol_gap > 0.10) {
				/*knapsack_atu *= (1.0+sol_gap);
				printf("Knapsack increased due high gap %f: %f\n", sol_gap, knapsack_atu);
				lambda = 2.0;*/
			}
			stable_sol_gap = -1; // disable
		}
		
		// lambda set on first iteration
		/*if (k == 0) {
			double gap = fabs(Zdk - Zheur) / Zheur;
			lambda = MIN(2, gap * 10.0);
		}*/

		// count not processed ou double processed pairs
		int not_processed = 0;
		int double_processed = 0;
		int processed = 0;

		for(int ij = 0; ij < pairs; ij++) {
			int count = 0;
			for(int k = 0; k < servers; k++) {
				if (x_ijk[k][ij] > 0)
					count++;
			}
			if (count == 0) {
				not_processed++;
			} else if (count >= 1) {
				processed++;
				if (count > 1)
					double_processed++;
			}
		}

		// compute subgrad vector
		double subgrad[pairs];
		#pragma omp parallel for
		for(int ij = 0; ij < pairs; ij++) {
			subgrad[ij] = 0.0;
			for(int k = 0; k < servers; k++) {
				subgrad[ij] += x_ijk[k][ij];
			}
			subgrad[ij] = subgrad[ij] - 1; //HERE
		}
		
		// compute norm
		double norm = 0;
		for(int ij = 0; ij < pairs; ij++)
			norm += subgrad[ij] * subgrad[ij];
		//norm = sqrt(norm);

		char choosed = ' ';
		if (best_processed < processed || (best_processed == processed && Zdk > zd_processed)) {
			choosed = '!';
			solutions++;

			// number of processed items improved. try to get a new upper bound
			lp_optimize_hr_round_decreasing_low_comm(hr, servers, opt_data, pairs, x_ijk, agg_server, true, md.f, md.g);
			double NewZheur = get_sm_objective(hr, opt_data, pairs, md.f, md.g, servers, multiplier, &heur_x0, NULL, NULL, NULL);
			if (NewZheur < Zheur) {
				choosed = 'x';
				best_processed = processed;
				zd_processed = Zdk;
				memcpy(best_x_ijk, x_ijk, sizeof(int)*servers*pairs);
				//printf("New feasible solution (UB): %'.0f, x0 %'f\n", NewZheur, heur_x0);
				Zheur = NewZheur;
				memcpy(best_u, md.u, sizeof(double)*pairs);
				memcpy(best_subgrad, subgrad, sizeof(double)*pairs);
			}
		}

		char improved = ' ';
		if (best_zd < Zdk) {
			improved = 'x';
			best_zd = Zdk;
			notimproved = 0;
		}
		else
			notimproved++;

		if (notimproved >= 50) {
			lambda /= lambda_reduce;
			decrease_lambda++;
			notimproved = 0;
		}

		//double varobj = MAX(.04*Zdk, (1.01*Zheur)-Zdk);
		//double varobj = 0.01 * (Zheur - Zdk);
		double varobj = (Zheur - Zdk);
		tk = (lambda * varobj) / norm;

		#pragma omp parallel for
		for(int ij = 0; ij < pairs; ij++) {
			md.u[ij] = md.u[ij] + tk * subgrad[ij];

			// update knapsack profits
			for(int s = 0; s < servers; s++) {
				profit[s][ij] = - (/*md.f * opt_data[ij].pnts +*/ (md.g * opt_data[ij].comm[s+1]) + md.u[ij]);
			}
		}

		if (choosed != ' ' || k%100==0) {
			fprintf(stderr, "\033[91m k=%4d, Zd: %'15.2f%c lamb: %5.4f, T[k]: %'15.4f, norm: %'10.2f, x0: %'10.0f np: %5d, dp: %5d v: %10.0f cp: %5d %c\r\033[0m", 
				k, Zdk, improved, lambda,
				tk, norm, x0, 
				not_processed, double_processed, knapsack_atu, processed, choosed);
		}

		// stop condition
		if ((fabs(tk) <= 1e-3 && lambda < 1e-4) || norm < 1e-10) {
			//printf("Exit. t[k] = %f, norm = %f\n", tk, norm);
			break;
		}
		
		k++;
	}

	printf("\nLast iteration k: %d\n", k);

	// debug md.u values	
	/*printf("-------\n");
	for(int ij = 0; ij < pairs; ij++) {
		printf("%f %f\n", best_u[ij], best_subgrad[ij]);
	}*/

	// round if needed
	double final_mkspan, final_comm;
	set_cell_place_from_partial_x(hr, servers, pairs, best_x_ijk, opt_data);
	Zheur = get_sm_objective(hr, opt_data, pairs, md.f, md.g, servers, multiplier, NULL, NULL, &final_mkspan, &final_comm);
	printf("Before LR rounding\nZ\tMkspan\tComm\n%.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);
	printf("Using solution with %d/%d correctly processed items\n", best_processed, pairs);

	lp_optimize_hr_round_decreasing_low_comm(hr, servers, opt_data, pairs, best_x_ijk, agg_server, true, md.f, md.g);
	Zheur = get_sm_objective(hr, opt_data, pairs, md.f, md.g, servers, multiplier, NULL, NULL, &final_mkspan, &final_comm);
	printf("After LR rounding\nZ\tMkspan\tComm\n%.2f\t%.2f\t%.2f\n", Zheur, final_mkspan, final_comm);

	*return_mkspan = final_mkspan;
	*return_comm = final_comm;

	/* print result for Prof. Les instances */
	/*for(int s = 1; s <= servers; s++) {
		for(int ij = 0; ij < pairs; ij++) {
			int xl = ij;
			int yl = 0;
			histogram_cell *cell = hr->get_cell(hr, xl, yl);
			printf("%d\t", cell->place == s ? 1 : 0);
		}
		printf("\n");
	}*/

	g_free(md.u);
	g_free(profit);
	g_free(weight);
	g_free(x_ijk);
	g_free(best_x_ijk);
}

