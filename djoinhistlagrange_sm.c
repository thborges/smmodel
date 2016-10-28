
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
	int *x0, int current_x0) {

    double result;

    // find makespan and communication cost
    double commcost = 0;
    int makespan[servers];
    memset(makespan, 0, sizeof makespan);
    for(int ij = 0; ij < pairs; ij++) {
    	for(int s = 0; s < servers; s++) {
    		if (x_ijk[s][ij] > 0) {
    			makespan[s] += opt_data[ij].pnts;
    			commcost += md->g * opt_data[ij].comm[s+1] + md->u[ij];
    		}
    	}
    }
    *x0 = 0;
    for(int s = 0; s < servers; s++) {
    	if (*x0 < makespan[s])
    		*x0 = makespan[s];
    }
	if (*x0 <= 0)
		printf("nenhum item colocado na knapsack\n");

	//printf("x0: %d\n", *x0);

    // find the constant term of the objective function
    double constantterm = 0;
    for(int ij = 0; ij < pairs; ij++)
    	constantterm += md->u[ij]; //HERE

    //return (md->f * (double)*x0) + commcost - constantterm;
	return (md->f * current_x0) + commcost - constantterm;
}

double lagrange_sm_get_heuristic_objective(dataset_histogram *hr, optimization_data_s *opt_data, int opt_atu,
    lagrange_model_data *md, int servers, double *x0_ub, double *x0_lb) {
    double x[servers+1];
    memset(x, 0, sizeof x);
    double netcost = 0.0;
    double cpucost = 0.0;

    for(int ij = 0; ij < opt_atu; ij++) {
        int xl = opt_data[ij].xl;
        int yl = opt_data[ij].yl;
        histogram_cell *cell = GET_HISTOGRAM_CELL(hr, xl, yl);
        netcost += opt_data[ij].comm[cell->place];
        x[cell->place] += opt_data[ij].pnts;
        cpucost += opt_data[ij].pnts;
    }

    // find x0: makespan
    *x0_ub = 0.0;
    for(int s = 1; s <= servers; s++) {
        if (x[s] > *x0_ub)
            *x0_ub = x[s];
    }

    *x0_lb = cpucost / (double)servers;

    return md->f * (*x0_ub) + md->g * netcost;
}

void lagrange_sm_optimize_hr(dataset_histogram *hr, int servers,
    optimization_data_s *opt_data, int pairs /*, multiway_histogram_estimate *agg_server*/) {

    lagrange_model_data md;
	md.u = g_new(double, pairs);

    printf("Pairs: %d\n", pairs);

    char *aux = getenv("LP_COMM");
    md.g = aux ? atof(aux) : 1.0;
    aux = getenv("LP_MKSP");
    md.f = aux ? atof(aux) : servers;
    printf("Tradeoff g: %f, f: %f\n", md.g, md.f);

    double x0_ub;
    double x0_lb;
    double Zheur = lagrange_sm_get_heuristic_objective(hr, opt_data, pairs, &md, servers, &x0_ub, &x0_lb);
	printf("Initial values from Greedy heuristic x0_ub: %f, x0_lb: %f\n", x0_ub, x0_lb);

    // 0-1 knapsack variables
    int profit[servers][pairs];
    int weight[pairs];
    int x_ijk[servers][pairs];

    for(int ij = 0; ij < pairs; ij++) {

    	// Initial values for md.u[]
    	md.u[ij] = - opt_data[ij].pnts;

        weight[ij] = opt_data[ij].pnts;

    	for(int s = 0; s < servers; s++) {
            // Initial profit values
            profit[s][ij] = - (/*md.f * opt_data[ij].pnts +*/ (md.g * opt_data[ij].comm[s+1]) + md.u[ij]);
    	}
    }

    const int stop = 10000;
    double lambda = 2.0;
    double t[stop];
    double Zd[stop];

    int notimproved = 0;
	int current_x0 = x0_lb;

	double best_zd = DBL_MAX;
	int best_x_ijk[servers][pairs];
	int zd_processed = 0;
	int best_processed = 0;

	printf("Z*: %f\n", Zheur);

    int k = 0;
    while (k < stop) {
	
		memset(x_ijk, 0, sizeof x_ijk);

		// run the knapsack algorithm for each server 
		for(int s = 0; s < servers; s++) {
			int z = minknap(pairs, profit[s], weight, x_ijk[s], current_x0);
		}

		int x0;
        Zd[k] = lagrange_sm_get_opt_value(&md, pairs, servers, x_ijk, opt_data, &x0, current_x0);

        // lambda set on k=3 iteration
        /*if (k == 0) {
            double gap = Zd[k] / Zheur;
            if (gap >= 0.1)
                lambda = 1.6;
            else if (gap >= 0.5)
                lambda = 0.75;
            else
                lambda = 0.5;
        }*/

        // count not processed ou double processed pairs
        int not_processed = 0;
        int double_processed = 0;
		int processed = 0;
		int weight_not_processed = 0;
        for(int ij = 0; ij < pairs; ij++) {
            int count = 0;
            for(int k = 0; k < servers; k++) {
                if (x_ijk[k][ij] > 0) {
                    count++;
                }
            }
            if (count == 0) {
				not_processed++;
				weight_not_processed += opt_data[ij].pnts;
			}
            else if (count > 1) double_processed++;
			else processed++;
        }

		/*if (not_processed < (int)(pairs*0.01) && double_processed <= 1) {
			int space_needed = current_x0 + (weight_not_processed/(double)servers);
			printf("Increased the knapsack from %d to %d\n", current_x0, space_needed);
			current_x0 = space_needed;
		}*/

		//if (x0 < current_x0)
		//	current_x0 = x0;

		// compute subgrad vector
        double subgrad[pairs];
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

		int choosed = 0;
		if (best_processed < processed || (best_processed == processed && Zd[k] < zd_processed)) {
			best_processed = processed;
			zd_processed = Zd[k];
			memcpy(best_x_ijk, x_ijk, sizeof x_ijk);
			choosed = 1;
		}

        if (best_zd > Zd[k]) {
            best_zd = Zd[k];
            notimproved = 0;
        }
        else
            notimproved++;

        if (notimproved >= 30) {
            lambda /= 2.0;
            notimproved = 0;
			if (lambda >= 1e-4)
				best_zd = DBL_MAX;
        }

        //double varobj = MAX(.04*Zd[k], (1.01*Zheur)-Zd[k]);
        double varobj = Zheur - Zd[k];
        t[k] = (lambda * varobj) / norm;

        //t[k] = (lambda * (Zheur - Zd[k])) / norm;

        for(int ij = 0; ij < pairs; ij++) {
            md.u[ij] = md.u[ij] + t[k] * subgrad[ij];

            // update knapsack profits
        	for(int s = 0; s < servers; s++) {
                profit[s][ij] = - (/*md.f * opt_data[ij].pnts +*/ (md.g * opt_data[ij].comm[s+1]) + md.u[ij]);
        	}
        }

        printf("\033[91m k=%4d, Zd: %'15.2f, lamb: %5.4f, T[k]: %'10.4f, norm: %'10.2f, x0: %'10d np: %5d, dp: %5d cp: %5d %c\n\033[0m", 
            k, Zd[k], lambda,
            t[k], norm, x0, 
            not_processed, double_processed, processed, choosed ? 'x' : ' ');

		// stop condition
        if ((fabs(t[k]) < 1e-4 && lambda < 1e-4) || norm < 1e-10) {

			/*if (not_processed > servers) {
				int space_needed = (weight_not_processed/(double)servers);
				space_needed = current_x0 + space_needed * 0.05;
				printf("Increased the knapsack from %d to %d\n", current_x0, space_needed);
				current_x0 = space_needed;
				lambda = 0.5;
			}
			else*/ {
				printf("Exit. t[k] = %f, norm = %f\n", t[k], norm);
    	        break;
			}
		}
        
        k++;
    }

    // debug md.u values    
    /*for(int ij = 0; ij < opt_atu; ij++) {
        printf("%d, %f\n", ij, md.u[ij]);
    }*/

	// round if needed
	printf("Found solution with %d/%d correctly processed items\n", best_processed, pairs);
   	lp_optimize_hr_round_decreasing_low_comm(hr, servers, opt_data, pairs, best_x_ijk);

	g_free(md.u);
}

