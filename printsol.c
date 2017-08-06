#include "structs.h"

void print_instance_and_solution_fo_file(optimization_data_s *opt_data, int opt_atu, int servers, int x[opt_atu][servers+1]) {
	FILE *mipsolf = fopen("mipsol.txt", "w");

	fprintf(mipsolf, "Computed costs c(i,j)\t");
	for(int j = 0; j < opt_atu; j++)
		fprintf(mipsolf, "j=%d\t", j+1);
	fprintf(mipsolf, "\n");

	for(int server = 1; server <= servers; server++) {
		fprintf(mipsolf, "i=%d\t", server);
		for(int cell = 0; cell < opt_atu; cell++)
			fprintf(mipsolf, "%f\t", opt_data[cell].comm[server]);
		fprintf(mipsolf, "\n");
	}

	fprintf(mipsolf, "\n\t");
	for(int j = 0; j < opt_atu; j++)
		fprintf(mipsolf, "j=%d\t", j+1);
	fprintf(mipsolf, "\n");

	fprintf(mipsolf, "w(j)\t");
	for(int cell = 0; cell < opt_atu; cell++)
		fprintf(mipsolf, "%f\t", opt_data[cell].pnts);
	fprintf(mipsolf, "\n");

	if (x) {
		fprintf(mipsolf, "x(i,j)\t");
		for(int j = 0; j < opt_atu; j++)
			fprintf(mipsolf, "j=%d\t", j+1);
		fprintf(mipsolf, "\n");

		for(int server = 1; server <= servers; server++) {
			fprintf(mipsolf, "i=%d\t", server);
			for(int cell = 0; cell < opt_atu; cell++)
				fprintf(mipsolf, "%d\t", x[cell][server]);
			fprintf(mipsolf, "\n");
		}
	}

	fclose(mipsolf);

}
