
#include <stdlib.h>
#include <math.h>

double maxdbl(double data[], size_t start, size_t n) {
	double max = data[start];
	for(int i = start+1; i<n; i++) {
		if (max < data[i])
			max = data[i];
	}
	return max;
}

double stdevd_ex(void *data, size_t start, size_t n, double (*getv)(const void *, const int n)) {
	double mean=0.0, sum_deviation=0.0;
	int i;
	for(i=start; i<n; ++i)
		mean+=getv(data, i);
	mean=mean/(n-start);
	for(i=start; i<n; ++i)
		sum_deviation += (getv(data, i)-mean)*(getv(data, i)-mean);
	return sqrt(sum_deviation/(n-start));

}
