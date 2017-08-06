
#ifndef UTILS_H
#define UTILS_H

double maxdbl(double data[], size_t start, size_t n);
double stdevd_ex(void *data, size_t start, size_t n, double (*getv)(const void *, const int n));

#endif

