/*
 * structs.h
 *
 *  Created on: 10/23/2016
 *      Author: Thiago Borges de Oliveira 
 */

#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>
#include <ogr_api.h>
#include <stdbool.h>
#include <assert.h>
#include <emmintrin.h>
#include "glibwrap.h"

#define GET_GRID_HISTOGRAM_CELL(h, x, y) (&((h)->hcells[(x)*(h)->yqtd + (y)]))
#define SET_IN_PLACE(var, place) (var |= (1<<(place-1)))
#define IS_IN_PLACE(var, place) (((var) >> (place-1)) & 1)
#define GET_ROD_ID(cell) ((cell)->xr + ((cell)->yr << 16))

#ifdef __MACH__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(void *thunk, const void *x, const void *y)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, thunk, compar)
#elif __linux__
#define decl_qsort_p_cmp(fname, x, y, thunk) int fname(const void *x, const void *y, void *thunk)
#define qsort_p(base, nmemb, size, compar, thunk) qsort_r(base, nmemb, size, compar, thunk)
#endif

#define g_new(ptype, count) (ptype*)malloc(sizeof(ptype) * (count))
#define g_new0(ptype, count) (ptype*)calloc(count, sizeof(ptype))
#define g_renew(ptype, pointer, count) (ptype*)realloc(pointer, sizeof(ptype) * (count))


struct Envelope { 
	union {
		struct {
			__m128d Min;
			__m128d Max;
		};
		struct {
			double MinX;
			double MinY;
			double MaxX;
			double MaxY;
		};
	};
};


typedef struct Envelope Envelope;

typedef struct {
	unsigned short x;
	unsigned short y;
	double cardin;
	double points;
	unsigned place;
	unsigned copies;
	void *extra_data;
	double avglength_x;
	double avglength_y;
	double objcount;
	double prep_cost;
} histogram_cell;

enum HistogramType {
	HIST_GRID,
	HIST_MINSKEW,
};

typedef struct {
	enum HistogramType htype;
	Envelope mbr;
	char extra_data[1];
} dataset_histogram_persist;

typedef struct dataset_head_ dataset;
typedef struct dataset_histogram dataset_histogram;

typedef	struct dataset_histogram {
	enum HistogramType htype;
	OGRwkbGeometryType geom_type;
	void *extra_data;

	int (*free)(dataset *ds);
	void (*print_geojson)(dataset *ds);
	void (*distribute)(dataset *ds);
	GList* (*intersects)(dataset *ds, Envelope ev);
	void (*clone_histogram)(dataset *dh, dataset_histogram *result);
	Envelope (*get_cell_envelope)(dataset_histogram *dh, histogram_cell *cell);
	histogram_cell *(*get_cell)(dataset_histogram *dh, unsigned short x, unsigned short y);
	histogram_cell *(*get_reference_cell)(dataset *ds, double x, double y);

	dataset_histogram_persist *(*get_data)(dataset *dh, int *size);
	void (*set_data)(dataset *dh, dataset_histogram_persist *data);
} dataset_histogram;

typedef struct {
	int xqtd;
	int yqtd;
	double xsize;
	double ysize;
	double *xtics;
	double *ytics;
	histogram_cell *hcells;
} grid_histogram_data;

typedef struct dataset_head_ {
	struct {
		Envelope mbr;
		struct dataset_histogram hist;
	} metadata;
} dataset;

typedef struct {
	int xr;
	int yr;
	histogram_cell *cell;
} right_opt_data;

typedef struct {
	int xl;
	int yl;
	histogram_cell *lcell;
	double pnts;
	double *comm;
	size_t rcells_size;
	right_opt_data *rcells;
} optimization_data_s;

typedef struct {
	double to_pnts;
	double io_objs;
	double io_pnts;
	double io_save;
	double inters;
} multiway_histogram_estimate;

void print_instance_and_solution_fo_file(optimization_data_s *opt_data, int opt_atu, int servers, int map[opt_atu][servers+1]);

#endif

