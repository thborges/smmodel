/*
 * structs.h
 *
 *  Created on: 10/23/2016
 *      Author: Thiago Borges de Oliveira 
 */

#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdio.h>

#define GET_HISTOGRAM_CELL(h, x, y) &(h)->hcells[(x)*(h)->yqtd + (y)]
#define SET_IN_PLACE(var, place) (var |= (1<<(place-1)))
#define IS_IN_PLACE(var, place) (((var) >> (place-1)) & 1)
#define GET_ROD_ID(cell) ((cell)->xr + ((cell)->yr << 16))

struct Envelope { 
		struct {
			double MinX;
			double MinY;
			double MaxX;
			double MaxY;
		};
};

typedef struct Envelope Envelope;

typedef struct {
	double cardin;
	double points;
	unsigned place;
	unsigned copies;
} histogram_cell;

typedef	struct {
	Envelope mbr;
	int xqtd;
	int yqtd;
	double xsize;
	double ysize;
	double *xtics;
	double *ytics;
	histogram_cell *hcells;
} dataset_histogram;

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

#endif

