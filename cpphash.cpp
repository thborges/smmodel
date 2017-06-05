
/*
 * main.c
 *
 *  Created on: 10/23/2016
 *      Author: Thiago Borges de Oliveira 
 */

#include <cstddef>
#include <map>
#include "cpphash.h"

using namespace std;

CppMap map_create() {
	map<int,void*> *result = new map<int,void*>();
	return result;
}

void map_destroy(CppMap *m) {
	map<int,void*> *o = (map<int,void*>*)m;
	delete o;
	m = NULL;
}

void *map_get(CppMap *m, int key) {
	map<int,void*> *o = (map<int,void*>*)m;
	return (*o)[key];
}

void map_put(CppMap *m, int key, CppMap *value) {
	map<int,void*> *o = (map<int,void*>*)m;
	(*o)[key] = value;
}

