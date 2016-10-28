
/*
 * main.c
 *
 *  Created on: 10/23/2016
 *      Author: Thiago Borges de Oliveira 
 */

#ifndef CPPHASH_H
#define CPPHASH_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void* CppMap;

CppMap map_create();
void map_destroy(CppMap *m);
void *map_get(CppMap *m, int key);
void map_put(CppMap *m, int key, CppMap *value);

#ifdef __cplusplus
}
#endif

#endif

