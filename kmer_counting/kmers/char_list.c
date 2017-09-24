/*
 *  char_list.c
 *
 *  Created on: Dec 12, 2011
 *
 *  Copyright (C) 2014  Ivan Gregor
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "char_list.h"

/**
 * Structure that represents the char list.
 * */
typedef struct {
	char * array;
	unsigned long maxSize;
	unsigned long currentSize;
} List;


/**
 * Creates a new char list.
 *
 * @param initSize initial size of the char list
 *
 * @return a structure that represents this char list or NULL
 * */
inline CList * clistCreate(unsigned long initSize){

	if (initSize < 1){
		return NULL;
	}

	List * list = (List *)malloc(sizeof(List));
	if (list == NULL){
		return NULL;
	}

	list->array = (char *)malloc((initSize + 1)*sizeof(char));
	if (list->array == NULL){
		free(list);
		return NULL;
	}
	//memset(list->array, 0, (initSize + 1)*sizeof(char));
	list->currentSize = 0;
	list->maxSize = initSize;
	return list;
}


/**
 * Gets the current number of chars stored in the array list.
 * */
inline unsigned long clistSize(CList * list){
		return (list != NULL)?(((List *)list)->currentSize):0;
}


/**
 * Adds a next char to the end of the char list.
 *
 * @return 0 if ok
 * */
inline int clistAdd(CList * list, char c){

	List * lst = (List *)list;
	if (lst == NULL){
		return 1;
	}

	if (lst->currentSize == lst->maxSize){

		char * newArray = (char *)malloc(((lst->maxSize)*2+1)*sizeof(char));
		if (newArray == NULL){
			return 1;
		}

        //memset((newArray + lst->maxSize),0,(lst->maxSize+1)*sizeof(char));
        memcpy(newArray,lst->array,lst->maxSize*sizeof(char));
		//int i;
		//for (i=0; i<lst->maxSize; i++){
		//	newArray[i] = lst->array[i];
		//}

        free(lst->array);
        lst->array = newArray;
        lst->maxSize = (lst->maxSize)*2;
	}
	lst->array[lst->currentSize] = c;
	++(lst->currentSize);

	return 0;
}


/**
 * Get the string represented by this char list.
 * */
inline char * clistGetString(CList * list){
	List * lst = (List *)list;
	if (lst == NULL || lst->array == NULL){
		return NULL;
	}

	lst->array[lst->currentSize] = '\0';
	return lst->array;
}


/**
 * Removes all elements from the char array.
 * */
inline void clistClear(CList *list){
	List * lst = (List *)list;
	//memset(lst->array, 0, (lst->currentSize + 1)*sizeof(char));
	lst->currentSize = 0;
}


/**
 * Deletes the list.
 * (free the memory)
 * */
inline void clistDestroy(CList *list){

	List * lst = (List *)list;
	if (lst == NULL){
		return;
	}
	if (lst->array != NULL){
		free(lst->array);
	}
	free(lst);
}

/**
 * Unit test
 * */
inline void clistUnitTest(){

	srand (time(NULL));

	unsigned long i,SIZE = 8300000;
	char randStr[SIZE+1];
	for (i=0; i<SIZE; i++){
		randStr[i] = (char)((rand() % 75) + 48);
	}
	randStr[SIZE] = '\0';
	CList * list = clistCreate(1);
	if (list == NULL){
		printf("clistUnitTest: Cannot create a list.\n");
		return;
	}
	for (i=0; i<SIZE; i++){
		if (clistAdd(list, randStr[i]) != 0){
			printf("clistUnitTest: Cannot add a char at %lu iter\n",i);
			return;
		}
	}
	char * str = clistGetString(list);
	if (str == NULL){
		printf("clistUnitTest: Can`t get the string\n");
		return;
	}
	if (strcmp(randStr, str) != 0){
		printf("clistUnitTest: Strings are not equal.\n");
		printf("str ref:   %s\n", randStr);
		printf("str clist: %s\n", str);
		return;
	}

	clistClear(list);
	clistDestroy(list);

	printf("clistUnitTest: success\n");
}

inline void clistSpeedTest(){
	srand (time(NULL));
	unsigned long i,j, SIZE = 10;
	char randStr[SIZE];
	for (i=0; i<SIZE; i++){
		randStr[i] = (char)((rand() % 26) + 65);
	}
	CList * list = clistCreate(100);
	for (i=0; i<250000; i++){
		for (j=0; j<1000; j++){
			if (clistAdd(list, randStr[i%SIZE]) != 0){
				printf("clistSpeedTest: can`t add char");
				return;
			}
		}
		clistClear(list);
	}

	clistDestroy(list);

	printf("clistSpeedTest: end\n");

}

