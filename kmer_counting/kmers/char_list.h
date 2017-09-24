/*
 *  char_list.h
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

#ifndef CHAR_LIST_H_
#define CHAR_LIST_H_

/**
 * Structure that represents the char list.
 * */
typedef void CList;

CList * clistCreate(unsigned long initSize);

unsigned long clistSize(CList * list);

int clistAdd(CList * list, char c);

char * clistGetString(CList * list);

void clistClear(CList *list);

void clistDestroy(CList * list);

void clistUnitTest();

void clistSpeedTest();

#endif /* CHAR_LIST_H_ */
