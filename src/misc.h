#ifndef MISC_H
#define MISC_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <system.h>
#include <stdlib.h>
#include <stdio.h>

void CreateGraph(Link** sys,unsigned int N);
void CalcHortonOrder(Link** sys,unsigned int N,unsigned int* order,unsigned short int* complete);
void CreateStrComplete(Link** sys,unsigned int N);

#endif //MISC_H
