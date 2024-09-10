#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <data_types.h>


typedef int IsEqualFunc(void* a, void* b);


//Comparison functions ***********************************************************************************************
int IsEqualInt(void* a, void* b)
{
    return *(int*)a == *(int*)b;
}

int IsEqualShort(void* a, void* b)
{
    return *(short int*)a == *(short int*)b;
}

int IsEqualDouble(void* a, void* b)
{
    return (fabs(*(double*)a - *(double*)b) < 1e-10);
}

int IsEqualFloat(void* a, void* b)
{
    return (fabsf(*(float*)a - *(float*)b) < 1e-10f);
}

int IsEqualChar(void* a, void* b)
{
    return *(char*)a == *(char*)b;
}

bool DataTypes_IsEqual(enum AsynchTypes type, void* lhs, void* rhs)
{

    static IsEqualFunc * const IsEqual[ASYNCH_NUM_DATA_TYPES] = {
        &IsEqualChar,
        &IsEqualShort,
        &IsEqualInt,
        &IsEqualFloat,
        &IsEqualDouble
    };

    return IsEqual[type](lhs, rhs);
}