#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdbool.h>


enum AsynchTypes
{
    ASYNCH_BAD_TYPE = -1,
    ASYNCH_CHAR,
    ASYNCH_SHORT,
    ASYNCH_INT,
    ASYNCH_FLOAT,
    ASYNCH_DOUBLE,
    ASYNCH_NUM_DATA_TYPES
};

bool DataTypes_IsEqual(enum AsynchTypes type, void* lhs , void* rhs);

#endif

