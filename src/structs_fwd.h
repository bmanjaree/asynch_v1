#if !defined(STRUCTS_FWD_H)
#define STRUCTS_FWD_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// Forward definitions
typedef struct ErrorData ErrorData;
typedef struct GlobalVars GlobalVars;

typedef struct Lookup Lookup;

typedef struct DataPoint DataPoint;
typedef struct TimeSerie TimeSerie;

typedef struct Link Link;
typedef struct TransData TransData;
typedef struct Workspace Workspace;
typedef struct ConnData ConnData;

typedef struct Forcing Forcing;
typedef struct ForcingData ForcingData;

typedef struct QVSData QVSData;

typedef struct OutputFunc OutputFunc;

typedef struct RKMethod RKMethod;
typedef struct RKSolutionNode RKSolutionNode;
typedef struct RKSolutionList RKSolutionList;

typedef struct AsynchModel AsynchModel;
typedef struct AsynchSolver AsynchSolver;

#endif //STRUCTS_FWD_H
