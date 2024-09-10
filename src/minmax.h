#if !defined(ASYNCH_MINMAX_H)
#define ASYNCH_MINMAX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(min)
static __inline double min(double a, double  b) {
    return a > b ? b : a;
}
#endif

#if !defined(max)
static __inline double max(double a, double b) {
    return a > b ? a : b;
}
#endif


#endif //ASYNCH_MINMAX_H