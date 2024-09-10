#if !defined (ASYNCH_DATE_MANIP_H)
#define ASYNCH_DATE_MANIP_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <time.h>


/// Copy a struct tm
void tmcpy(struct tm *source,struct tm *dest);

/// Returns the number of days in a month
int days_in_month(int mon,int year);


#endif //ASYNCH_DATE_MANIP_H

