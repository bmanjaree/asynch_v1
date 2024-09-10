#if !defined(CONFIG_H)
#define CONFIG_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <structs_fwd.h>


/// Reads in the data from a .gbl file. All data will be global data for the entire river system.
///
/// \param globalfilename String with the filename of the .gbl file.
/// \param errors (set by this method) Will contain the error data for the entire river system, if the error data is global.
/// \param forcings 
/// \param db_connections NULL pointer that will be set to an SQL database, if needed.
/// \param rkdfilename (set by this method) Will be the filename of the .rkd file, if the error data is not global.
/// \param model Pointer to a custom model
/// \return Configuration data read in from the file globalfilename.
GlobalVars* Read_Global_Data(
    char *globalfilename,
    ErrorData *errors,
    Forcing *forcings,
    ConnData *db_connections,
    char *rkdfilename,
    AsynchModel const *model,
    void *external);


#endif //CONFIG_H

