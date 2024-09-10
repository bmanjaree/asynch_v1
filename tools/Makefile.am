if USE_POSTGRESQL

bin_PROGRAMS = asynch_createtables asynch_deletetables

asynch_createtables_SOURCES = createtables.c
asynch_createtables_LDADD = $(POSTGRESQL_LIBS)

asynch_deletetables_SOURCES = deletetables.c
asynch_deletetables_LDADD = $(POSTGRESQL_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS)

endif
