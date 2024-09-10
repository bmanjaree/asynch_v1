noinst_LIBRARIES = libforecaster.a libassim.a

libforecaster_a_SOURCES = forecaster_methods.c
libassim_a_SOURCES = assim_ls_methods.c assim_models.c

bin_PROGRAMS = assim #assim_persis assimpersis_end assim_least_square_linear

include_HEADERS = forecaster_methods.h assim_models.h assim_ls_methods.h

assim_SOURCES = assim.c optparse.c
assim_LDADD = libassim.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS) $(METIS_LIBS) $(PETSC_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS) $(HDF5_LDFLAGS) $(METIS_LDFLAGS) $(PETSC_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS) $(HDF5_CFLAGS) $(METIS_CFLAGS) $(PETSC_CFLAGS)
