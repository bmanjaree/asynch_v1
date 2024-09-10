asynchdir = $(pkgpythondir)
asynch_PYTHON = asynch_interface.py

asynch_LTLIBRARIES = _asynchpy.la

_asynchpy_la_SOURCES = asynch_interface_py.c
_asynchpy_la_LDFLAGS = -module -avoid-version -shared #-export-symbols-regex initmyext
_asynchpy_la_LIBADD = $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS) $(PYTHON_EXTRA_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS) $(HDF5_LDFLAGS) $(SSH2_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS) $(HDF5_CFLAGS) $(SSH2_CPPFLAGS) $(PYTHON_CPPFLAGS)

EXTRA_DIST = asynch_interface_py.h

