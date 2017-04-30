ACLOCAL_AMFLAGS = -I m4

opt=
if HAVE_DOXYGEN
	opt = docs
endif
SUBDIRS = src $(opt)
EXTRA_DIST = docs
dist-hook:
	rm $(distdir)/docs/Makefile.in
