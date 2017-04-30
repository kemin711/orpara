ACLOCAL_AMFLAGS = -I m4

opt=
if HAVE_DOXYGEN
	opt = docs
endif
SUBDIRS = src $(opt)

EXTRA_DIST = docs
dist-hook:
	rm -f $(distdir)/docs/doxyfile.stamp
	rm -f $(distdir)/docs/Makefile
	rm -rf $(distdir)/docs/html
	rm -rf $(distdir)/docs/latex
	rm -rf $(distdir)/docs/man
