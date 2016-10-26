

EXTRA_DIST = defs ss.chr17test ss.her2 \
             align.test.ok coords1.test.ok \
             setup.genomecomp.ok setup.ref123positions.ok \
             map.test.ok \
             fa.iittest iittest.iit.ok iit_get.out.ok

CLEANFILES = align.test.out \
             coords.chr17test \
             gmap_setup gmapindex chr17test/chr17test.* Makefile.chr17test \
             map.test.out \
             iittest.iit iit_get.out iit_dump.out

if MAINTAINER
TESTS = align.test coords1.test setup1.test iit.test \
        bioinfo.splicing.test
else
TESTS = align.test coords1.test setup1.test iit.test 
endif

TESTS_ENVIRONMENT = top_builddir=$(top_builddir)

testsubdir = testSubDir

distclean-local:
	rm -rf $(testsubdir)

okay:
	../src/gmap -A -g ss.chr17test ss.her2 > align.test.ok
#	./fa_coords -o coords.chr17test ss.chr17test
#	./gmap_setup -D . -d chr17test ss.chr17test && make -f Makefile.chr17test coords && make -f Makefile.chr17test gmapdb
	./gmap_build -D . -d chr17test ss.chr17test
	../src/gmap -D . -d chr17test ss.her2 > map.test.ok
	mv -f coords.chr17test coords1.test.ok
	mv -f chr17test.genomecomp setup.genomecomp.ok
	mv -f chr17test.ref123positions setup.ref123positions.ok
	../src/iit_store -o iittest fa.iittest
	../src/iit_get iittest 20000 > iit_get.out.ok
	mv iittest.iit iittest.iit.ok
