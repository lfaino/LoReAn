# Copied from ../util so "gmap_setup -B src" works during "make check"
CLEANFILES = fa_coords gmap_process

# Include master.c and master.h, so they get included in distribution
EXTRA_DIST = mpidebug.c mpidebug.h master.c master.h

# -DUTILITYP=1 needed for iit_store, iit_dump, iit_get, gmapindex, and
# get-genome so they can handle both small and large genomes at run
# time

bin_PROGRAMS = gmap gmapl get-genome gmapindex iit_store iit_get iit_dump \
               gsnap gsnapl uniqscan uniqscanl snpindex cmetindex atoiindex \
               sam_sort


GMAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 gbuffer.c gbuffer.h genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 genome-write.c genome-write.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h uinttable.c uinttable.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 diagdef.h diag.c diag.h diagpool.c diagpool.h \
 cmet.c cmet.h atoi.c atoi.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h \
 scores.h intron.c intron.h maxent.c maxent.h maxent_hr.c maxent_hr.h samflags.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h cellpool.c cellpool.h stage2.c stage2.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 boyer-moore.c boyer-moore.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h sense.h fastlog.h stage3.c stage3.h \
 request.c request.h result.c result.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h outbuffer.c outbuffer.h \
 chimera.c chimera.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gmap_CC = $(PTHREAD_CC)
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_gmap_SOURCES = $(GMAP_FILES)


GMAPL_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 gbuffer.c gbuffer.h genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 genome-write.c genome-write.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h uinttable.c uinttable.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 diagdef.h diag.c diag.h diagpool.c diagpool.h \
 cmet.c cmet.h atoi.c atoi.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h \
 scores.h intron.c intron.h maxent.c maxent.h maxent_hr.c maxent_hr.h samflags.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h cellpool.c cellpool.h stage2.c stage2.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 boyer-moore.c boyer-moore.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h sense.h fastlog.h stage3.c stage3.h \
 request.c request.h result.c result.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h outbuffer.c outbuffer.h \
 chimera.c chimera.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gmapl_CC = $(PTHREAD_CC)
gmapl_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1
gmapl_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_gmapl_SOURCES = $(GMAPL_FILES)



GSNAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h samflags.h samprint.c samprint.h \
 mapq.c mapq.h shortread.c shortread.h substring.c substring.h junction.c junction.h stage3hr.c stage3hr.h \
 spanningelt.c spanningelt.h cmet.c cmet.h atoi.c atoi.h \
 comp.h maxent.c maxent.h pairdef.h pair.c pair.h pairpool.c pairpool.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h cellpool.c cellpool.h stage2.c stage2.h \
 intron.c intron.h boyer-moore.c boyer-moore.h changepoint.c changepoint.h pbinom.c pbinom.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 gbuffer.c gbuffer.h translation.c translation.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 chimera.c chimera.h sense.h fastlog.h stage3.c stage3.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 splice.c splice.h indel.c indel.h bitpack64-access.c bitpack64-access.h \
 bytecoding.c bytecoding.h univdiagdef.h univdiag.c univdiag.h sarray-read.c sarray-read.h \
 stage1hr.c stage1hr.h \
 request.c request.h resulthr.c resulthr.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h outbuffer.c outbuffer.h \
 datadir.c datadir.h mode.h \
 getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags
# Previously included -lrt for shm_open, but we are not calling that

gsnap_CC = $(PTHREAD_CC)
gsnap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1
gsnap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_gsnap_SOURCES = $(GSNAP_FILES)


GSNAPL_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h uint8list.c uint8list.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h samflags.h samprint.c samprint.h \
 mapq.c mapq.h shortread.c shortread.h substring.c substring.h junction.c junction.h stage3hr.c stage3hr.h \
 spanningelt.c spanningelt.h cmet.c cmet.h atoi.c atoi.h \
 comp.h maxent.c maxent.h pairdef.h pair.c pair.h pairpool.c pairpool.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h cellpool.c cellpool.h stage2.c stage2.h \
 intron.c intron.h boyer-moore.c boyer-moore.h changepoint.c changepoint.h pbinom.c pbinom.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 gbuffer.c gbuffer.h translation.c translation.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 chimera.c chimera.h sense.h fastlog.h stage3.c stage3.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 splice.c splice.h indel.c indel.h bitpack64-access.c bitpack64-access.h \
 stage1hr.c stage1hr.h \
 request.c request.h resulthr.c resulthr.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h outbuffer.c outbuffer.h \
 datadir.c datadir.h mode.h \
 getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gsnapl_CC = $(PTHREAD_CC)
gsnapl_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1 -DLARGE_GENOMES=1
gsnapl_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_gsnapl_SOURCES = $(GSNAPL_FILES)


UNIQSCAN_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h \
 mapq.c mapq.h shortread.c shortread.h substring.c substring.h junction.c junction.h stage3hr.c stage3hr.h \
 spanningelt.c spanningelt.h cmet.c cmet.h atoi.c atoi.h \
 comp.h maxent.c maxent.h pairdef.h pair.c pair.h pairpool.c pairpool.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h cellpool.c cellpool.h stage2.c stage2.h \
 intron.c intron.h boyer-moore.c boyer-moore.h changepoint.c changepoint.h pbinom.c pbinom.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 chimera.c chimera.h sense.h fastlog.h stage3.c stage3.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 splice.c splice.h indel.c indel.h bitpack64-access.c bitpack64-access.h \
 bytecoding.c bytecoding.h univdiagdef.h univdiag.c univdiag.h sarray-read.c sarray-read.h \
 stage1hr.c stage1hr.h resulthr.c resulthr.h \
 datadir.c datadir.h mode.h \
 getopt.c getopt1.c getopt.h uniqscan.c

uniqscan_CC = $(PTHREAD_CC)
uniqscan_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1
uniqscan_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
uniqscan_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_uniqscan_SOURCES = $(UNIQSCAN_FILES)


UNIQSCANL_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h uint8list.c uint8list.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h \
 mapq.c mapq.h shortread.c shortread.h substring.c substring.h junction.c junction.h stage3hr.c stage3hr.h \
 spanningelt.c spanningelt.h cmet.c cmet.h atoi.c atoi.h \
 comp.h maxent.c maxent.h pairdef.h pair.c pair.h pairpool.c pairpool.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h cellpool.c cellpool.h stage2.c stage2.h \
 intron.c intron.h boyer-moore.c boyer-moore.h changepoint.c changepoint.h pbinom.c pbinom.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 chimera.c chimera.h sense.h fastlog.h stage3.c stage3.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 splice.c splice.h indel.c indel.h bitpack64-access.c bitpack64-access.h \
 stage1hr.c stage1hr.h resulthr.c resulthr.h \
 datadir.c datadir.h mode.h \
 getopt.c getopt1.c getopt.h uniqscan.c

uniqscanl_CC = $(PTHREAD_CC)
uniqscanl_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1 -DLARGE_GENOMES=1
uniqscanl_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
uniqscanl_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

dist_uniqscanl_SOURCES = $(UNIQSCANL_FILES)


GMAPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iit-write-univ.c iit-write-univ.h \
 iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h genome.c genome.h \
 genomicpos.c genomicpos.h compress-write.c compress-write.h genome-write.c genome-write.h \
 compress.c compress.h popcount.c popcount.h genome128_hr.c genome128_hr.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-write.c bitpack64-write.h \
 indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 table.c table.h tableuint.c tableuint.h tableuint8.c tableuint8.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 saca-k.c saca-k.h uint8list.c uint8list.h bytecoding.c bytecoding.h sarray-write.c sarray-write.h \
 gmapindex.c

gmapindex_CC = $(PTHREAD_CC)
gmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DUTILITYP=1
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapindex_SOURCES = $(GMAPINDEX_FILES)


SNPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h bitpack64-write.c bitpack64-write.h \
 indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 chrom.c chrom.h \
 complement.h md5.c md5.h bzip2.c bzip2.h sequence.c sequence.h genome.c genome.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h snpindex.c

snpindex_CC = $(PTHREAD_CC)
snpindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
snpindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
snpindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_snpindex_SOURCES = $(SNPINDEX_FILES)


CMETINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h access.c access.h \
 univinterval.c univinterval.h interval.h interval.c \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.h iit-read.c \
 complement.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h bitpack64-write.c bitpack64-write.h \
 indexdbdef.h indexdb.h indexdb.c indexdb-write.c indexdb-write.h \
 cmet.c cmet.h \
 complement.h md5.c md5.h bzip2.c bzip2.h sequence.c sequence.h genome.c genome.h \
 uintlist.c uintlist.h intlist.c intlist.h \
 saca-k.c saca-k.h uint8list.c uint8list.h bytecoding.c bytecoding.h sarray-write.c sarray-write.h \
 list.c list.h listdef.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h cmetindex.c

cmetindex_CC = $(PTHREAD_CC)
cmetindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
cmetindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
cmetindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_cmetindex_SOURCES = $(CMETINDEX_FILES)


ATOIINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h access.c access.h \
 univinterval.c univinterval.h interval.h interval.c \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.h iit-read.c \
 complement.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h genome128_hr.c genome128_hr.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h bitpack64-write.c bitpack64-write.h \
 indexdbdef.h indexdb.h indexdb.c indexdb-write.c indexdb-write.h \
 atoi.c atoi.h \
 complement.h md5.c md5.h bzip2.c bzip2.h sequence.c sequence.h genome.c genome.h \
 uintlist.c uintlist.h intlist.c intlist.h \
 saca-k.c saca-k.h uint8list.c uint8list.h bytecoding.c bytecoding.h sarray-write.c sarray-write.h \
 list.c list.h listdef.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h atoiindex.c

atoiindex_CC = $(PTHREAD_CC)
atoiindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
atoiindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
atoiindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_atoiindex_SOURCES = $(ATOIINDEX_FILES)


GET_GENOME_FILES = fopen.h bool.h types.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h sequence.c sequence.h \
 genome.c genome.h \
 genomicpos.c genomicpos.h chrom.c chrom.h \
 chrnum.c chrnum.h \
 datadir.c datadir.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h get-genome.c

get_genome_CC = $(PTHREAD_CC)
get_genome_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
get_genome_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_get_genome_SOURCES = $(GET_GENOME_FILES)


IIT_STORE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 doublelist.c doublelist.h \
 iit-write-univ.c iit-write-univ.h iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h table.c table.h chrom.c chrom.h \
 getopt.c getopt1.c getopt.h iit_store.c

iit_store_CC = $(PTHREAD_CC)
iit_store_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_store_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_store_LDADD = $(PTHREAD_LIBS)
dist_iit_store_SOURCES = $(IIT_STORE_FILES)


IIT_GET_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 complement.h parserange.c parserange.h getopt.c getopt1.c getopt.h iit_get.c

iit_get_CC = $(PTHREAD_CC)
iit_get_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_get_LDADD = $(PTHREAD_LIBS)
dist_iit_get_SOURCES = $(IIT_GET_FILES)

IIT_DUMP_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 intlistdef.h intlist.c intlist.h list.c list.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 getopt.c getopt1.c getopt.h iit_dump.c

iit_dump_CC = $(PTHREAD_CC)
iit_dump_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_dump_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_dump_LDADD = $(PTHREAD_LIBS)
dist_iit_dump_SOURCES = $(IIT_DUMP_FILES)


SAM_SORT_FILES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.h chrnum.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 stopwatch.c stopwatch.h access.c access.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 samflags.h samheader.c samheader.h samread.c samread.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h sam_sort.c

sam_sort_CC = $(PTHREAD_CC)
sam_sort_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
sam_sort_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
sam_sort_LDADD = $(PTHREAD_LIBS)
dist_sam_sort_SOURCES = $(SAM_SORT_FILES)
