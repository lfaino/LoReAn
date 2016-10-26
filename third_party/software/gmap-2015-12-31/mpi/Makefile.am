

CC = $(MPICC)
srcdir = $(top_srcdir)/src

# This is a configure-time conditional
if MPI_FOUND
bin_PROGRAMS = mpi_gsnap
endif


MPI_GSNAP_FILES = $(srcdir)/mpidebug.c $(srcdir)/mpidebug.h \
 $(srcdir)/fopen.h $(srcdir)/bool.h $(srcdir)/types.h $(srcdir)/separator.h $(srcdir)/comp.h \
 $(srcdir)/except.c $(srcdir)/except.h $(srcdir)/assert.c $(srcdir)/assert.h $(srcdir)/mem.c $(srcdir)/mem.h \
 $(srcdir)/intlistdef.h $(srcdir)/intlist.c $(srcdir)/intlist.h $(srcdir)/listdef.h $(srcdir)/list.c $(srcdir)/list.h \
 $(srcdir)/littleendian.c $(srcdir)/littleendian.h $(srcdir)/bigendian.c $(srcdir)/bigendian.h \
 $(srcdir)/univinterval.c $(srcdir)/univinterval.h $(srcdir)/interval.c $(srcdir)/interval.h \
 $(srcdir)/uintlist.c $(srcdir)/uintlist.h \
 $(srcdir)/stopwatch.c $(srcdir)/stopwatch.h $(srcdir)/access.c $(srcdir)/access.h \
 $(srcdir)/filestring.c $(srcdir)/filestring.h \
 $(srcdir)/iit-read-univ.c $(srcdir)/iit-read-univ.h $(srcdir)/iitdef.h $(srcdir)/iit-read.c $(srcdir)/iit-read.h \
 $(srcdir)/md5.c $(srcdir)/md5.h $(srcdir)/complement.h $(srcdir)/bzip2.c $(srcdir)/bzip2.h $(srcdir)/sequence.c $(srcdir)/sequence.h $(srcdir)/reader.c $(srcdir)/reader.h \
 $(srcdir)/genomicpos.c $(srcdir)/genomicpos.h $(srcdir)/compress.c $(srcdir)/compress.h \
 $(srcdir)/genome.c $(srcdir)/genome.h \
 $(srcdir)/popcount.c $(srcdir)/popcount.h $(srcdir)/genome128_hr.c $(srcdir)/genome128_hr.h $(srcdir)/genome_sites.c $(srcdir)/genome_sites.h \
 $(srcdir)/bitpack64-read.c $(srcdir)/bitpack64-read.h $(srcdir)/bitpack64-readtwo.c $(srcdir)/bitpack64-readtwo.h \
 $(srcdir)/indexdbdef.h $(srcdir)/indexdb.c $(srcdir)/indexdb.h $(srcdir)/indexdb_hr.c $(srcdir)/indexdb_hr.h \
 $(srcdir)/oligo.c $(srcdir)/oligo.h \
 $(srcdir)/chrom.c $(srcdir)/chrom.h $(srcdir)/segmentpos.c $(srcdir)/segmentpos.h \
 $(srcdir)/chrnum.c $(srcdir)/chrnum.h \
 $(srcdir)/maxent_hr.c $(srcdir)/maxent_hr.h $(srcdir)/samflags.h $(srcdir)/samprint.c $(srcdir)/samprint.h \
 $(srcdir)/mapq.c $(srcdir)/mapq.h $(srcdir)/shortread.c $(srcdir)/shortread.h $(srcdir)/substring.c $(srcdir)/substring.h $(srcdir)/junction.c $(srcdir)/junction.h $(srcdir)/stage3hr.c $(srcdir)/stage3hr.h \
 $(srcdir)/spanningelt.c $(srcdir)/spanningelt.h $(srcdir)/cmet.c $(srcdir)/cmet.h $(srcdir)/atoi.c $(srcdir)/atoi.h \
 $(srcdir)/comp.h $(srcdir)/maxent.c $(srcdir)/maxent.h $(srcdir)/pairdef.h $(srcdir)/pair.c $(srcdir)/pair.h $(srcdir)/pairpool.c $(srcdir)/pairpool.h $(srcdir)/diag.c $(srcdir)/diag.h $(srcdir)/diagpool.c $(srcdir)/diagpool.h \
 $(srcdir)/orderstat.c $(srcdir)/orderstat.h $(srcdir)/oligoindex_hr.c $(srcdir)/oligoindex_hr.h $(srcdir)/cellpool.c $(srcdir)/cellpool.h $(srcdir)/stage2.c $(srcdir)/stage2.h \
 $(srcdir)/intron.c $(srcdir)/intron.h $(srcdir)/boyer-moore.c $(srcdir)/boyer-moore.h $(srcdir)/changepoint.c $(srcdir)/changepoint.h $(srcdir)/pbinom.c $(srcdir)/pbinom.h \
 $(srcdir)/dynprog.c $(srcdir)/dynprog.h $(srcdir)/dynprog_simd.c $(srcdir)/dynprog_simd.h \
 $(srcdir)/dynprog_single.c $(srcdir)/dynprog_single.h $(srcdir)/dynprog_genome.c $(srcdir)/dynprog_genome.h $(srcdir)/dynprog_cdna.c $(srcdir)/dynprog_cdna.h $(srcdir)/dynprog_end.c $(srcdir)/dynprog_end.h \
 $(srcdir)/gbuffer.c $(srcdir)/gbuffer.h $(srcdir)/translation.c $(srcdir)/translation.h \
 $(srcdir)/doublelist.c $(srcdir)/doublelist.h $(srcdir)/smooth.c $(srcdir)/smooth.h \
 $(srcdir)/chimera.c $(srcdir)/chimera.h $(srcdir)/sense.h $(srcdir)/fastlog.h $(srcdir)/stage3.c $(srcdir)/stage3.h \
 $(srcdir)/splicestringpool.c $(srcdir)/splicestringpool.h $(srcdir)/splicetrie_build.c $(srcdir)/splicetrie_build.h $(srcdir)/splicetrie.c $(srcdir)/splicetrie.h \
 $(srcdir)/splice.c $(srcdir)/splice.h $(srcdir)/indel.c $(srcdir)/indel.h $(srcdir)/bitpack64-access.c $(srcdir)/bitpack64-access.h \
 $(srcdir)/bytecoding.c $(srcdir)/bytecoding.h $(srcdir)/univdiagdef.h $(srcdir)/univdiag.c $(srcdir)/univdiag.h $(srcdir)/sarray-read.c $(srcdir)/sarray-read.h \
 $(srcdir)/stage1hr.c $(srcdir)/stage1hr.h \
 $(srcdir)/request.c $(srcdir)/request.h $(srcdir)/resulthr.c $(srcdir)/resulthr.h $(srcdir)/output.c $(srcdir)/output.h \
 $(srcdir)/master.c $(srcdir)/master.h \
 $(srcdir)/inbuffer.c $(srcdir)/inbuffer.h $(srcdir)/samheader.c $(srcdir)/samheader.h $(srcdir)/outbuffer.c $(srcdir)/outbuffer.h \
 $(srcdir)/datadir.c $(srcdir)/datadir.h $(srcdir)/mode.h \
 $(srcdir)/getopt.c $(srcdir)/getopt1.c $(srcdir)/getopt.h $(srcdir)/gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

mpi_gsnap_CC = $(MPICC)
mpi_gsnap_CFLAGS = $(MPI_CFLAGS) $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1 -DUSE_MPI=1
mpi_gsnap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
mpi_gsnap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

# Need nodist_ to prevent sources from being copied from src/ subdirectory
nodist_mpi_gsnap_SOURCES = $(MPI_GSNAP_FILES)


MPI_GMAP_FILES = $(srcdir)/mpidebug.c $(srcdir)/mpidebug.h \
 $(srcdir)/= $(srcdir)/fopen.h $(srcdir)/bool.h $(srcdir)/types.h $(srcdir)/separator.h $(srcdir)/comp.h \
 $(srcdir)/except.c $(srcdir)/except.h $(srcdir)/assert.c $(srcdir)/assert.h $(srcdir)/mem.c $(srcdir)/mem.h \
 $(srcdir)/intlistdef.h $(srcdir)/intlist.c $(srcdir)/intlist.h $(srcdir)/listdef.h $(srcdir)/list.c $(srcdir)/list.h \
 $(srcdir)/littleendian.c $(srcdir)/littleendian.h $(srcdir)/bigendian.c $(srcdir)/bigendian.h \
 $(srcdir)/univinterval.c $(srcdir)/univinterval.h $(srcdir)/interval.c $(srcdir)/interval.h $(srcdir)/uintlist.c $(srcdir)/uintlist.h \
 $(srcdir)/stopwatch.c $(srcdir)/stopwatch.h $(srcdir)/access.c $(srcdir)/access.h \
 $(srcdir)/filestring.c $(srcdir)/filestring.h \
 $(srcdir)/iit-read-univ.c $(srcdir)/iit-read-univ.h $(srcdir)/iitdef.h $(srcdir)/iit-read.c $(srcdir)/iit-read.h \
 $(srcdir)/md5.c $(srcdir)/md5.h $(srcdir)/complement.h $(srcdir)/bzip2.c $(srcdir)/bzip2.h $(srcdir)/sequence.c $(srcdir)/sequence.h $(srcdir)/reader.c $(srcdir)/reader.h \
 $(srcdir)/genomicpos.c $(srcdir)/genomicpos.h $(srcdir)/compress.c $(srcdir)/compress.h $(srcdir)/compress-write.c $(srcdir)/compress-write.h \
 $(srcdir)/gbuffer.c $(srcdir)/gbuffer.h $(srcdir)/genome.c $(srcdir)/genome.h \
 $(srcdir)/popcount.c $(srcdir)/popcount.h $(srcdir)/genome128_hr.c $(srcdir)/genome128_hr.h $(srcdir)/genome_sites.c $(srcdir)/genome_sites.h \
 $(srcdir)/genome-write.c $(srcdir)/genome-write.h \
 $(srcdir)/bitpack64-read.c $(srcdir)/bitpack64-read.h $(srcdir)/bitpack64-readtwo.c $(srcdir)/bitpack64-readtwo.h \
 $(srcdir)/indexdbdef.h $(srcdir)/indexdb.c $(srcdir)/indexdb.h $(srcdir)/indexdb_hr.c $(srcdir)/indexdb_hr.h \
 $(srcdir)/oligo.c $(srcdir)/oligo.h $(srcdir)/block.c $(srcdir)/block.h \
 $(srcdir)/chrom.c $(srcdir)/chrom.h $(srcdir)/segmentpos.c $(srcdir)/segmentpos.h \
 $(srcdir)/chrnum.c $(srcdir)/chrnum.h $(srcdir)/uinttable.c $(srcdir)/uinttable.h $(srcdir)/gregion.c $(srcdir)/gregion.h \
 $(srcdir)/matchdef.h $(srcdir)/match.c $(srcdir)/match.h $(srcdir)/matchpool.c $(srcdir)/matchpool.h \
 $(srcdir)/diagnostic.c $(srcdir)/diagnostic.h $(srcdir)/stage1.c $(srcdir)/stage1.h \
 $(srcdir)/diagdef.h $(srcdir)/diag.c $(srcdir)/diag.h $(srcdir)/diagpool.c $(srcdir)/diagpool.h \
 $(srcdir)/cmet.c $(srcdir)/cmet.h $(srcdir)/atoi.c $(srcdir)/atoi.h \
 $(srcdir)/orderstat.c $(srcdir)/orderstat.h $(srcdir)/oligoindex_hr.c $(srcdir)/oligoindex_hr.h \
 $(srcdir)/scores.h $(srcdir)/intron.c $(srcdir)/intron.h $(srcdir)/maxent.c $(srcdir)/maxent.h $(srcdir)/maxent_hr.c $(srcdir)/maxent_hr.h $(srcdir)/samflags.h $(srcdir)/pairdef.h $(srcdir)/pair.c $(srcdir)/pair.h \
 $(srcdir)/pairpool.c $(srcdir)/pairpool.h $(srcdir)/cellpool.c $(srcdir)/cellpool.h $(srcdir)/stage2.c $(srcdir)/stage2.h \
 $(srcdir)/doublelist.c $(srcdir)/doublelist.h $(srcdir)/smooth.c $(srcdir)/smooth.h \
 $(srcdir)/splicestringpool.c $(srcdir)/splicestringpool.h $(srcdir)/splicetrie_build.c $(srcdir)/splicetrie_build.h $(srcdir)/splicetrie.c $(srcdir)/splicetrie.h \
 $(srcdir)/boyer-moore.c $(srcdir)/boyer-moore.h \
 $(srcdir)/dynprog.c $(srcdir)/dynprog.h $(srcdir)/dynprog_simd.c $(srcdir)/dynprog_simd.h \
 $(srcdir)/dynprog_single.c $(srcdir)/dynprog_single.h $(srcdir)/dynprog_genome.c $(srcdir)/dynprog_genome.h $(srcdir)/dynprog_cdna.c $(srcdir)/dynprog_cdna.h $(srcdir)/dynprog_end.c $(srcdir)/dynprog_end.h \
 $(srcdir)/translation.c $(srcdir)/translation.h \
 $(srcdir)/pbinom.c $(srcdir)/pbinom.h $(srcdir)/changepoint.c $(srcdir)/changepoint.h $(srcdir)/sense.h $(srcdir)/fastlog.h $(srcdir)/stage3.c $(srcdir)/stage3.h \
 $(srcdir)/request.c $(srcdir)/request.h $(srcdir)/result.c $(srcdir)/result.h $(srcdir)/output.c $(srcdir)/output.h \
 $(srcdir)/inbuffer.c $(srcdir)/inbuffer.h $(srcdir)/samheader.c $(srcdir)/samheader.h $(srcdir)/outbuffer.c $(srcdir)/outbuffer.h \
 $(srcdir)/chimera.c $(srcdir)/chimera.h $(srcdir)/datadir.c $(srcdir)/datadir.h \
 $(srcdir)/getopt.c $(srcdir)/getopt1.c $(srcdir)/getopt.h $(srcdir)/gmap.c

# Note: dist_ commands get read by bootstrap, and don't follow the flags

mpi_gmap_CC = $(MPICC)
mpi_gmap_CFLAGS = $(MPI_CFLAGS) $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(SIMD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DUSE_MPI=1
mpi_gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
mpi_gmap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)

nodist_mpi_gmap_SOURCES = $(MPI_GMAP_FILES)


