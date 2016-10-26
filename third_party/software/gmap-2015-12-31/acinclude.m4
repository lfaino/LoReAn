M4 Source Code

m4_define([DATE],['`date +%Y-%m-%d`'])

m4_include([config/pagesize.m4])
m4_include([config/madvise-flags.m4])
m4_include([config/mmap-flags.m4])
m4_include([config/acx_mmap_fixed.m4])
m4_include([config/acx_mmap_variable.m4])
m4_include([config/shm-flags.m4])

m4_include([config/ax_mpi.m4])
m4_include([config/acx_pthread.m4])

m4_include([config/builtin-popcount.m4])
m4_include([config/struct-stat64.m4])
m4_include([config/expand.m4])
m4_include([config/perl.m4])

m4_include([config/fopen.m4])
m4_include([config/asm-bsr.m4])
m4_include([config/sse2_shift_defect.m4])

m4_include([config/ax_gcc_x86_cpuid.m4])
m4_include([config/ax_gcc_x86_avx_xgetbv.m4])
m4_include([config/ax_check_compile_flag.m4])
m4_include([config/ax_ext.m4])
