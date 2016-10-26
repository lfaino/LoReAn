/* $Id: datadir.h 157221 2015-01-22 18:38:57Z twu $ */
#ifndef DATADIR_INCLUDED
#define DATADIR_INCLUDED

#include <stdio.h>

extern char *
Datadir_find_genomedir (char *user_genomedir);
extern char *
Datadir_find_genomesubdir (char **fileroot, char **dbversion,
			   char *user_genomedir, char *dbroot);
extern char *
Datadir_find_mapdir (char *user_mapdir, char *genomesubdir, char *fileroot);
extern void
Datadir_list_directory_multicol (FILE *fp, char *directory);

extern void
Datadir_list_directory (FILE *fp, char *directory);
extern void
Datadir_avail_gmap_databases (FILE *fp, char *user_genomedir);
extern void
Datadir_avail_maps (FILE *fp, char *user_mapdir, char *genomesubdir, char *fileroot);

#endif

