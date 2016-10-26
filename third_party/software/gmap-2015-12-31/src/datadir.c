static char rcsid[] = "$Id: datadir.c 73988 2012-09-13 23:55:05Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "datadir.h"
#include <stdio.h>
#include <stdlib.h>		/* For getenv */
#include <string.h>
#include <strings.h>		/* For rindex */
#include <pwd.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed for dirent.h */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include <math.h>		/* for qsort */
#include <string.h>		/* for strcmp */
#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "list.h"

/* Note: GMAPDB is defined externally by configure */
#ifndef GMAPDB
#error A default value for GMAPDB was not provided to configure.  Please do so, or edit the Makefile
#endif

static char *
read_config_file (FILE *fp, char *tag) {
  char *directory, seentag[1024], dirbuffer[1024], Buffer[1024];

  while (fgets(Buffer,1024,fp) != NULL) {
    if (sscanf(Buffer,"%s=%s",seentag,dirbuffer) > 0 && !strcmp(seentag,tag)) {
      directory = (char *) CALLOC(strlen(dirbuffer)+1,sizeof(char));
      strcpy(directory,dirbuffer);
      return directory;
    }
  }
  return NULL;
}


static FILE *
find_homedir_config () {
  FILE *fp = NULL;
  struct passwd *p;
  char *user, *configfile;

  if ((user = getenv("USER")) != NULL) {
    if ((p = getpwnam(user)) != NULL) {
      configfile = (char *) CALLOC(strlen(p->pw_dir)+strlen("/")+strlen(".gmaprc")+1,sizeof(char));
      sprintf(configfile,"%s/.gmaprc",p->pw_dir);
      fp = FOPEN_READ_TEXT(configfile);
      FREE(configfile);
    }
  }
  return fp;
}
    


static char *
find_fileroot (char *genomesubdir, char *genomedir, char *dbroot) {
  char *fileroot, *filename, *p;
  struct dirent *entry;
  DIR *dp;

  if ((dp = opendir(genomesubdir)) == NULL) {
    /* Problem found.  Try to diagnose */
    if ((dp = opendir(genomedir)) == NULL) {
      fprintf(stderr,"Unable to find genome directory %s.  Either recompile the GMAP package\n",genomedir);
      fprintf(stderr,"to have the correct default directory (seen by doing gmap --version),\n");
      fprintf(stderr,"or use the -D flag to gmap to specify the correct genome directory.\n");
      exit(9);
    } else {
      fprintf(stderr,"Unable to find genome %s in directory %s.\n",dbroot,genomedir);
      fprintf(stderr,"Make sure you have typed the genome correctly, or use the -D flag\n");
      fprintf(stderr,"(or the -F flag for cmetindex or atoiindex) to specify a directory.\n");
      fprintf(stderr,"For example, '-D .' specifies this directory.\n");
      exit(9);
    }
  }
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if ((p = rindex(filename,'.')) != NULL) {
      if (!strcmp(p,".version")) {
	fileroot = (char *) CALLOC(p - &(filename[0]) + 1,sizeof(char));
	strncpy(fileroot,filename,p-&(filename[0]));
	if (closedir(dp) < 0) {
	  fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
	}
	return fileroot;
      }
    }
  }

  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }
  fprintf(stderr,"Unable to find file ending with .version in directory %s\n",genomesubdir);
  exit(9);
}




static char *
get_dbversion (char *filename) {
  char *dbversion = NULL, Buffer[100], *p;
  FILE *fp;

  fp = FOPEN_READ_TEXT(filename);
  if (!fp) {
    return NULL;
  } else if (fgets(Buffer,100,fp) == NULL) {
    fclose(fp);
    return NULL;
  } else {
    if ((p = rindex(Buffer,'\n')) != NULL) {
      *p = '\0';
    }
    fclose(fp);
  }
  
  dbversion = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
  strcpy(dbversion,Buffer);
  return dbversion;
}


char *
Datadir_find_genomedir (char *user_genomedir) {
  FILE *fp;
  char *genomedir;

  if (user_genomedir != NULL) {
    genomedir = (char *) CALLOC(strlen(user_genomedir)+1,sizeof(char));
    strcpy(genomedir,user_genomedir);

  } else if (getenv("GMAPDB") != NULL) {
    /* Use genomedir provided by environment variable */
    genomedir = (char *) CALLOC(strlen(getenv("GMAPDB"))+1,sizeof(char));
    strcpy(genomedir,getenv("GMAPDB"));

  } else if ((fp = FOPEN_READ_TEXT("./.gmaprc")) != NULL) {
    genomedir = read_config_file(fp,"GMAPDB");
    fclose(fp);

  } else if ((fp = find_homedir_config()) != NULL) {
    genomedir = read_config_file(fp,"GMAPDB");
    fclose(fp);

  } else {
    genomedir = (char *) CALLOC(strlen(GMAPDB)+1,sizeof(char));
    strcpy(genomedir,GMAPDB);
  }

  return genomedir;
}


/* Allocates space for genomesubdir, fileroot, and dbversion */
char *
Datadir_find_genomesubdir (char **fileroot, char **dbversion,
			   char *user_genomedir, char *dbroot) {
  FILE *fp;
  char *genomesubdir, *genomedir, *filename, *p, *dbrootdir, *newgenomedir;

  /* First get genomedir */
  genomedir = Datadir_find_genomedir(user_genomedir);

  /* Append directory part of dbroot */
  if ((p = rindex(dbroot,'/')) != NULL) {
    *p = '\0';
    p++;
    dbrootdir = dbroot;
    dbroot = p;
    
    newgenomedir = (char *) CALLOC(strlen(genomedir)+strlen("/")+strlen(dbrootdir)+1,sizeof(char));
    sprintf(newgenomedir,"%s/%s",genomedir,dbrootdir);
    FREE(genomedir);
    genomedir = newgenomedir;
  }

  /* Find version file */
  filename = (char *) CALLOC(strlen(genomedir) + strlen("/") + strlen(dbroot) + 
			     strlen(".version") + 1,sizeof(char));
  sprintf(filename,"%s/%s.version",genomedir,dbroot);
  if ((fp = FOPEN_READ_TEXT(filename)) != NULL) {
    /* Found in top-level genomedir */
    fclose(fp);
    FREE(filename);
    genomesubdir = genomedir;
    *fileroot = (char *) CALLOC(strlen(dbroot)+1,sizeof(char));
    strcpy(*fileroot,dbroot);

  } else {
    FREE(filename);
    genomesubdir = (char *) CALLOC(strlen(genomedir) + strlen("/") + strlen(dbroot) + 1,sizeof(char));
    sprintf(genomesubdir,"%s/%s",genomedir,dbroot);

    if ((*fileroot = find_fileroot(genomesubdir,genomedir,dbroot)) != NULL) {
      /* Found in subdirectory */
       FREE(genomedir);
    } else {
      fprintf(stderr,"Error: Can't open genome files in %s or %s.\n",genomedir,genomesubdir);
      fprintf(stderr,"       Please specify directory using -D flag, GMAPDB environment variable,\n");
      fprintf(stderr,"       or a configuration file .gmaprc with the line GMAPDB=<directory>;\n");
      fprintf(stderr,"       or recompile the GMAP package using the --with-gmapdb flag to configure.\n");
      Datadir_avail_gmap_databases(stderr,user_genomedir);
      exit(9);
    }
  }
    
  filename = (char *) CALLOC(strlen(genomesubdir) + strlen("/") + strlen(*fileroot) + strlen(".version") + 1,
			     sizeof(char));
  sprintf(filename,"%s/%s.version",genomesubdir,*fileroot);
  if ((*dbversion = get_dbversion(filename)) == NULL) {
    /* Something wrong with version file.  Use dbroot instead */
    *dbversion = (char *) CALLOC(strlen(dbroot)+1,sizeof(char));
    strcpy(*dbversion,dbroot);
  }

  FREE(filename);

  return genomesubdir;
}


char *
Datadir_find_mapdir (char *user_mapdir, char *genomesubdir, char *fileroot) {
  char *mapdir;

  if (user_mapdir != NULL) {
    mapdir = (char *) CALLOC(strlen(user_mapdir)+1,sizeof(char));
    strcpy(mapdir,user_mapdir);
  } else {
    mapdir = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".maps")+1,sizeof(char));
    sprintf(mapdir,"%s/%s.maps",genomesubdir,fileroot);
  }

  return mapdir;
}


void
Datadir_list_directory_multicol (FILE *fp, char *directory) {
  char *filename;
  struct dirent *entry;
  DIR *dp;
  int pos = 0;

  if ((dp = opendir(directory)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",directory);
    exit(9);
  }
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (filename[0] != '.') {
      if (pos == 0) {
	fprintf(fp,"     ");
	pos += strlen("     ");
      } else {
	fprintf(fp," ");
	pos++;
	while (pos % 10 != 0) {
	  printf(" ");
	  pos++;
	}
      }
      fprintf(fp,"%s",filename);
      pos += strlen(filename);
      if (pos > 60) {
	fprintf(fp,"\n");
	pos = 0;
      }
    }
  }
  fprintf(fp,"\n");
  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",directory);
  }

  return;
}

void
Datadir_list_directory (FILE *fp, char *directory) {
  char *filename;
  struct dirent *entry;
  DIR *dp;

  if ((dp = opendir(directory)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",directory);
    exit(9);
  }
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    fprintf(fp,"%s\n",filename);
  }
  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",directory);
  }

  return;
}

static int
strcmp_cmp (const void *a, const void *b) {
  return strcmp(* (char * const *) a, * (char * const *) b);
}


void
Datadir_avail_gmap_databases (FILE *fp, char *user_genomedir) {
  char *genomedir;
  struct dirent *entry;
  char *filename;
  DIR *dp;
  List_T databases = NULL;
  char **array;
  int n, i;

  genomedir = Datadir_find_genomedir(user_genomedir);
  fprintf(fp,"Available gmap databases in directory %s:\n",genomedir);

  if ((dp = opendir(genomedir)) == NULL) {
    fprintf(stderr,"Unable to open genomedir %s\n",genomedir);
    exit(9);
  }
  while ((entry = readdir(dp)) != NULL) {
    filename = (char *) CALLOC(strlen(genomedir)+strlen("/")+strlen(entry->d_name)+strlen("/")+
			       strlen(entry->d_name)+strlen(".version")+1,sizeof(char));
    sprintf(filename,"%s/%s/%s.version",genomedir,entry->d_name,entry->d_name);
    if (Access_file_exists_p(filename) == true) {
      FREE(filename);
      filename = (char *) CALLOC(strlen(entry->d_name)+1,sizeof(char));
      strcpy(filename,entry->d_name);
      databases = List_push(databases,(void *) filename);
    } else {
      FREE(filename);
    }
  }
  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close genomedir %s\n",genomedir);
  }

  if ((n = List_length(databases)) == 0) {
    fprintf(fp,"  (none found)\n");
  } else {
    array = (char **) List_to_array(databases,NULL);
    qsort(array,n,sizeof(char *),strcmp_cmp);
    for (i = 0; i < n; i++) {
      fprintf(fp,"%s\n",array[i]);
      FREE(array[i]);
    }
    FREE(array);
    List_free(&databases);
  }

  FREE(genomedir);
  return;
}

void
Datadir_avail_maps (FILE *fp, char *user_mapdir, char *genomesubdir, char *fileroot) {
  char *mapdir;
  struct dirent *entry;
  char *filename;
  DIR *dp;
  List_T maps = NULL;
  char **array;
  int n, i;

  mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
  fprintf(fp,"Available maps in directory %s:\n",mapdir);

  if ((dp = opendir(mapdir)) == NULL) {
    fprintf(stderr,"Unable to open mapdir %s\n",mapdir);
    exit(9);
  }
  while ((entry = readdir(dp)) != NULL) {
    if (entry->d_name[0] != '.') {
      filename = (char *) CALLOC(strlen(mapdir)+strlen("/")+strlen(entry->d_name)+1,
				 sizeof(char));
      sprintf(filename,"%s/%s",mapdir,entry->d_name);
      
      if (Access_file_exists_p(filename) == true) {
	FREE(filename);
	filename = (char *) CALLOC(strlen(entry->d_name)+1,sizeof(char));
	strcpy(filename,entry->d_name);
	maps = List_push(maps,(void *) filename);
      } else {
	FREE(filename);
      }
    }
  }
  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close mapdir %s\n",mapdir);
  }

  if ((n = List_length(maps)) == 0) {
    fprintf(fp,"  (none found)\n");
  } else {
    array = (char **) List_to_array(maps,NULL);
    qsort(array,n,sizeof(char *),strcmp_cmp);
    for (i = 0; i < n; i++) {
      fprintf(fp,"%s\n",array[i]);
      FREE(array[i]);
    }
    FREE(array);
    List_free(&maps);
  }

  FREE(mapdir);
  return;

}

