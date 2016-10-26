static char rcsid[] = "$Id: samheader.c 157094 2015-01-21 00:33:35Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samheader.h"

#include <stdlib.h>
#include <string.h>
#include "mem.h"


#define CHUNK 1024

#ifdef USE_MPI
MPI_File
#else
FILE *
#endif
SAM_header_open_file (SAM_split_output_type split_output, char *split_output_root, bool appendp) {
#ifdef USE_MPI
  MPI_File output;
#else
  FILE *output;
  char *write_mode;
#endif
  char *filename, *suffix;

  if (split_output == OUTPUT_NONE) {

#ifdef USE_MPI
    /* output file name is passed in through split_output_root */
    if (appendp == true) {
      MPI_File_open(MPI_COMM_WORLD,split_output_root,MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
                    MPI_INFO_NULL,&output);
    } else {
      /* Need to remove existing file, if any */
      MPI_File_open(MPI_COMM_WORLD,split_output_root,MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                    MPI_INFO_NULL,&output);
      MPI_File_close(&output);
      MPI_File_open(MPI_COMM_WORLD,split_output_root,MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL,&output);
    }
    return output;

#else
    if (appendp == true) {
      write_mode = "a";
    } else {
      write_mode = "w";
    }

    filename = (char *) CALLOC(strlen(split_output_root)+1,sizeof(char));
    sprintf(filename,"%s",split_output_root);

    if ((output = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
      return (FILE *) NULL;
    } else {
      FREE(filename);
      return output;
    }
#endif

  } else {
    switch (split_output) {
    case OUTPUT_NONE: /* Handled above */ abort();

    case OUTPUT_NM: suffix = "nomapping"; break;

#ifdef GSNAP
    case OUTPUT_UU: suffix = "unpaired_uniq"; break;
    case OUTPUT_UT: suffix = "unpaired_transloc"; break;
    case OUTPUT_UM: suffix = "unpaired_mult"; break;
#else
    case OUTPUT_UU: suffix = "uniq"; break;
    case OUTPUT_UT: suffix = "transloc"; break;
    case OUTPUT_UM: suffix = "mult"; break;
#endif

#ifdef GSNAP
    case OUTPUT_UC: suffix = "unpaired_circular"; break;
    case OUTPUT_UX: suffix = "unpaired_mult_xs"; break;
#else
    case OUTPUT_UC: suffix = "circular"; break;
    case OUTPUT_UX: suffix = "mult_xs"; break;
#endif

    case OUTPUT_HU: suffix = "halfmapping_uniq"; break;
    case OUTPUT_HT: suffix = "halfmapping_transloc"; break;
    case OUTPUT_HM: suffix = "halfmapping_mult"; break;

    case OUTPUT_PI: suffix = "paired_uniq_inv"; break;
    case OUTPUT_PS: suffix = "paired_uniq_scr"; break;
    case OUTPUT_PL: suffix = "paired_uniq_long"; break;
    case OUTPUT_PM: suffix = "paired_mult"; break;

    case OUTPUT_CU: suffix = "concordant_uniq"; break;
    case OUTPUT_CT: suffix = "concordant_transloc"; break;
    case OUTPUT_CM: suffix = "concordant_mult"; break;

    case OUTPUT_HC: suffix = "halfmapping_circular"; break;
    case OUTPUT_PC: suffix = "paired_uniq_circular"; break;
    case OUTPUT_CC: suffix = "concordant_circular"; break;

    case OUTPUT_HX: suffix = "halfmapping_mult_xs"; break;
    case OUTPUT_PX: suffix = "paired_mult_xs"; break;
    case OUTPUT_CX: suffix = "concordant_mult_xs"; break;

    default:
      fprintf(stderr,"Cannot handle split output type %d\n",split_output);
      abort();
    }

    filename = (char *) CALLOC(strlen(split_output_root)+strlen(".")+strlen(suffix)+1,sizeof(char));
    sprintf(filename,"%s.%s",split_output_root,suffix);

#ifdef USE_MPI
    if (appendp == true) {
      MPI_File_open(MPI_COMM_WORLD,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
		    MPI_INFO_NULL,&output);
    } else {
      /* Need to remove existing file, if any */
      MPI_File_open(MPI_COMM_WORLD,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
		    MPI_INFO_NULL,&output);
      MPI_File_close(&output);
      MPI_File_open(MPI_COMM_WORLD,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		    MPI_INFO_NULL,&output);
    }
    FREE(filename);
    return output;

#else
    if (appendp == true) {
      write_mode = "a";
    } else {
      write_mode = "w";
    }

    if ((output = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
      return (FILE *) NULL;
    } else {
      FREE(filename);
      return output;
    }
#endif
  }
}


/* Called only by sam_sort */
Filestring_T
SAM_header_change_HD_tosorted (FILE *input, int headerlen) {
  Filestring_T fp;
  char buffer[CHUNK], c, c0, c1, c2;

  fp = Filestring_new(/*id*/0);

  /* @HD */
  while (headerlen > 0 && (c = fgetc(input)) != '\t') {
    PUTC(c,fp);
    headerlen--;
  }
  if (headerlen > 0) {
    PUTC('\t',fp);
    headerlen--;
  }

  /* VN */
  while (headerlen > 0 && (c = fgetc(input)) != '\t') {
    PUTC(c,fp);
    headerlen--;
  }
  if (headerlen > 0) {
    PUTC('\t',fp);
    headerlen--;
  }

  if (headerlen > 3) {
    /* SO: */
    c0 = fgetc(input);
    c1 = fgetc(input);
    c2 = fgetc(input);
    FPRINTF(fp,"%c%c%c",c0,c1,c2);
    headerlen -= 3;

    if (c0 == 'S' && c1 == 'O' && c2 == ':') {
      FPRINTF(fp,"coordinate\n");
      while (headerlen > 0 && fgetc(input) != '\n') {
	/* Skip given SO value */
	headerlen--;
      }
      headerlen--;
    }
  }

  while (headerlen > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,input);
    Filestring_puts(fp,buffer,/*strlength*/CHUNK);
    headerlen -= CHUNK;
  }
  if (headerlen > 0) {
    fread(buffer,sizeof(char),headerlen,input);
    Filestring_puts(fp,buffer,/*strlength*/headerlen);
  }

  return fp;
}


#ifdef USE_MPI
void
SAM_header_print_HD (MPI_File fp, int nworkers, bool orderedp) {

  MPI_File_write_shared(fp,"@HD",strlen("@HD"),MPI_CHAR,MPI_STATUS_IGNORE);
  MPI_File_write_shared(fp,"\tVN:1.0",strlen("\tVN:1.0"),MPI_CHAR,MPI_STATUS_IGNORE);
  if (nworkers > 1 && orderedp == false) {
    MPI_File_write_shared(fp,"\tSO:unsorted",strlen("\tSO:unsorted"),MPI_CHAR,MPI_STATUS_IGNORE);
  } else {
    /* Picard does not recognize type unknown */
    /* fprintf(fp,"\tSO:unknown"); */
    MPI_File_write_shared(fp,"\tSO:unsorted",strlen("\tSO:unsorted"),MPI_CHAR,MPI_STATUS_IGNORE);
  }
  MPI_File_write_shared(fp,"\n",1,MPI_CHAR,MPI_STATUS_IGNORE);

  return;
}

#else
void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp) {

  fprintf(fp,"@HD");
  fprintf(fp,"\tVN:1.0");	/* or 1.4 */
  if (nworkers > 1 && orderedp == false) {
    fprintf(fp,"\tSO:unsorted");
  } else {
    /* Picard does not recognize type unknown */
    /* fprintf(fp,"\tSO:unknown"); */
    fprintf(fp,"\tSO:unsorted");
  }
  fprintf(fp,"\n");

  return;
}
#endif


#ifdef USE_MPI
void
SAM_header_print_PG (MPI_File fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  MPI_File_write_shared(fp,"@PG",strlen("@PG"),MPI_CHAR,MPI_STATUS_IGNORE);
#ifdef GSNAP
  MPI_File_write_shared(fp,"\tID:GSNAP",strlen("\tID:GSNAP"),MPI_CHAR,MPI_STATUS_IGNORE);
  MPI_File_write_shared(fp,"\tPN:gsnap",strlen("\tPN:gsnap"),MPI_CHAR,MPI_STATUS_IGNORE);
#elif defined(PMAP)
  MPI_File_write_shared(fp,"\tID:PMAP",strlen("\tID:PMAP"),MPI_CHAR,MPI_STATUS_IGNORE);
  MPI_File_write_shared(fp,"\tPN:pmap",strlen("\tPN:pmap"),MPI_CHAR,MPI_STATUS_IGNORE);
#else
  MPI_File_write_shared(fp,"\tID:GMAP",strlen("\tID:GMAP"),MPI_CHAR,MPI_STATUS_IGNORE);
  MPI_File_write_shared(fp,"\tPN:gmap",strlen("\tPN:gmap"),MPI_CHAR,MPI_STATUS_IGNORE);
#endif
  MPI_File_write_shared(fp,PACKAGE_VERSION,strlen(PACKAGE_VERSION),MPI_CHAR,MPI_STATUS_IGNORE);

  MPI_File_write_shared(fp,"\tCL:",strlen("\tCL:"),MPI_CHAR,MPI_STATUS_IGNORE);
  argstart = &(argv[-optind]);
  MPI_File_write_shared(fp,argstart[0],strlen(argstart[0]),MPI_CHAR,MPI_STATUS_IGNORE);
  for (c = 1; c < argc + optind; c++) {
    MPI_File_write_shared(fp," ",1,MPI_CHAR,MPI_STATUS_IGNORE);
    MPI_File_write_shared(fp,argstart[c],strlen(argstart[c]),MPI_CHAR,MPI_STATUS_IGNORE);
  }
  MPI_File_write_shared(fp,"\n",1,MPI_CHAR,MPI_STATUS_IGNORE);

  return;
}

#else
void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  fprintf(fp,"@PG");
#ifdef GSNAP
  fprintf(fp,"\tID:GSNAP");
  fprintf(fp,"\tPN:gsnap");
#elif defined(PMAP)
  fprintf(fp,"\tID:PMAP");
  fprintf(fp,"\tPN:pmap");
#else
  fprintf(fp,"\tID:GMAP");
  fprintf(fp,"\tPN:gmap");
#endif
  fprintf(fp,"\tVN:%s",PACKAGE_VERSION);

  fprintf(fp,"\tCL:");
  argstart = &(argv[-optind]);
  fprintf(fp,"%s",argstart[0]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(fp," %s",argstart[c]);
  }
  fprintf(fp,"\n");

  return;
}
#endif


/* Called only by sam_sort */
int
SAM_header_length (int *lastchar, FILE *fp) {
  int headerlen = 0;
  int c;

  while (!feof(fp) && (c = getc(fp)) == '@') {
    headerlen++;
    while (!feof(fp) && (c = getc(fp)) != '\n') {
      headerlen++;
    }
    headerlen++;
  }
  /* headerlen++; -- Don't count in header, but as part of first SAM line */

  *lastchar = c;
  return headerlen;
}

