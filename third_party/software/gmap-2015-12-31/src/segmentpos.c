static char rcsid[] = "$Id: segmentpos.c 155282 2014-12-12 19:42:54Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "segmentpos.h"
#include <stdio.h>
#include <string.h>		/* For strcmp */
#include "mem.h"
#include "intlist.h"
#include "separator.h"
#include "univinterval.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define ONEBASEDP 1		/* 1-based coordinates.  Also defined in pair.c */

#define MAXACCESSIONS 10

#define T Segmentpos_T
struct T {
  Chrom_T chrom;		/* 4 bytes */
  Chrpos_T chrpos1;		/* 4 bytes */
  Chrpos_T chrpos2;		/* 4 bytes */
  Chrpos_T length;		/* 4 bytes */
  int type;			/* 4 bytes */
  bool revcompp;
};


Chrom_T
Segmentpos_chrom (T this) {
  return this->chrom;
}

Chrpos_T
Segmentpos_chrpos1 (T this) {
  return this->chrpos1;
}

Chrpos_T
Segmentpos_chrpos2 (T this) {
  return this->chrpos2;
}

Chrpos_T
Segmentpos_length (T this) {
  return this->length;
}

int
Segmentpos_type (T this) {
  return this->type;
}

bool
Segmentpos_revcompp (T this) {
  return this->revcompp;
}


T
Segmentpos_new (Chrom_T chrom, Chrpos_T chrpos1, Chrpos_T chrpos2, 
		bool revcompp, Chrpos_T length, int type) {
  T new = (T) MALLOC(sizeof(*new));

  new->chrom = chrom;
  new->chrpos1 = chrpos1;
  new->chrpos2 = chrpos2;
  new->revcompp = revcompp;
  new->length = length;
  new->type = type;
  return new;
}

void
Segmentpos_free (T *old) {
  Chrom_free(&(*old)->chrom);
  FREE(*old);
  return;
}

void
Segmentpos_print (Filestring_T fp, T this, char *acc, Univcoord_T chroffset) {
  FPRINTF(fp,"%s\t%u\t%s\t%u\t%u\n",acc,chroffset+this->chrpos1,Chrom_string(this->chrom),this->chrpos1,this->length);
  return;
}

int
Segmentpos_compare_alpha (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  if ((cmp = Chrom_cmp_alpha(a->chrom,b->chrom)) != 0) {
    return cmp;
  } else if (a->chrpos1 < b->chrpos1) {
    return -1;
  } else if (b->chrpos1 < a->chrpos1) {
    return 1;
  } else {
    return 0;
  }
}

int
Segmentpos_compare_numeric_alpha (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  if ((cmp = Chrom_cmp_numeric_alpha(a->chrom,b->chrom)) != 0) {
    return cmp;
  } else if (a->chrpos1 < b->chrpos1) {
    return -1;
  } else if (b->chrpos1 < a->chrpos1) {
    return 1;
  } else {
    return 0;
  }
}

int
Segmentpos_compare_chrom (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  if ((cmp = Chrom_cmp_chrom(a->chrom,b->chrom)) != 0) {
    return cmp;
  } else if (a->chrpos1 < b->chrpos1) {
    return -1;
  } else if (b->chrpos1 < a->chrpos1) {
    return 1;
  } else {
    return 0;
  }
}

int
Segmentpos_compare_order (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  if ((cmp = Chrom_cmp_order(a->chrom,b->chrom)) != 0) {
    return cmp;
  } else if (a->chrpos1 < b->chrpos1) {
    return -1;
  } else if (b->chrpos1 < a->chrpos1) {
    return 1;
  } else {
    return 0;
  }
}

static bool
altstrain_sufficient_p (int *indices, int nindices, Univ_IIT_T contig_iit, Univcoord_T position1,
			Univcoord_T position2, char *align_strain) {
  int i;
  Univcoord_T contig_start, contig_end;
  Chrpos_T contig_length;
  int index, contig_straintype;
  Univinterval_T interval;

  i = 0;
  while (i < nindices) {
    index = indices[i];
    interval = Univ_IIT_interval(contig_iit,index);
    contig_straintype = Univinterval_type(interval);
    if (!strcmp(Univ_IIT_typestring(contig_iit,contig_straintype),align_strain)) {
      contig_start = Univinterval_low(interval);
      contig_length = Univinterval_length(interval);
      contig_end = contig_start + contig_length;

      if (contig_start <= position1 && contig_end >= position2) {
	debug(printf("Altstrain is sufficient: %u..%u subsumes %u..%u\n",
		     contig_start,contig_end,position1,position2));
	return true;
      }
    }
    i++;
  }

  debug(printf("Altstrain is not sufficient\n"));
  return false;
}


static bool
contig_print_p (Univ_IIT_T contig_iit, int contig_straintype, bool referencealignp,
		char *align_strain, bool printreferencep, bool printaltp) {
  if (contig_straintype == 0) {
    /* Contig is from reference strain */
    if (printreferencep == true) {
      return true;
    } else {
      return false;
    }

  } else if (referencealignp == true || align_strain == NULL) {
    /* Contig is from an alternate strain */
    if (printaltp == true) {
      return true;
    } else {
      return false;
    }

  } else if (strcmp(Univ_IIT_typestring(contig_iit,contig_straintype),align_strain)) {
    /* Contig is from a non-relevant alternate strain */
    return false;

  } else {
    /* Contig is from the aligned alternate strain */
    if (printaltp == true) {
      return true;
    } else {
      return false;
    }
  }
}


void
Segmentpos_print_accessions (Filestring_T fp, Univ_IIT_T contig_iit, Univcoord_T position1,
			     Univcoord_T position2, bool referencealignp, 
                             char *align_strain) {
  Univcoord_T contig_start;
  Chrpos_T contig_length;
  int relstart, relend;		/* Need to be signed int, not long or unsigned long */
  int index, contig_straintype, i = 0;
  char *label, *comma1, *comma2, firstchar;
  int *indices, nindices, j;
  Univinterval_T interval;
  bool printreferencep, printaltp, firstprintp = false, allocp;

  FPRINTF(fp,"    Accessions: ");

  indices = Univ_IIT_get(&nindices,contig_iit,position1,position2);
  if (referencealignp == true) {
    printreferencep = true;
    printaltp = false;
  } else if (altstrain_sufficient_p(indices,nindices,contig_iit,position1,position2,align_strain) == true) {
    printreferencep = false;
    printaltp = true;
  } else {
    printreferencep = true;
    printaltp = true;
  }

  j = 0;
  while (j < nindices && i < MAXACCESSIONS) {
    index = indices[j];
    interval = Univ_IIT_interval(contig_iit,index);
    contig_straintype = Univinterval_type(interval);
    if (contig_print_p(contig_iit,contig_straintype,referencealignp,align_strain,
		       printreferencep,printaltp) == true) {
      contig_start = Univinterval_low(interval);
      contig_length = Univinterval_length(interval);

      relstart = position1 - contig_start;
      if (relstart < 0) {
	relstart = 0;
      }
      relend = position2 - contig_start;
      if ((Chrpos_T) relend > contig_length) {
	relend = contig_length;
      }

      comma1 = Genomicpos_commafmt((Univcoord_T) (relstart + ONEBASEDP));
      comma2 = Genomicpos_commafmt((Univcoord_T) (relend + ONEBASEDP));
      
      if (firstprintp == true) {
	FPRINTF(fp,"; ");
      } else {
	firstprintp = true;
      }

      firstchar = Univ_IIT_annotation_firstchar(contig_iit,index);
      if (firstchar == '-') {
	FPRINTF(fp,"[-]");
      }

      label = Univ_IIT_label(contig_iit,index,&allocp);
      FPRINTF(fp,"%s",label);
      if (allocp == true) {
	FREE(label);
      }

      if (referencealignp == false && contig_straintype == 0) {
	FPRINTF(fp,"[reference strain]");
      }
      FPRINTF(fp,":%s%s%s (out of %u bp)",comma1,SEPARATOR,comma2,contig_length);

      FREE(comma2);
      FREE(comma1);

      i++;
    }
    j++;
  }
  FPRINTF(fp,"\n");

  if (indices != NULL) {
    FREE(indices);
  }

  return;
}


