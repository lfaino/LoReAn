static char rcsid[] = "$Id: resulthr.c 155282 2014-12-12 19:42:54Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "resulthr.h"
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "stage3hr.h"


/* Assignment of resulttype */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



#define T Result_T
struct T {
  Resulttype_T resulttype;
  int id;
  void **array;
  int npaths;
  int first_absmq;
  int second_absmq;
  void **array2;
  int npaths2;
  int first_absmq2;
  int second_absmq2;
  double worker_runtime;

  SAM_split_output_type split_output;
  int strlength;
  char *string;
};


Resulttype_T
Result_resulttype (T this) {
  return this->resulttype;
}


char *
Pairtype_string (Pairtype_T pairtype) {
  switch (pairtype) {

  case CONCORDANT: return "concordant";
  case PAIRED_UNSPECIFIED: return "paired_unspecified";
  case PAIRED_INVERSION: return "paired_scramble";
  case PAIRED_SCRAMBLE: return "paired_scramble";
  case PAIRED_TOOLONG: return "paired_toolong";
  case CONCORDANT_TRANSLOCATIONS: return "concordant_translocations";
  case CONCORDANT_TERMINAL: return "concordant_terminal";
  case UNPAIRED: return "unpaired";
  case UNSPECIFIED: return "unspecified";
  default: 
    fprintf(stderr,"Unknown pairtype %d\n",pairtype);
    abort();
  }
  return "";
}


char *
Resulttype_string (Resulttype_T resulttype) {
  switch (resulttype) {
  case SINGLEEND_NOMAPPING: return "singleend_nomapping";
  case PAIREDEND_NOMAPPING: return "pairedend_nomapping";
  case SINGLEEND_UNIQ: return "singleend_uniq";
  case SINGLEEND_TRANSLOC: return "singleend_transloc";
  case SINGLEEND_MULT: return "singleend_mult";
  case PAIRED_UNIQ: return "paired_uniq";
  case PAIRED_MULT: return "paired_mult";
  case CONCORDANT_UNIQ: return "concordant_uniq";
  case CONCORDANT_TRANSLOC: return "concordant_transloc";
  case CONCORDANT_MULT: return "concordant_mult";
  case HALFMAPPING_UNIQ: return "halfmapping_uniq";
  case HALFMAPPING_TRANSLOC: return "halfmapping_transloc";
  case HALFMAPPING_MULT: return "halfmapping_mult";
  case UNPAIRED_UNIQ: return "unpaired_uniq";
  case UNPAIRED_TRANSLOC: return "unpaired_transloc";
  case UNPAIRED_MULT: return "unpaired_mult";
  default: 
    fprintf(stderr,"Unknown resulttype %d\n",resulttype);
    abort();
  }
  return "";
}


int
Result_id (T this) {
  return this->id;
}


void **
Result_array (int *npaths, int *first_absmq, int *second_absmq, T this) {
  *npaths = this->npaths;
  *first_absmq = this->first_absmq;
  *second_absmq = this->second_absmq;
  return this->array;
}


/* For second end, when not paired with first end */
void **
Result_array2 (int *npaths, int *first_absmq, int *second_absmq, T this) {
  *npaths = this->npaths2;
  *first_absmq = this->first_absmq2;
  *second_absmq = this->second_absmq2;
  return this->array2;
}

double
Result_worker_runtime (T this) {
  return this->worker_runtime;
}


T
Result_single_read_new (int id, void **resultarray, int npaths, int first_absmq, int second_absmq) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Stage3end_T stage3end;

  if (npaths == 0) {
    new->resulttype = SINGLEEND_NOMAPPING;
  } else if (Stage3end_chrnum((Stage3end_T) resultarray[0]) == 0) {
#if 0
    if (npaths == 1) {
      new->resulttype = SINGLEEND_TRANSLOC;
    } else {
      fprintf(stderr,"Unexpected: multiple singleend transloc\n");
      abort();
    }
#else
    new->resulttype = SINGLEEND_TRANSLOC;
#endif

  } else {
    if (npaths > 1) {
      new->resulttype = SINGLEEND_MULT;
    } else {
      stage3end = ((Stage3end_T *) resultarray)[0];
      if (Stage3end_anomalous_splice_p(stage3end) == true) {
	new->resulttype = SINGLEEND_TRANSLOC;
      } else {
	new->resulttype = SINGLEEND_UNIQ;
      }
    }
  }

  new->id = id;
  new->array = resultarray;
  new->npaths = npaths;
  new->first_absmq = first_absmq;
  new->second_absmq = second_absmq;

  return new;
}

T
Result_paired_read_new (int id, void **resultarray, int npaths, int first_absmq, int second_absmq,
			Pairtype_T final_pairtype) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Stage3pair_T stage3pair;

  if (npaths == 0) {
    abort();

  } else if (final_pairtype == CONCORDANT_TRANSLOCATIONS) {
#if 0
    if (npaths == 1) {
      new->resulttype = CONCORDANT_TRANSLOC;
    } else {
      fprintf(stderr,"Unexpected: multiple singleend transloc\n");
      abort();
    }
#else
    new->resulttype = CONCORDANT_TRANSLOC;
#endif

  } else if (final_pairtype == PAIRED_UNSPECIFIED) {
    if (npaths == 1) {
      new->resulttype = PAIRED_UNIQ;
    } else {
      new->resulttype = PAIRED_MULT;
    }

  } else if (final_pairtype == CONCORDANT) {
    if (npaths > 1) {
      new->resulttype = CONCORDANT_MULT;
    } else {
      stage3pair = ((Stage3pair_T *) resultarray)[0];
      if (Stage3pair_anomalous_splice_p(stage3pair) == true) {
	new->resulttype = CONCORDANT_TRANSLOC;
      } else {
	new->resulttype = CONCORDANT_UNIQ;
      }
    }

  } else {
    fprintf(stderr,"final_pairtype %d not recognized\n",final_pairtype);
    abort();
  }

  new->id = id;
  new->array = resultarray;
  new->npaths = npaths;
  new->first_absmq = first_absmq;
  new->second_absmq = second_absmq;

  return new;
}

T
Result_paired_as_singles_new (int id, void **hits5, int npaths5, int first_absmq5, int second_absmq5,
			      void **hits3, int npaths3, int first_absmq3, int second_absmq3) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Stage3end_T stage3end_5, stage3end_3;

  debug(printf("npaths5 = %d, npaths3 = %d\n",npaths5,npaths3));
  if (npaths5 == 0 && npaths3 == 0) {
    new->resulttype = PAIREDEND_NOMAPPING;
  } else if (npaths5 == 0 && npaths3 == 1) {
    stage3end_3 = (Stage3end_T) hits3[0];
    if (Stage3end_chrnum(stage3end_3) == 0) {
      new->resulttype = HALFMAPPING_TRANSLOC;
    } else if (Stage3end_anomalous_splice_p(stage3end_3) == true) {
      new->resulttype = HALFMAPPING_TRANSLOC;
    } else {
      new->resulttype = HALFMAPPING_UNIQ;
    }
  } else if (npaths5 == 1 && npaths3 == 0) {
    stage3end_5 = (Stage3end_T) hits5[0];
    if (Stage3end_chrnum(stage3end_5) == 0) {
      new->resulttype = HALFMAPPING_TRANSLOC;
    } else if (Stage3end_anomalous_splice_p(stage3end_5) == true) {
      new->resulttype = HALFMAPPING_TRANSLOC;
    } else {
      new->resulttype = HALFMAPPING_UNIQ;
    }
  } else if (npaths5 == 0 || npaths3 == 0) {
    new->resulttype = HALFMAPPING_MULT;
  } else if (npaths5 == 1 && npaths3 == 1) {
    stage3end_5 = (Stage3end_T) hits5[0];
    stage3end_3 = (Stage3end_T) hits3[0];
    if (Stage3end_chrnum(stage3end_5) == 0 ||
	Stage3end_chrnum(stage3end_3) == 0) {
      new->resulttype = UNPAIRED_TRANSLOC;
    } else if (Stage3end_anomalous_splice_p(stage3end_5) == true ||
	       Stage3end_anomalous_splice_p(stage3end_3) == true) {
      new->resulttype = UNPAIRED_TRANSLOC;
    } else {
      new->resulttype = UNPAIRED_UNIQ;
    }
  } else {
    new->resulttype = UNPAIRED_MULT;
  }

  new->id = id;
  new->array = hits5;
  new->npaths = npaths5;
  new->first_absmq = first_absmq5;
  new->second_absmq = second_absmq5;
  new->array2 = hits3;
  new->npaths2 = npaths3;
  new->first_absmq2 = first_absmq3;
  new->second_absmq2 = second_absmq3;

  return new;
}

void
Result_free (T *old) {
  int i;
  Stage3end_T stage3;
  Stage3pair_T stage3pair;

  switch ((*old)->resulttype) {
  case SINGLEEND_NOMAPPING: case PAIREDEND_NOMAPPING: break;

  case SINGLEEND_UNIQ: case SINGLEEND_TRANSLOC: case SINGLEEND_MULT:
    for (i = 0; i < (*old)->npaths; i++) {
      stage3 = (*old)->array[i];
      Stage3end_free(&stage3);
    }
    FREE_OUT((*old)->array);
    break;

  case PAIRED_UNIQ: case PAIRED_MULT: case CONCORDANT_UNIQ: case CONCORDANT_TRANSLOC: case CONCORDANT_MULT:
    for (i = 0; i < (*old)->npaths; i++) {
      stage3pair = (*old)->array[i];
      Stage3pair_free(&stage3pair);
    }
    FREE_OUT((*old)->array);
    break;

  default:
    for (i = 0; i < (*old)->npaths2; i++) {
      stage3 = (*old)->array2[i];
      Stage3end_free(&stage3);
    }
    FREE_OUT((*old)->array2);

    for (i = 0; i < (*old)->npaths; i++) {
      stage3 = (*old)->array[i];
      Stage3end_free(&stage3);
    }
    FREE_OUT((*old)->array);
  }

  FREE_OUT(*old);

  return;
}



