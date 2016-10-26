static char rcsid[] = "$Id: junction.c 166641 2015-05-29 21:13:04Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "junction.h"
#include "mem.h"
#include "complement.h"


#define T Junction_T
struct T {
  Junctiontype_T type;
  int nindels;
  Univcoord_T deletionpos;

  Chrpos_T splice_distance;
  int sensedir;
  double donor_prob;
  double acceptor_prob;
};


void
Junction_print (T this) {
  if (this == NULL) {
    printf("No junction\n");
  } else if (this->type == INS_JUNCTION) {
    printf("Insertion of %d\n",this->nindels);
  } else if (this->type == DEL_JUNCTION) {
    printf("Deletion of %d at %u\n",this->nindels,this->deletionpos);
  } else if (this->type == SPLICE_JUNCTION) {
    if (this->splice_distance == 0) {
      printf("Splice ambiguous with sense %d, prob %f and %f\n",
	     this->sensedir,this->donor_prob,this->acceptor_prob);
    } else {
      printf("Splice with sense %d of %u, prob %f and %f\n",
	     this->sensedir,this->splice_distance,this->donor_prob,this->acceptor_prob);
    }
  }
  return;
}

void
Junction_free (T *old) {
  FREE(*old);
  return;
}

void
Junction_gc (List_T *list) {
  List_T p;
  T old;

  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Junction_free(&old);
  }
  List_free(&(*list));
  return;
}

T
Junction_new_insertion (int nindels) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = INS_JUNCTION;
  new->nindels = nindels;
  new->deletionpos = 0;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_deletion (int nindels, Univcoord_T deletionpos) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = DEL_JUNCTION;
  new->nindels = nindels;
  new->deletionpos = deletionpos;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = SPLICE_JUNCTION;
  new->nindels = 0;
  new->deletionpos = 0;

  new->splice_distance = splice_distance;
  new->sensedir = sensedir;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}


T
Junction_new_chimera (int sensedir, double donor_prob, double acceptor_prob) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = CHIMERA_JUNCTION;
  new->nindels = 0;
  new->deletionpos = 0;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}

T
Junction_copy (T old) {
  T new = (T) MALLOC(sizeof(*new));

  new->type = old->type;
  new->nindels = old->nindels;
  new->deletionpos = old->deletionpos;

  new->splice_distance = old->splice_distance;
  new->sensedir = old->sensedir;
  new->donor_prob = old->donor_prob;
  new->acceptor_prob = old->acceptor_prob;

  return new;
}




Junctiontype_T
Junction_type (T this) {
  return this->type;
}

double
Junction_prob (T this) {
  return this->donor_prob + this->acceptor_prob;
}

int
Junction_sensedir (T this) {
  return this->sensedir;
}

double
Junction_donor_prob (T this) {
  return this->donor_prob;
}

double
Junction_acceptor_prob (T this) {
  return this->acceptor_prob;
}


int
Junction_nindels (T this) {
  return this->nindels;
}

int
Junction_adj (T this) {
  if (this->type == DEL_JUNCTION) {
    return +this->nindels;
  } else if (this->type == INS_JUNCTION) {
    return -this->nindels;
  } else {
    return 0;
  }
}




static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return sequence;
}

char *
Junction_deletion_string (T this, Genome_T genome, bool plusp) {
  char *deletion_string;
  
  deletion_string = (char *) MALLOC((this->nindels+1)*sizeof(char));
  /* printf("deletionpos = %u\n",this->deletionpos); */
  Genome_fill_buffer_simple(genome,this->deletionpos,this->nindels,deletion_string);
  if (plusp == false) {
    make_complement_inplace(deletion_string,this->nindels);
  }

  return deletion_string;
}


Chrpos_T
Junction_splice_distance (T this) {
  return this->splice_distance;
}

void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob) {
  this->splice_distance = distance;
  this->donor_prob = donor_prob;
  this->acceptor_prob = acceptor_prob;

  return;
}


