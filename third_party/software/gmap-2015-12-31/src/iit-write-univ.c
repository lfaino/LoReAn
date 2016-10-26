static char rcsid[] = "$Id: iit-write-univ.c 153948 2014-11-24 17:46:46Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "iit-write-univ.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strlen and strcmp */
#include "assert.h"
#include "mem.h"
#include "fopen.h"
#include "types.h"		/* For HAVE_64_BIT and UINT8 */
#include "univinterval.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


typedef struct Node_T *Node_T;
struct Node_T {
  int index;
  Univcoord_T value;
  int a;
  int b;
  Node_T left;
  Node_T right;
};


static Node_T
Node_new () {
  Node_T new = (Node_T) MALLOC(sizeof(*new));

  new->index = -1;
  new->left = new->right = NULL;
  return new;
}

static void 
Node_gc (Node_T *old) {
  if (*old) {
    Node_gc(&((*old)->left));
    Node_gc(&((*old)->right));
    FREE(*old);
  }
  return;
}


#define T IIT_T


/************************************************************************/

static int 
is_right_split (int i, int j, int r, Univcoord_T value, int *sigmas, 
		struct Univinterval_T *intervals) {
  int iota, lambda;

  iota = 0;
  for (lambda = i; lambda <= j; lambda++) {
    if (Univinterval_is_contained(value,intervals,sigmas[lambda]) == true) {
      iota = lambda;
    }
  }
  return (iota == r);
}  

static bool
is_sorted (int array[], int i, int j, Univcoord_T (*endpoint)(), struct Univinterval_T *intervals) {
  int lambda;

  for (lambda = i; lambda <= j - 1; lambda++) {
    if (endpoint(intervals,array[lambda]) > endpoint(intervals,array[lambda + 1])) {
      return false;
    }
  }
  return true;
}

static bool
is_empty (int array[], int i, int j) {
  int lambda;

  for (lambda = i; lambda <= j; lambda++) {
    if (array[lambda] != 0) {
      return false;
    }
  }
  return true;
}

static bool
is_valid_input (int i, int j, int *sigmas, int *omegas, struct Univinterval_T *intervals) {
  int lambda, iota;

  assert(is_sorted (sigmas, i, j, Univinterval_array_low, intervals));
  assert(is_empty (omegas, i, j));
  for (lambda = i; lambda <= j; lambda++) {
    iota = sigmas[lambda];
    assert(Univinterval_is_contained(Univinterval_array_low(intervals,iota),intervals,iota)
	   || Univinterval_is_contained(Univinterval_array_high(intervals,iota),intervals,iota));
	   
  }
  return true;
}

/************************************************************************/

static bool
Node_is_valid_output (Node_T node, int i, int j, int *sigmas, int *omegas, 
		      struct Univinterval_T *intervals) {
  int lambda;

  assert((i <= node->a) && (node->b <= j));
  assert(is_sorted (sigmas, node->a, node->b, Univinterval_array_low, intervals));
  assert(is_sorted (omegas, node->a, node->b, Univinterval_array_high, intervals));
  assert(is_right_split (i, j, node->b, node->value, sigmas, intervals));

  for (lambda = node->a; lambda <= node->b; lambda++) {
    assert(Univinterval_is_contained(node->value,intervals,sigmas[lambda])
	   && Univinterval_is_contained(node->value,intervals,omegas[lambda]));
  }

  for (lambda = i; lambda <= node->a - 1; lambda++) {
    assert(Univinterval_is_contained(node->value,intervals,sigmas[lambda]) == false);
  }
  for (lambda = node->b + 1; lambda <= j; lambda++) {
    assert(Univinterval_is_contained(node->value,intervals,sigmas[lambda]) == false);
  }

  return true;
}



/* Refinement of node_make().
   Select proper value and split off "right of" triangles. */
static void 
node_select (int *index, Univcoord_T *value, int i, int j, 
	     int *sigmas, int *omegas, struct Univinterval_T *intervals) {
  int r = j - (j - i) / 3;
  Univcoord_T k = Univinterval_array_low(intervals,sigmas[r]);

  while ((r < j) && (Univinterval_array_low(intervals,sigmas[r + 1]) == k)) {
    r ++;
  }

  if (Univinterval_is_contained(k,intervals,sigmas[r]) == false) {
    /* adjust r to the left for "open" intervals */
    while ((r > i) && (Univinterval_is_contained(k,intervals,sigmas[r-1]) == false)) {
      r --;
      printf(" basic_iit: (-)\n");
    }
    if (Univinterval_is_contained(k,intervals,sigmas[r]) == false) {
      r --;
      printf(" basic_iit: [-]\n");
      assert(r == i - 1);
      printf(" basic_iit WARNING: empty NODE!?!\n");
    }
  }
  assert(is_right_split (i, j, r, k, sigmas, intervals));
  *index = r;
  *value = k;
  return;
}


/* Makes node out of sigmas[i..j], and recurses. */
static Node_T
Node_make (int *nnodes, int i, int j, int *sigmas, int *omegas, struct Univinterval_T *intervals) {
  Node_T node;
  int lambda, iota;
  int q, r;

  assert(is_valid_input (i, j, sigmas, omegas, intervals));
  if (i > j) {
    return (Node_T) NULL;
  } else {
    node = Node_new();
    *nnodes += 1;

    /* select value & get "right of" intervals sigma[r+1..j] */
    node_select(&r,&node->value,i,j,sigmas,omegas,intervals);

    /* mark "contains" intervals in sigma[i..r] to omega[i<q+1..r] */
    q = r;
    for (lambda = r; lambda >= i; lambda--) {
      if (Univinterval_is_contained(node->value,intervals,sigmas[lambda]) == true) {
	omegas[q] = sigmas[lambda];
	sigmas[lambda] = 0;
	q --;
      }
    }

    /* move remaining "left of" intervals from sigma[i..r] to sigma[i..q] */
    iota = i;
    for (lambda = i; lambda <= r; lambda++) {
      if (sigmas[lambda] != 0) {
	sigmas[iota] = sigmas[lambda];
	iota ++;
      }
    }
    assert(iota == q + 1);
    
    /* copy omega[q+1..r] back to sigma[q+1..r] & sort omega[q+1..r] */
    for (lambda = q+1; lambda <= r; lambda++) {
      sigmas[lambda] = omegas[lambda];
    }
    Univinterval_qsort_by_omega(omegas,q+1,r,intervals);
    node->a = q + 1;
    node->b = r;

    debug(printf(" NODE=%llu [%d..%d], left: %d, cont: %d, right: %d\n",
		 (unsigned long long) node->value, i, j, q - i + 1, r - q, j - r));

    assert(Node_is_valid_output (node, i, j, sigmas, omegas, intervals));

    /* recurse */
    node->left  = Node_make(&(*nnodes),i,q,sigmas,omegas,intervals);
    node->right = Node_make(&(*nnodes),r+1,j,sigmas,omegas,intervals);
    
    return node;
  }
}



static void 
Node_index (Node_T node, int *index) {
  if (node != NULL) {
    node->index = (*index)++;
    Node_index(node->left,&(*index));
    Node_index(node->right,&(*index));
  }
  return;
}


static int
IIT_count_nnodes (List_T intervallist) {
  int nnodes;
  Node_T root;
  int nintervals, i;
  List_T p;
  struct Univinterval_T *intervals;
  int *sigmas;
  int *omegas;

  if ((nintervals = List_length(intervallist)) == 0) {
    return 0;
  } else {
    intervals = (struct Univinterval_T *) CALLOC(nintervals,sizeof(struct Univinterval_T));
    for (p = intervallist, i = 0; i < nintervals; i++, p = List_next(p)) {
      memcpy(&(intervals[i]),List_head(p),sizeof(struct Univinterval_T));
    }

    /* IIT ordering of intervals */
    nnodes = 0;
    sigmas = (int *) CALLOC(nintervals+1,sizeof(int));
    for (i = 1; i <= nintervals; i++) {
      sigmas[i] = i;
    }

    /* Sort sigmas with respect to Univinterval_array_low */
    Univinterval_qsort_by_sigma(sigmas,1,nintervals,intervals);
    
    omegas = (int *) CALLOC(nintervals+1,sizeof(int));

    /* make first node, and recurse... */
    root = Node_make(&nnodes,1,nintervals,sigmas,omegas,intervals);

    Node_gc(&root);
    FREE(omegas);
    FREE(sigmas);
    FREE(intervals);

    return nnodes;
  }
}


static void
IIT_build_univ (Node_T *root, struct Univinterval_T **intervals, int **sigmas, int **omegas,
		List_T intervallist, int nintervals) {
  int index = 0;		/* Must be initialized to 0 */
  int i;
  List_T p;
  int nnodes;

  if (nintervals == 0) {
    *intervals = (struct Univinterval_T *) NULL;
  } else {
    *intervals = (struct Univinterval_T *) CALLOC(nintervals,sizeof(struct Univinterval_T));
    for (p = intervallist, i = 0; i < nintervals; i++, p = List_next(p)) {
      memcpy(&((*intervals)[i]),List_head(p),sizeof(struct Univinterval_T));
    }
  }

  /* IIT ordering of intervals */
  nnodes = 0;
  *sigmas = (int *) CALLOC(nintervals+1,sizeof(int));
  for (i = 1; i <= nintervals; i++) {
    (*sigmas)[i] = i;
  }

  /* Sort sigmas with respect to Univinterval_array_low */
  Univinterval_qsort_by_sigma(*sigmas,1,nintervals,*intervals);

  *omegas = (int *) CALLOC(nintervals+1,sizeof(int));

  /* make first node, and recurse... */
  *root = Node_make(&nnodes,1,nintervals,*sigmas,*omegas,*intervals);
  Node_index(*root,&index);

  return;
}


/************************************************************************
 *   Output procedures
 ************************************************************************/

/************************************************************************
 * File format for version 1 (universal coordinates)
 *
 *   total_nintervals: sizeof(int)  [so maximum is 2 billion]
 *      If >= 0, then using 4-byte quantities (unsigned int) for interval coordinates
 *      If < 0, then using 8-byte quantities for interval coordinates (generally only used for v1)
 *   ntypes: sizeof(int)
 *   nnodes: sizeof(int)
 *
 *   sigmas: (total_nintervals+1)*sizeof(int)
 *   omegas: (total_nintervals+1)*sizeof(int)
 *   nodes: nnodes*[uint + 4*int or uint8 + 4*int, depending on whether total_nintervals > 0]
 *
 *   intervals: total_nintervals*[uint+uint+int or uint8+uint8+int, depending on whether total_nintervals > 0]
 *
 *   typepointers: (ntypes+1)*sizeof(UINT4)
 *   types: ntypes*(variable length strings, including '\0')
 *
 *   labelorder: total_nintervals*sizeof(int);
 *   labelpointers: (total_nintervals+1)*sizeof(UINT4)  [changes to long unsigned int in v4 or if label_pointer_size = 8 in v5]
 *   labels: total_nintervals*(variable length strings, including '\0')
 *
 *   annotpointers: (total_nintervals+1)*sizeof(UINT4)  [changes to long unsigned int in v4 or if annot_pointer_size = 8 in v5]
 *   annotations: total_nintervals*(variable length strings, including '\0')
 *
 ************************************************************************/


/* For making labelorder */
struct Sortitem_T {
  int recno;
  char *label;
};

static int
Sortitem_cmp (const void *x, const void *y) {
  struct Sortitem_T *a = (struct Sortitem_T *) x;
  struct Sortitem_T *b = (struct Sortitem_T *) y;
  
  return strcmp(a->label,b->label);
}

static int *
get_labelorder (List_T divlist, Table_T labeltable, int total_nintervals) {
  int *labelorder, recno, i, k = 0;
  struct Sortitem_T *sortitems;
  char *divstring;
  List_T labellist, p;

  labelorder = (int *) CALLOC(total_nintervals,sizeof(int));
  sortitems = (struct Sortitem_T *) CALLOC(total_nintervals,sizeof(struct Sortitem_T));

  divstring = (char *) List_head(divlist);
  labellist = (List_T) Table_get(labeltable,(void *) divstring);
  recno = 0;
  for (p = labellist; p != NULL; p = List_next(p)) {
    sortitems[k].recno = recno;
    sortitems[k].label = (char *) List_head(p);
    k++;
    recno++;
  }

  qsort(sortitems,total_nintervals,sizeof(struct Sortitem_T),Sortitem_cmp);
  for (i = 0; i < total_nintervals; i++) {
    labelorder[i] = sortitems[i].recno;
  }
  FREE(sortitems);
  return labelorder;
}


/* Prints in DFS order */
static void
Node_fwrite (FILE *fp, Node_T node, bool coord_values_8p) {
  int leftindex, rightindex;
  UINT4 uint4;

  if (node != NULL) {
    if (node->left == NULL) {
      leftindex = -1;
    } else {
      leftindex = node->left->index;
    }

    if (node->right == NULL) {
      rightindex = -1;
    } else {
      rightindex = node->right->index;
    }

    if (coord_values_8p == true) {
      FWRITE_UINT8(node->value,fp);
    } else {
      uint4 = (UINT4) node->value;
      FWRITE_UINT(uint4,fp);
    }
    FWRITE_INT(node->a,fp);
    FWRITE_INT(node->b,fp);
    FWRITE_INT(leftindex,fp);
    FWRITE_INT(rightindex,fp);

    Node_fwrite(fp,node->left,coord_values_8p);
    Node_fwrite(fp,node->right,coord_values_8p);
  }
  return;
}


static void
IIT_write_univ_header (FILE *fp, int total_nintervals, int ntypes, int nnodes,
		       bool coord_values_8p) {
  int neg_total_nintervals;

  /* Writing before ...done */
  if (coord_values_8p == true) {
    fprintf(stderr,"coordinates require 8 bytes each...");
    neg_total_nintervals = -total_nintervals;
    FWRITE_INT(neg_total_nintervals,fp);
  } else {
    fprintf(stderr,"coordinates require 4 bytes each...");
    FWRITE_INT(total_nintervals,fp);
  }
  FWRITE_INT(ntypes,fp);
  FWRITE_INT(nnodes,fp);

  return;
}


static void
IIT_write_univ_footer (FILE *fp, List_T divlist, List_T typelist, Table_T intervaltable,
		       Table_T labeltable, Table_T annottable, int total_nintervals,
		       bool coord_values_8p, bool label_pointers_8p, bool annot_pointers_8p) {
  List_T intervallist, labellist, annotlist, p;
  Univinterval_T interval;
  UINT4 uint4;
#ifdef HAVE_64_BIT
  UINT8 pointer8 = 0LU;
#endif
  UINT4 pointer = 0U;
  int *labelorder;
  char *divstring, *typestring, *label, *annot;

#ifndef HAVE_64_BIT
  if (label_pointers_8p == true || annot_pointers_8p == true) {
    fprintf(stderr,"Machine is not 64-bit, so cannot write 8-byte pointers.\n");
    exit(9);
  }
#endif

  divstring = (char *) List_head(divlist);
  intervallist = (List_T) Table_get(intervaltable,(void *) divstring);
  if (coord_values_8p == true) {
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Univinterval_T) List_head(p);
      FWRITE_UINT8(interval->low,fp);
      FWRITE_UINT8(interval->high,fp);
      FWRITE_INT(interval->type,fp);
    }
  } else {
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Univinterval_T) List_head(p);
      uint4 = (UINT4) interval->low;
      FWRITE_UINT(uint4,fp);
      uint4 = (UINT4) interval->high;
      FWRITE_UINT(uint4,fp);
      FWRITE_INT(interval->type,fp);
    }
  }

  /* Write type pointers */
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
  for (p = typelist; p != NULL; p = List_next(p)) {
    typestring = (char *) List_head(p);
    if (typestring == NULL) {
      abort();			/* Even type 0 has an empty string */
    } else {
      pointer += (UINT4) strlen(typestring)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }
  }

  /* Write types */
  for (p = typelist; p != NULL; p = List_next(p)) {
    typestring = (char *) List_head(p);
    FWRITE_CHARS(typestring,strlen(typestring)+1,fp); /* Write '\0' */
  }      

  /* Write labelorder */
  if (total_nintervals > 0) {
    labelorder = get_labelorder(divlist,labeltable,total_nintervals);
    FWRITE_INTS(labelorder,total_nintervals,fp);
    FREE(labelorder);
  }

  /* Write label pointers */
#ifdef HAVE_64_BIT
  if (label_pointers_8p == true) {
    pointer8 = 0LU;
    FWRITE_UINT8(pointer8,fp);
  } else {
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
  }
#else
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
#endif

  divstring = (char *) List_head(divlist);
  labellist = (List_T) Table_get(labeltable,(void *) divstring);
  for (p = labellist; p != NULL; p = List_next(p)) {
    label = (char *) List_head(p);
#ifdef HAVE_64_BIT
    if (label_pointers_8p == true) {
      pointer8 += (UINT8) strlen(label)+1LU;	/* Add count for '\0' */
      FWRITE_UINT8(pointer8,fp);
    } else {
      pointer += (UINT4) strlen(label)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }
#else
    pointer += (UINT4) strlen(label)+1U;	/* Add count for '\0' */
    FWRITE_UINT(pointer,fp);
#endif
  }

  /* Write labels */
  divstring = (char *) List_head(divlist);
  labellist = (List_T) Table_get(labeltable,(void *) divstring);
  for (p = labellist; p != NULL; p = List_next(p)) {
    label = (char *) List_head(p);
    FWRITE_CHARS(label,strlen(label)+1,fp); /* Write '\0' */
  }

  /* Write annot pointers */
#ifdef HAVE_64_BIT
  if (annot_pointers_8p == true) {
    pointer8 = 0LU;
    FWRITE_UINT8(pointer8,fp);
  } else {
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
  }
#else
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
#endif

  divstring = (char *) List_head(divlist);
  annotlist = (List_T) Table_get(annottable,(void *) divstring);
  for (p = annotlist; p != NULL; p = List_next(p)) {
    annot = (char *) List_head(p);
#ifdef HAVE_64_BIT
    if (annot_pointers_8p == true) {
      pointer8 += (UINT8) strlen(annot)+1LU;	/* Add count for '\0' */
      FWRITE_UINT8(pointer8,fp);
    } else {
      pointer += (UINT4) strlen(annot)+1;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }
#else
    pointer += (UINT4) strlen(annot)+1;	/* Add count for '\0' */
    FWRITE_UINT(pointer,fp);
#endif
  }

  /* Write annotations */
  divstring = (char *) List_head(divlist);
  annotlist = (List_T) Table_get(annottable,(void *) divstring);
  for (p = annotlist; p != NULL; p = List_next(p)) {
    annot = (char *) List_head(p);
    FWRITE_CHARS(annot,strlen(annot)+1,fp); /* Write '\0' */
  }

  return;
}


/* If annotlist is NULL, X's are written */
void
IIT_write_univ (char *iitfile, List_T divlist, List_T typelist, Table_T intervaltable,
		Table_T labeltable, Table_T annottable,
		bool coord_values_8p, bool label_pointers_8p, bool annot_pointers_8p) {
  Node_T root;
  FILE *fp;
  List_T intervallist;
  char *divstring;
  int total_nintervals, nnodes;
  struct Univinterval_T *intervals;
  int *sigmas, *omegas;

  if ((fp = FOPEN_WRITE_BINARY(iitfile)) == NULL) {
    fprintf(stderr,"Error: can't open file %s\n",iitfile);
    exit(9);
  } else {
    assert(List_length(divlist) == 1);
    divstring = (char *) List_head(divlist);
    intervallist = (List_T) Table_get(intervaltable,(void *) divstring);
    total_nintervals = List_length(intervallist);
    nnodes = IIT_count_nnodes(intervallist);

#if 0
    if (total_nintervals == 0) {
      fprintf(stderr,"Error: No intervals were seen in input file\n");
      exit(9);
    }
#endif

    fprintf(stderr,"Writing IIT file header information...");
    IIT_write_univ_header(fp,total_nintervals,List_length(typelist),nnodes,coord_values_8p);
    fprintf(stderr,"done\n");

    divstring = (char *) List_head(divlist);
    intervallist = (List_T) Table_get(intervaltable,(void *) divstring);

    if (divstring[0] == '\0') {
      fprintf(stderr,"Processing null division/chromosome...sorting...");
    } else {
      fprintf(stderr,"Processing division/chromosome %s...sorting...",divstring);
    }
    IIT_build_univ(&root,&intervals,&sigmas,&omegas,intervallist,total_nintervals);

    fprintf(stderr,"writing...");
    FWRITE_INTS(sigmas,total_nintervals + 1,fp);
    FWRITE_INTS(omegas,total_nintervals + 1,fp);
    Node_fwrite(fp,root,coord_values_8p);
    fprintf(stderr,"done (%d intervals)\n",total_nintervals);

    Node_gc(&root);
    FREE(omegas);
    FREE(sigmas);
    FREE(intervals);

    fprintf(stderr,"Writing IIT file footer information...");
    IIT_write_univ_footer(fp,divlist,typelist,intervaltable,
			  labeltable,annottable,total_nintervals,
			  coord_values_8p,label_pointers_8p,annot_pointers_8p);
    fprintf(stderr,"done\n");

    fclose(fp);
    return;
  }
}



