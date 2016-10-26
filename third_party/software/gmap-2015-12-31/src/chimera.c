static char rcsid[] = "$Id: chimera.c 173190 2015-09-01 18:59:44Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chimera.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For sqrt */
#include <string.h>		/* For memset */
#include "mem.h"
#include "genomicpos.h"
#include "types.h"
#include "maxent.h"
#include "intron.h"
#include "comp.h"
#include "complement.h"


#define GBUFFERLEN 1024

/* Chimera assembly/bestpath matrix */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Chimera detection */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Exon-exon boundary detection */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* local_join_p and distant_join_p */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Chimera_bestpath */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


#define T Chimera_T
struct T {
  Stage3_T from;
  Stage3_T to;

  int chimerapos;
  int equivpos;

  int cdna_direction;
  int exonexonpos;

  char donor1;
  char donor2;
  char acceptor2;
  char acceptor1;

  bool donor_watsonp;
  bool acceptor_watsonp;

  double donor_prob;
  double acceptor_prob;
};


Stage3_T
Chimera_left_part (T this) {
  return this->from;
}

Stage3_T
Chimera_right_part (T this) {
  return this->to;
}


int 
Chimera_pos (T this) {
  return this->chimerapos;
}

int
Chimera_equivpos (T this) {
  return this->equivpos;
}

int
Chimera_cdna_direction (T this) {
  return this->cdna_direction;
}


void
Chimera_print_sam_tag (Filestring_T fp, T this, Univ_IIT_T chromosome_iit) {
  char donor_strand, acceptor_strand;
  char *donor_chr, *acceptor_chr;
  bool alloc1p, alloc2p;

  if (this->cdna_direction >= 0) {
    if (this->donor_watsonp == true) {
      donor_strand = '+';
    } else {
      donor_strand = '-';
    }
    if (this->acceptor_watsonp == true) {
      acceptor_strand = '+';
    } else {
      acceptor_strand = '-';
    }
  } else {
    if (this->donor_watsonp == true) {
      donor_strand = '-';
    } else {
      donor_strand = '+';
    }
    if (this->acceptor_watsonp == true) {
      acceptor_strand = '-';
    } else {
      acceptor_strand = '+';
    }
  }

  FPRINTF(fp,"%c%c-%c%c,%.2f,%.2f",
	  this->donor1,this->donor2,this->acceptor2,this->acceptor1,this->donor_prob,this->acceptor_prob);
  donor_chr = Univ_IIT_label(chromosome_iit,Stage3_chrnum(this->from),&alloc1p);
  acceptor_chr = Univ_IIT_label(chromosome_iit,Stage3_chrnum(this->to),&alloc2p);
  FPRINTF(fp,",%c%s@%u..%c%s@%u",
	  donor_strand,donor_chr,Stage3_chrend(this->from),
	  acceptor_strand,acceptor_chr,Stage3_chrstart(this->to));
  FPRINTF(fp,",%d..%d",this->chimerapos+1,this->equivpos+1);
  if (alloc2p == true) {
    FREE(acceptor_chr);
  }
  if (alloc1p == true) {
    FREE(donor_chr);
  }

  return;
}


double
Chimera_donor_prob (T this) {
  if (this->exonexonpos < 0) {
    return 0.0;
  } else {
    return this->donor_prob;
  }
}

double
Chimera_acceptor_prob (T this) {
  if (this->exonexonpos < 0) {
    return 0.0;
  } else {
    return this->acceptor_prob;
  }
}


T
Chimera_new (Stage3_T from, Stage3_T to, int chimerapos, int chimeraequivpos,
	     int exonexonpos, int cdna_direction,
	     char donor1, char donor2, char acceptor2, char acceptor1,
	     bool donor_watsonp, bool acceptor_watsonp,
	     double donor_prob, double acceptor_prob) {
  T new = (T) MALLOC(sizeof(*new));

  new->from = from;
  new->to = to;

  new->chimerapos = chimerapos;
  new->equivpos = chimeraequivpos;
  new->exonexonpos = exonexonpos;
  new->cdna_direction = cdna_direction;
  new->donor1 = donor1;
  new->donor2 = donor2;
  new->acceptor2 = acceptor2;
  new->acceptor1 = acceptor1;
  new->donor_watsonp = donor_watsonp;
  new->acceptor_watsonp = acceptor_watsonp;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}

void
Chimera_free (T *old) {
  FREE(*old);
  return;
}


void
Chimera_print (Filestring_T fp, T this) {
  if (this->exonexonpos > 0) {
    FPRINTF(fp," *** Possible chimera with exon-exon boundary");
    if (this->cdna_direction > 0) {
      FPRINTF(fp," (sense)");
    } else if (this->cdna_direction < 0) {
      FPRINTF(fp," (antisense)");
    }
    FPRINTF(fp," at %d (dinucl = %c%c-%c%c, donor_prob = %.3f, acceptor_prob = %.3f)",
	    this->exonexonpos+1,this->donor1,this->donor2,this->acceptor2,this->acceptor1,
	    this->donor_prob,this->acceptor_prob);
  } else if (this->equivpos == this->chimerapos) {
    FPRINTF(fp," *** Possible chimera with breakpoint at %d",this->chimerapos+1);
  } else {
    FPRINTF(fp," *** Possible chimera with breakpoint at %d..%d",this->chimerapos+1,this->equivpos+1);
  }

  return;
}


#define NPSEUDO 10.0

int
Chimera_alignment_break (int *newstart, int *newend, Stage3_T stage3, int queryntlength, double fthreshold) {
  int breakpoint;
  int *matchscores;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int start, end, pos, x = 0, y, n, x_left, y_left, x_right, y_right, n_left, n_right;
  double theta, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep, best_pos = -1, best_theta_left, best_theta_right;

  /*
  start = Sequence_trim_start(queryseq);
  end = Sequence_trim_end(queryseq);
  */

  start = Stage3_querystart(stage3);
  end = Stage3_queryend(stage3);

  if (queryntlength < MAX_QUERYLENGTH_STACK) {
    matchscores = (int *) CALLOCA(queryntlength,sizeof(int));
  } else {
    matchscores = (int *) CALLOC(queryntlength,sizeof(int));
  }
  Pair_matchscores(matchscores,Stage3_pairarray(stage3),Stage3_npairs(stage3),queryntlength);

  x = 0;
  for (pos = start; pos < end; pos++) {
    x += matchscores[pos];
  }
  n = end - start;
  y = n - x;

  /* when rss_sep == rss, fscore == 0 */
  min_rss_sep = rss = (double) x * (double) y/(double) n;
  if (rss == 0.0) {
    if (queryntlength < MAX_QUERYLENGTH_STACK) {
      FREEA(matchscores);
    } else {
      FREE(matchscores);
    }
    return 0;
  }

  theta = (double) x/(double) n;
  x_pseudo = NPSEUDO * theta;
  debug1(printf("%d %d %d %f\n",x,y,n,theta));
  
  x_left = y_left = n_left = 0;
  x_right = x;
  y_right = y;
  n_right = n;

  debug1(printf("%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","match","x.left","y.left","n.left","x.right","y.right","n.right",
		"theta.left","theta.right","rss.left","rss.right","fscore"));

  for (pos = start; pos < end-1; pos++) {
    if (matchscores[pos] == 1) {
      x_left++;
      x_right--;
    } else {
      y_left++;
      y_right--;
    }
    n_left++;
    n_right--;
    
    theta_left = ((double) x_left + x_pseudo)/((double) n_left + NPSEUDO);
    theta_right = ((double) x_right + x_pseudo)/((double) n_right + NPSEUDO);
    rss_left = x_left*(1.0-theta_left)*(1.0-theta_left) + y_left*theta_left*theta_left;
    rss_right = x_right*(1.0-theta_right)*(1.0-theta_right) + y_right*theta_right*theta_right;
    rss_sep = rss_left + rss_right;

    if (rss_sep == 0) {
      debug1(printf("%d %d %d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,matchscores[pos],x_left,y_left,n_left,x_right,y_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (n - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %d %d %d %f %f %f %f %f\n",
		    pos,matchscores[pos],x_left,y_left,n_left,x_right,y_right,n_right,
		    theta_left,theta_right,rss_left,rss_right,fscore);
	     );

      /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
	 1) is maximized when rss_sep is minimized */

      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	best_pos = pos;
	best_theta_left = theta_left;
	best_theta_right = theta_right;
      }
    }
  }
  if (queryntlength < MAX_QUERYLENGTH_STACK) {
    FREEA(matchscores);
  } else {
    FREE(matchscores);
  }

  fscore = ((double) (n - 2))*(rss - min_rss_sep)/min_rss_sep;
  if (fscore < fthreshold) {
    return 0;
  } else {
    breakpoint = best_pos;
    debug1(printf("at %d, fscore = %f\n",breakpoint,fscore));
    if (best_theta_left < best_theta_right) {
      /* trim left */
      *newstart = breakpoint;
      *newend = end;
      return breakpoint - start;
    } else {
      /* trim right */
      *newstart = start;
      *newend = breakpoint;
      return end - breakpoint;
    }
  }
}


bool
Chimera_local_join_p (Stage3_T from, Stage3_T to, int chimera_slop) {
  debug3(printf("? local_join_p from [%p] %d..%d (%u..%u) -> to [%p] %d..%d (%u..%u) => ",
		from,Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		to,Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug3(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug3(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if (Stage3_chrnum(from) != Stage3_chrnum(to)) {
    debug3(printf("false, because different chromosomes\n"));
    return false;

  } else if (Stage3_watsonp(from) != Stage3_watsonp(to)) {
    debug3(printf("false, because different strands\n"));
    return false;

  } else if (Stage3_querystart(from) >= Stage3_querystart(to) &&
	     Stage3_queryend(from) <= Stage3_queryend(to)) {
    debug3(printf("false, because from %d..%d is subsumed by to %d..%d\n",
		  Stage3_querystart(from),Stage3_queryend(from),
		  Stage3_querystart(to),Stage3_queryend(to)));
    return false;

  } else if (Stage3_querystart(to) >= Stage3_querystart(from) &&
	     Stage3_queryend(to) <= Stage3_queryend(from)) {
    debug3(printf("false, because to %d..%d is subsumed by from %d..%d\n",
		  Stage3_querystart(to),Stage3_queryend(to),
		  Stage3_querystart(from),Stage3_queryend(from)));
    return false;

  } else if (Stage3_queryend(from) - Stage3_querystart(to) > chimera_slop ||
	     Stage3_querystart(to) - Stage3_queryend(from) > chimera_slop) {
    debug3(printf("false, because %d - %d > chimera_slop %d or %d - %d > chimera_slop %d\n",
		  Stage3_queryend(from),Stage3_querystart(to),chimera_slop,
		  Stage3_querystart(to),Stage3_queryend(from),chimera_slop));
    return false;

  } else {
    if (Stage3_watsonp(from) == true) {
      if (Stage3_genomicend(from) > Stage3_genomicstart(to) + chimera_slop) {
	debug3(printf("false, because genomic %u > %u + %d\n",
		      Stage3_genomicend(from),Stage3_genomicstart(to),chimera_slop));
	return false;
      } else {
	return true;
      }

    } else {
      if (Stage3_genomicend(from) + chimera_slop < Stage3_genomicstart(to)) {
	debug3(printf("false, because genomic %u + %d < %u\n",
		      Stage3_genomicend(from),chimera_slop,Stage3_genomicstart(to)));
	return false;
      } else {
	return true;
      }
      
    }

  }
}


bool
Chimera_distant_join_p (Stage3_T from, Stage3_T to, int chimera_slop) {
  debug3(printf("? chimeric_join_p from %d..%d (%u..%u) -> to %d..%d (%u..%u) => ",
		Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug3(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug3(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if (Stage3_querystart(from) >= Stage3_querystart(to) &&
      Stage3_queryend(from) <= Stage3_queryend(to)) {
    debug3(printf("false, because from %d..%d is subsumed by to %d..%d\n",
		  Stage3_querystart(from),Stage3_queryend(from),
		  Stage3_querystart(to),Stage3_queryend(to)));
    return false;
  } else if (Stage3_querystart(to) >= Stage3_querystart(from) &&
	     Stage3_queryend(to) <= Stage3_queryend(from)) {
    debug3(printf("false, because to %d..%d is subsumed by from %d..%d\n",
		  Stage3_querystart(to),Stage3_queryend(to),
		  Stage3_querystart(from),Stage3_queryend(from)));
    return false;
  } else if (Stage3_queryend(from) - Stage3_querystart(to) <= chimera_slop &&
	     Stage3_querystart(to) - Stage3_queryend(from) <= chimera_slop) {
    debug3(printf("true, because %d - %d <= %d and %d - %d <= %d\n",
		  Stage3_queryend(from),Stage3_querystart(to),chimera_slop,
		  Stage3_querystart(to),Stage3_queryend(from),chimera_slop));
    return true;
  } else {
    debug3(printf(" %d and %d not within chimera_slop %d",
		  Stage3_queryend(from) - Stage3_querystart(to),Stage3_querystart(to) - Stage3_queryend(from),
		  chimera_slop));
    debug3(printf("false\n"));
    return false;
  }
}


#define NEG_INFINITY -1000000
#define PRE_EXTENSION_SLOP 6

bool
Chimera_bestpath (int *five_score, int *three_score, int *chimerapos, int *chimeraequivpos, int *bestfrom, int *bestto, 
		  Stage3_T *stage3array_sub1, int npaths_sub1, Stage3_T *stage3array_sub2, int npaths_sub2, 
		  int queryntlength, int chimera_slop, bool localp) {
  int **matrix_sub1, **matrix_sub2, *from, *to, *bestscoreatpos, i, j, pos, score, 
    bestscore = NEG_INFINITY;
  bool **gapp_sub1, **gapp_sub2;
  bool foundp = false;
  
  debug4(printf("Chimera_bestpath called\n"));

  from = (int *) CALLOC(queryntlength,sizeof(int));
  to = (int *) CALLOC(queryntlength,sizeof(int));
  bestscoreatpos = (int *) CALLOC(queryntlength,sizeof(int));

  matrix_sub1 = (int **) CALLOC(npaths_sub1,sizeof(int *));
  gapp_sub1 = (bool **) CALLOC(npaths_sub1,sizeof(bool *));
  debug4(printf("sub1:"));
  for (i = 0; i < npaths_sub1; i++) {
    debug4(printf(" %p",stage3array_sub1[i]));
    matrix_sub1[i] = (int *) CALLOC(queryntlength,sizeof(int));
    gapp_sub1[i] = (bool *) CALLOC(queryntlength,sizeof(bool));
    debug4(Pair_dump_array(Stage3_pairarray(stage3array_sub1[i]),Stage3_npairs(stage3array_sub1[i]),true));
    /* Allow pre_extension_slop, in case the parts need extensions to merge */
    Pair_pathscores(gapp_sub1[i],matrix_sub1[i],Stage3_pairarray(stage3array_sub1[i]),
		    Stage3_npairs(stage3array_sub1[i]),Stage3_cdna_direction(stage3array_sub1[i]),
		    queryntlength,FIVE,PRE_EXTENSION_SLOP);
  }
  debug4(printf("\n"));

  debug4(printf("sub2:"));
  matrix_sub2 = (int **) CALLOC(npaths_sub2,sizeof(int *));
  gapp_sub2 = (bool **) CALLOC(npaths_sub2,sizeof(bool *));
  for (i = 0; i < npaths_sub2; i++) {
    debug4(printf(" %p",stage3array_sub2[i]));
    matrix_sub2[i] = (int *) CALLOC(queryntlength,sizeof(int));
    gapp_sub2[i] = (bool *) CALLOC(queryntlength,sizeof(bool));
    debug4(Pair_dump_array(Stage3_pairarray(stage3array_sub2[i]),Stage3_npairs(stage3array_sub2[i]),true));
    /* Allow pre_extension_slop, in case the parts need extensions to merge */
    Pair_pathscores(gapp_sub2[i],matrix_sub2[i],Stage3_pairarray(stage3array_sub2[i]),
		    Stage3_npairs(stage3array_sub2[i]),Stage3_cdna_direction(stage3array_sub2[i]),
		    queryntlength,THREE,PRE_EXTENSION_SLOP);
  }
  debug4(printf("\n"));

  for (pos = 0; pos < queryntlength; pos++) {
    bestscoreatpos[pos] = NEG_INFINITY;
  }
  debug4(printf("npaths_sub1 = %d, npaths_sub2 = %d\n",npaths_sub1,npaths_sub2));
  for (i = 0; i < npaths_sub1; i++) {
    for (j = 0; j < npaths_sub2; j++) {
      if (stage3array_sub1[i] == stage3array_sub2[j]) {
	/* Same stage3 object, so not joinable */
      } else if (localp == true && Chimera_local_join_p(stage3array_sub1[i],stage3array_sub2[j],chimera_slop) == false) {
	/* Not joinable */
      } else {
	for (pos = 0; pos < queryntlength - 1; pos++) {
	  debug4(printf("pos %d, gapp %d and %d\n",pos,gapp_sub1[i][pos],gapp_sub2[j][pos]));
	  if (gapp_sub1[i][pos] == false && gapp_sub2[j][pos+1] == false) {
#if 0
	    score = matrix_sub2[j][queryntlength-1] - matrix_sub2[j][pos] + matrix_sub1[i][pos] /* - 0 */;
#else
	    /* For new Pair_pairscores computation */
	    score = matrix_sub1[i][pos] + matrix_sub2[j][pos];
#endif
	    debug4(printf("score %d\n",score));
	    if (score > bestscoreatpos[pos]) {
	      bestscoreatpos[pos] = score;
	      from[pos] = i;
	      to[pos] = j;
	    }
	  }
	}
      }
    }
  }

  for (pos = 0; pos < queryntlength - 1; pos++) {
    if (bestscoreatpos[pos] > bestscore) {
      bestscore = bestscoreatpos[pos];
      *chimerapos = *chimeraequivpos = pos;
      *bestfrom = from[pos];
      *bestto = to[pos];
      foundp = true;
    } else if (bestscoreatpos[pos] == bestscore) {
      *chimeraequivpos = pos;
    }
  }

  if (foundp == true) {
    *five_score = matrix_sub1[*bestfrom][*chimerapos] /* - 0 */;
#if 0
    *three_score = matrix_sub2[*bestto][queryntlength-1] - matrix_sub2[*bestto][*chimerapos];
#else
    *three_score = matrix_sub2[*bestto][*chimerapos];
#endif

    debug4(
	  for (pos = 0; pos < queryntlength - 1; pos++) {
	    printf("%d:",pos);
	    for (i = 0; i < npaths_sub1; i++) {
	      printf("\t%d",matrix_sub1[i][pos]);
	      if (gapp_sub1[i][pos] == true) {
		printf("X");
	      }
	    }
	    printf("\t|");
	    for (i = 0; i < npaths_sub2; i++) {
	      printf("\t%d",matrix_sub2[i][pos]);
	      if (gapp_sub2[i][pos] == true) {
		printf("X");
	      }
	    }
	    printf("\t||");
	    printf("%d (%d->%d)",bestscoreatpos[pos],from[pos],to[pos]);
	    if (pos >= *chimerapos && pos <= *chimeraequivpos) {
	      printf(" ** ");
	    }
	    printf("\n");
	  }
	  printf("From path %d to path %d at pos %d..%d.  5 score = %d, 3 score = %d\n",
		 *bestfrom,*bestto,*chimerapos,*chimeraequivpos,*five_score,*three_score);
	  fflush(stdout);
	  );
  }
  
  for (i = 0; i < npaths_sub2; i++) {
    FREE(gapp_sub2[i]);
    FREE(matrix_sub2[i]);
  }
  FREE(gapp_sub2);
  FREE(matrix_sub2);

  for (i = 0; i < npaths_sub1; i++) {
    FREE(gapp_sub1[i]);
    FREE(matrix_sub1[i]);
  }
  FREE(gapp_sub1);
  FREE(matrix_sub1);

  FREE(bestscoreatpos);
  FREE(to);
  FREE(from);

  debug4(printf("Chimera_bestpath returning %d\n",foundp));
  return foundp;
}

static char *complCode = COMPLEMENT_UC;

/* Modeled after Chimera_bestpath */
/* Called if Chimera_find_exonexon fails */
int
Chimera_find_breakpoint (int *chimeraequivpos, char *donor1, char *donor2, char *acceptor2, char *acceptor1,
			 Stage3_T left_part, Stage3_T right_part, int queryntlength, Genome_T genome,
			 Chrpos_T left_chrlength, Chrpos_T right_chrlength) {
  int chimerapos = 0, breakpoint;
  int *matrix_sub1, *matrix_sub2, pos, score, bestscore;
  bool *gapp_sub1, *gapp_sub2;
  Univcoord_T left;

  /* Don't allow pre_extension_slop here, because the ends have already been extended */
  matrix_sub1 = (int *) CALLOC(queryntlength,sizeof(int));
  gapp_sub1 = (bool *) CALLOC(queryntlength,sizeof(bool));
  debug4(Pair_dump_array(Stage3_pairarray(left_part),Stage3_npairs(left_part),true));
  Pair_pathscores(gapp_sub1,matrix_sub1,Stage3_pairarray(left_part),Stage3_npairs(left_part),
		  Stage3_cdna_direction(left_part),queryntlength,FIVE,/*pre_extension_slop*/0);

  matrix_sub2 = (int *) CALLOC(queryntlength,sizeof(int));
  gapp_sub2 = (bool *) CALLOC(queryntlength,sizeof(bool));
  debug4(Pair_dump_array(Stage3_pairarray(right_part),Stage3_npairs(right_part),true));
  Pair_pathscores(gapp_sub2,matrix_sub2,Stage3_pairarray(right_part),Stage3_npairs(right_part),
		  Stage3_cdna_direction(right_part),queryntlength,THREE,/*pre_extension_slop*/0);


  bestscore = -100000;
  for (pos = 0; pos < queryntlength - 1; pos++) {
    debug(
	  printf("%d:",pos);
	  printf("\t%d",matrix_sub1[pos]);
	  if (gapp_sub1[pos] == true) {
	    printf("X");
	  }
	  printf("\t|");
	  printf("\t%d",matrix_sub2[pos]);
	  if (gapp_sub2[pos] == true) {
	    printf("X");
	  }
	  printf("\t||");
	  );

    if (gapp_sub1[pos] == false) {
      if (gapp_sub2[pos+1] == false) {
	/* Check for the same stage3 object on both lists */
#if 0
	/* ? Old formula for use before Pair_pathscores had cdnaend argument */
	score = matrix_sub2[queryntlength-1] - matrix_sub2[pos] + matrix_sub1[pos] /* - 0 */;
#else
	score = matrix_sub1[pos] + matrix_sub2[pos+1];
#endif

	if (score > bestscore) {
	  bestscore = score;
	  chimerapos = *chimeraequivpos = pos;
	} else if (score == bestscore) {
	  *chimeraequivpos = pos;
	}

	debug(
	      printf("%d = %d + %d",score,matrix_sub1[pos],matrix_sub2[pos+1]);
	      if (pos >= chimerapos && pos <= *chimeraequivpos) {
		printf(" ** chimerapos %d, chimeraequivpos %d",chimerapos,*chimeraequivpos);
	      }
	      );

      }
    }
    debug(printf("\n"));
  }
  debug(printf("chimerapos %d, chimeraequivpos %d\n",chimerapos,*chimeraequivpos));

#if 0
  *five_score = matrix_sub1[*chimerapos] /* - 0 */;
  *three_score = matrix_sub2[queryntlength-1] - matrix_sub2[*chimerapos];
#endif

  FREE(gapp_sub2);
  FREE(matrix_sub2);

  FREE(gapp_sub1);
  FREE(matrix_sub1);

  if (chimerapos == 0) {
    /* Never found a breakpoint */
    return -1;
  } else {
    breakpoint = (chimerapos + (*chimeraequivpos))/2;

    if (Stage3_watsonp(left_part) == true) {
      if ((left = Stage3_genomicpos(left_part,breakpoint,/*headp*/false)) >= left_chrlength - 2) {
	debug(printf("left %u >= left_chrlength %u - 2, so not finding donor dinucleotides\n",left,left_chrlength));
	*donor1 = *donor2 = 'N';
      } else {
	debug(printf("left %u < left_chrlength %u - 2, so okay\n",left,left_chrlength));
	*donor1 = Genome_get_char(genome,left+1);
	*donor2 = Genome_get_char(genome,left+2);
      }
    } else {
      if ((left = Stage3_genomicpos(left_part,breakpoint,/*headp*/false)) < 2) {
	debug(printf("left %u < 2, so not finding donor dinucleotides\n",left));
	*donor1 = *donor2 = 'N';
      } else {
	debug(printf("left %u >= 2, so okay\n",left));
	*donor1 = complCode[(int) Genome_get_char(genome,left-1)];
	*donor2 = complCode[(int) Genome_get_char(genome,left-2)];
      }
    }

    if (Stage3_watsonp(right_part) == true) {
      if ((left = Stage3_genomicpos(right_part,breakpoint+1,/*headp*/true)) < 2) {
	debug(printf("left %u < 2, so not finding acceptor dinucleotides\n",left));
	*acceptor1 = *acceptor2 = 'N';
      } else {
	debug(printf("left %u >= 2, so okay\n",left));
	*acceptor2 = Genome_get_char(genome,left-2);
	*acceptor1 = Genome_get_char(genome,left-1);
      }
    } else {
      if ((left = Stage3_genomicpos(right_part,breakpoint+1,/*headp*/true)) >= right_chrlength - 2) {
	debug(printf("left %u >= right_chrlength %u - 2, so not finding acceptor dinucleotides\n",left,right_chrlength));
	*acceptor1 = *acceptor2 = 'N';
      } else {
	debug(printf("left %u <right_chrlength %u - 2, so okay\n",left,right_chrlength));
	*acceptor2 = complCode[(int) Genome_get_char(genome,left+2)];
	*acceptor1 = complCode[(int) Genome_get_char(genome,left+1)];
      }
    }

    return chimerapos;
  }
}


static double
find_exonexon_fwd (int *exonexonpos, char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		   char *comp, bool *donor_watsonp, bool *acceptor_watsonp, double *donor_prob, double *acceptor_prob,
		   Stage3_T left_part, Stage3_T right_part, Genome_T genome, Genome_T genomealt,
		   Univ_IIT_T chromosome_iit, int breakpoint_start, int breakpoint_end) {
  Sequence_T donor_genomicseg, acceptor_genomicseg, donor_genomicalt, acceptor_genomicalt;
  char *donor_ptr, *acceptor_ptr, *donor_altptr, *acceptor_altptr;
  int i, j;
  Univcoord_T left;
  int donor_length, acceptor_length;
  bool revcomp;
  char left1, left2, right2, right1;
  char left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;
  double donor_prob_1, acceptor_prob_1, donor_altprob_1, acceptor_altprob_1, bestproduct = 0.0, product;

  *exonexonpos = -1;

  donor_length = breakpoint_end - breakpoint_start + DONOR_MODEL_LEFT_MARGIN + DONOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(left_part,breakpoint_start,/*headp*/false);
  if ((*donor_watsonp = Stage3_watsonp(left_part)) == true) {
    left -= DONOR_MODEL_LEFT_MARGIN;
    revcomp = false;
  } else {
    left += DONOR_MODEL_LEFT_MARGIN;
    left -= donor_length;
    revcomp = true;
  }

  debug2(printf("Getting donor at left %u\n",left));
  donor_genomicseg = Genome_get_segment(genome,left,donor_length+1,chromosome_iit,revcomp);
  donor_ptr = Sequence_fullpointer(donor_genomicseg);
  donor_genomicalt = Genome_get_segment_alt(genomealt,left,donor_length+1,chromosome_iit,revcomp);
  donor_altptr = Sequence_fullpointer(donor_genomicalt);

  debug2(
	 for (i = 0; i < ACCEPTOR_MODEL_LEFT_MARGIN - DONOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(stdout,donor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  acceptor_length = breakpoint_end - breakpoint_start + ACCEPTOR_MODEL_LEFT_MARGIN + ACCEPTOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(right_part,breakpoint_end+1,/*headp*/true);
  if ((*acceptor_watsonp = Stage3_watsonp(right_part)) == true) {
    left += ACCEPTOR_MODEL_RIGHT_MARGIN;
    left -= acceptor_length;
    revcomp = false;
  } else {
    left -= ACCEPTOR_MODEL_RIGHT_MARGIN;
    revcomp = true;
  }

  debug2(printf("Getting acceptor at left %u\n",left));
  acceptor_genomicseg = Genome_get_segment(genome,left,acceptor_length+1,chromosome_iit,revcomp);
  acceptor_ptr = Sequence_fullpointer(acceptor_genomicseg);
  acceptor_genomicalt = Genome_get_segment(genomealt,left,acceptor_length+1,chromosome_iit,revcomp);
  acceptor_altptr = Sequence_fullpointer(acceptor_genomicalt);
  debug2(
	 for (i = 0; i < DONOR_MODEL_LEFT_MARGIN - ACCEPTOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(stdout,acceptor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  *donor_prob = 0.0;
  *acceptor_prob = 0.0;
  for (i = DONOR_MODEL_LEFT_MARGIN, j = ACCEPTOR_MODEL_LEFT_MARGIN; 
       i <= donor_length - DONOR_MODEL_RIGHT_MARGIN && 
	 j <= acceptor_length - ACCEPTOR_MODEL_RIGHT_MARGIN;
       i++, j++) {

    left1 = donor_ptr[i+1];
    left2 = donor_ptr[i+2];
    right2 = acceptor_ptr[j-2];
    right1 = acceptor_ptr[j-1];

    left1_alt = donor_altptr[i+1];
    left2_alt = donor_altptr[i+2];
    right2_alt = acceptor_altptr[j-2];
    right1_alt = acceptor_altptr[j-1];

    debug2(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     /*cdna_direction*/+1);
    debug2(printf("  Introntype is %s\n",Intron_type_string(introntype)));

    donor_prob_1 = Maxent_donor_prob(&(donor_ptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
    acceptor_prob_1 = Maxent_acceptor_prob(&(acceptor_ptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));
    donor_altprob_1 = Maxent_donor_prob(&(donor_altptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
    acceptor_altprob_1 = Maxent_acceptor_prob(&(acceptor_altptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));

    debug2(printf("%d %c%c %c%c %.2f %.2f\n",
		  breakpoint_start - DONOR_MODEL_LEFT_MARGIN + i,
		  donor_ptr[i+1],donor_ptr[i+2],acceptor_ptr[j-2],acceptor_ptr[j-1],
		  donor_prob_1,acceptor_prob_1));

    if (donor_prob_1 < 0.50 && acceptor_prob_1 < 0.50 && donor_altprob_1 < 0.50 && acceptor_altprob_1 < 0.50) {
      /* Skip */
    } else if (introntype != NONINTRON || donor_prob_1 > 0.90 || acceptor_prob_1 > 0.90 || donor_altprob_1 > 0.90 || acceptor_altprob_1 > 0.90) {
      if ((product = donor_prob_1*acceptor_prob_1) > bestproduct) {
	bestproduct = product;
	*donor1 = donor_ptr[i+1];
	*donor2 = donor_ptr[i+2];
	*acceptor2 = acceptor_ptr[j-2];
	*acceptor1 = acceptor_ptr[j-1];
	if (donor_prob_1 >= donor_altprob_1) {
	  *donor_prob = donor_prob_1;
	} else {
	  *donor_prob = donor_altprob_1;
	}
	if (acceptor_prob_1 >= acceptor_altprob_1) {
	  *acceptor_prob = acceptor_prob_1;
	} else {
	  *acceptor_prob = acceptor_altprob_1;
	}
	*exonexonpos = breakpoint_start - DONOR_MODEL_LEFT_MARGIN + i;

	switch (introntype) {
	case GTAG_FWD: *comp = FWD_CANONICAL_INTRON_COMP; break;
	case GCAG_FWD: *comp = FWD_GCAG_INTRON_COMP; break;
	case ATAC_FWD: *comp = FWD_ATAC_INTRON_COMP; break;
	default: *comp = NONINTRON_COMP; break;
	}

      }
    }
  }

  Sequence_free(&acceptor_genomicalt);
  Sequence_free(&donor_genomicalt);
  Sequence_free(&acceptor_genomicseg);
  Sequence_free(&donor_genomicseg);

  return bestproduct;
}

static double
find_exonexon_rev (int *exonexonpos, char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		   char *comp, bool *donor_watsonp, bool *acceptor_watsonp, double *donor_prob, double *acceptor_prob,
		   Stage3_T left_part, Stage3_T right_part, Genome_T genome, Genome_T genomealt,
		   Univ_IIT_T chromosome_iit, int breakpoint_start, int breakpoint_end) {
  Sequence_T donor_genomicseg, acceptor_genomicseg, donor_genomicalt, acceptor_genomicalt;
  char *donor_ptr, *acceptor_ptr, *donor_altptr, *acceptor_altptr;
  int i, j;
  Univcoord_T left;
  int donor_length, acceptor_length;
  bool revcomp;
  char left1, left2, right2, right1;
  char left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;
  double donor_prob_1, acceptor_prob_1, donor_altprob_1, acceptor_altprob_1, bestproduct = 0.0, product;

  *exonexonpos = -1;

  donor_length = breakpoint_end - breakpoint_start + DONOR_MODEL_LEFT_MARGIN + DONOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(right_part,breakpoint_end+1,/*headp*/true);
  if ((*donor_watsonp = Stage3_watsonp(right_part)) == true) {
    left += DONOR_MODEL_LEFT_MARGIN;
    left -= donor_length;
    revcomp = true;
  } else {
    left -= DONOR_MODEL_LEFT_MARGIN;
    revcomp = false;
  }

  debug2(printf("Getting donor at left %u\n",left));
  donor_genomicseg = Genome_get_segment(genome,left,donor_length+1,chromosome_iit,revcomp);
  donor_ptr = Sequence_fullpointer(donor_genomicseg);
  donor_genomicalt = Genome_get_segment(genomealt,left,donor_length+1,chromosome_iit,revcomp);
  donor_altptr = Sequence_fullpointer(donor_genomicalt);
  debug2(
	 for (i = 0; i < ACCEPTOR_MODEL_LEFT_MARGIN - DONOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(stdout,donor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  acceptor_length = breakpoint_end - breakpoint_start + ACCEPTOR_MODEL_LEFT_MARGIN + ACCEPTOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(left_part,breakpoint_start,/*headp*/false);
  if ((*acceptor_watsonp = Stage3_watsonp(left_part)) == true) {
    left -= ACCEPTOR_MODEL_RIGHT_MARGIN;
    revcomp = true;
  } else {
    left += ACCEPTOR_MODEL_RIGHT_MARGIN;
    left -= acceptor_length;
    revcomp = false;
  }

  debug2(printf("Getting acceptor at left %u\n",left));
  acceptor_genomicseg = Genome_get_segment(genome,left,acceptor_length+1,chromosome_iit,revcomp);
  acceptor_ptr = Sequence_fullpointer(acceptor_genomicseg);
  acceptor_genomicalt = Genome_get_segment(genomealt,left,acceptor_length+1,chromosome_iit,revcomp);
  acceptor_altptr = Sequence_fullpointer(acceptor_genomicalt);
  debug2(
	 for (i = 0; i < DONOR_MODEL_LEFT_MARGIN - ACCEPTOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(stdout,acceptor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  *donor_prob = 0.0;
  *acceptor_prob = 0.0;
  for (i = DONOR_MODEL_LEFT_MARGIN, j = ACCEPTOR_MODEL_LEFT_MARGIN; 
       i <= donor_length - DONOR_MODEL_RIGHT_MARGIN && 
	 j <= acceptor_length - ACCEPTOR_MODEL_RIGHT_MARGIN;
       i++, j++) {

    left1 = donor_ptr[i+1];
    left2 = donor_ptr[i+2];
    right2 = acceptor_ptr[j-2];
    right1 = acceptor_ptr[j-1];

    left1_alt = donor_altptr[i+1];
    left2_alt = donor_altptr[i+2];
    right2_alt = acceptor_altptr[j-2];
    right1_alt = acceptor_altptr[j-1];

    /* Use cdna_direction == +1, because revcomp already applied */
    debug2(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     /*cdna_direction*/+1);
    debug2(printf("  Introntype is %s\n",Intron_type_string(introntype)));

    donor_prob_1 = Maxent_donor_prob(&(donor_ptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
    acceptor_prob_1 = Maxent_acceptor_prob(&(acceptor_ptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));
    donor_altprob_1 = Maxent_donor_prob(&(donor_altptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
    acceptor_altprob_1 = Maxent_acceptor_prob(&(acceptor_altptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));

    debug2(printf("%d %c%c %c%c %.2f %.2f\n",
		  breakpoint_end + DONOR_MODEL_LEFT_MARGIN - i,
		  donor_ptr[i+1],donor_ptr[i+2],acceptor_ptr[j-2],acceptor_ptr[j-1],
		  donor_prob_1,acceptor_prob_1));

    if (donor_prob_1 < 0.50 && acceptor_prob_1 < 0.50 && donor_altprob_1 < 0.50 && acceptor_altprob_1 < 0.50) {
      /* Skip */
    } else if (introntype != NONINTRON || donor_prob_1 > 0.90 || acceptor_prob_1 > 0.90 || donor_altprob_1 > 0.90 || acceptor_altprob_1 > 0.90) {
      if ((product = donor_prob_1*acceptor_prob_1) > bestproduct) {
	bestproduct = product;
	*donor1 = donor_ptr[i+1];
	*donor2 = donor_ptr[i+2];
	*acceptor2 = acceptor_ptr[j-2];
	*acceptor1 = acceptor_ptr[j-1];
	if (donor_prob_1 >= donor_altprob_1) {
	  *donor_prob = donor_prob_1;
	} else {
	  *donor_prob = donor_altprob_1;
	}
	if (acceptor_prob_1 >= acceptor_altprob_1) {
	  *acceptor_prob = acceptor_prob_1;
	} else {
	  *acceptor_prob = acceptor_altprob_1;
	}
	*exonexonpos = breakpoint_end + DONOR_MODEL_LEFT_MARGIN - i;

	/* Have to look for forward intron types, but return the revcomp comp */
	switch (introntype) {
	case GTAG_FWD: *comp = REV_CANONICAL_INTRON_COMP; break;
	case GCAG_FWD: *comp = REV_GCAG_INTRON_COMP; break;
	case ATAC_FWD: *comp = REV_ATAC_INTRON_COMP; break;
	default: *comp = NONINTRON_COMP; break;
	}
      }
    }
  }

  Sequence_free(&acceptor_genomicalt);
  Sequence_free(&donor_genomicalt);
  Sequence_free(&acceptor_genomicseg);
  Sequence_free(&donor_genomicseg);

  return bestproduct;
}


int
Chimera_find_exonexon (int *found_cdna_direction, int *try_cdna_direction,
		       char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		       char *comp, bool *donor_watsonp, bool *acceptor_watsonp, double *donor_prob, double *acceptor_prob,
		       Stage3_T left_part, Stage3_T right_part, Genome_T genome, Genome_T genomealt,
		       Univ_IIT_T chromosome_iit, int breakpoint_start, int breakpoint_end) {
  int exonexonpos_fwd, exonexonpos_rev;
  char donor1_fwd, donor2_fwd, acceptor2_fwd, acceptor1_fwd,
    donor1_rev, donor2_rev, acceptor2_rev, acceptor1_rev;
  char comp_fwd, comp_rev;
  bool donor_watsonp_fwd, acceptor_watsonp_fwd, donor_watsonp_rev, acceptor_watsonp_rev;
  double bestproduct_fwd, bestproduct_rev, donor_prob_fwd, donor_prob_rev, acceptor_prob_fwd, acceptor_prob_rev;
  int left_cdna_direction, right_cdna_direction;

  debug2(printf("Starting Chimera_find_exonexon with breakpoint %d..%d\n",breakpoint_start,breakpoint_end));
  debug2(printf("left part covers query %d to %d\n",Stage3_querystart(left_part),Stage3_queryend(left_part)));
  debug2(printf("right part covers query %d to %d\n",Stage3_querystart(right_part),Stage3_queryend(right_part)));

#if 0
  if (Stage3_queryend(left_part) < Stage3_querystart(right_part)) {
    breakpoint_start = Stage3_queryend(left_part);
    breakpoint_end = Stage3_querystart(right_part);
  } else {
    breakpoint_start = Stage3_querystart(right_part);
    breakpoint_end = Stage3_queryend(left_part);
  }
  breakpoint_start -= 10;
  if (breakpoint_start < 0) {
    breakpoint_start = 0;
  }
  breakpoint_end += 10;
  if (breakpoint_end > querylength-1) {
    breakpoint_end = querylength-1;
  }
#endif

  if (breakpoint_end < breakpoint_start) {
    debug2(printf("Breakpoints do not make sense, so not computing\n"));
    *found_cdna_direction = *try_cdna_direction = 0;
    return -1;
  }

  debug2(printf("Starting search for exon-exon boundary at breakpoint_start %d to breakpoint_end %d\n",
		breakpoint_start,breakpoint_end));

  left_cdna_direction = Stage3_cdna_direction(left_part);
  right_cdna_direction = Stage3_cdna_direction(right_part);

  if (left_cdna_direction == 0 && right_cdna_direction == 0) {
    *try_cdna_direction = 0;
  } else if (left_cdna_direction >= 0 && right_cdna_direction >= 0) {
    *try_cdna_direction = +1;
  } else if (left_cdna_direction <= 0 && right_cdna_direction <= 0) {
    *try_cdna_direction = -1;
  } else {
    *try_cdna_direction = 0;
  }

  if (*try_cdna_direction == +1) {
    *found_cdna_direction = +1;
    bestproduct_fwd = find_exonexon_fwd(&exonexonpos_fwd,&donor1_fwd,&donor2_fwd,&acceptor2_fwd,&acceptor1_fwd,
					&comp_fwd,&donor_watsonp_fwd,&acceptor_watsonp_fwd,&donor_prob_fwd,&acceptor_prob_fwd,
					left_part,right_part,genome,genomealt,chromosome_iit,breakpoint_start,breakpoint_end);
    bestproduct_rev = 0.0;

  } else if (*try_cdna_direction == -1) {
    *found_cdna_direction = -1;
    bestproduct_rev = find_exonexon_rev(&exonexonpos_rev,&donor1_rev,&donor2_rev,&acceptor2_rev,&acceptor1_rev,
					&comp_rev,&donor_watsonp_rev,&acceptor_watsonp_rev,&donor_prob_rev,&acceptor_prob_rev,
					left_part,right_part,genome,genomealt,chromosome_iit,breakpoint_start,breakpoint_end);
    bestproduct_fwd = 0.0;

  } else {
    bestproduct_fwd = find_exonexon_fwd(&exonexonpos_fwd,&donor1_fwd,&donor2_fwd,&acceptor2_fwd,&acceptor1_fwd,
					&comp_fwd,&donor_watsonp_fwd,&acceptor_watsonp_fwd,&donor_prob_fwd,&acceptor_prob_fwd,
					left_part,right_part,genome,genomealt,chromosome_iit,breakpoint_start,breakpoint_end);
    bestproduct_rev = find_exonexon_rev(&exonexonpos_rev,&donor1_rev,&donor2_rev,&acceptor2_rev,&acceptor1_rev,
					&comp_rev,&donor_watsonp_rev,&acceptor_watsonp_rev,&donor_prob_rev,&acceptor_prob_rev,
					left_part,right_part,genome,genomealt,chromosome_iit,breakpoint_start,breakpoint_end);
  }

  if (bestproduct_fwd == 0.0 && bestproduct_rev == 0.0) {
    *found_cdna_direction = 0;
    *donor1 = 'N';
    *donor2 = 'N';
    *acceptor2 = 'N';
    *acceptor1 = 'N';
    *comp = NONINTRON_COMP;
    *donor_watsonp = true;
    *acceptor_watsonp = true;
    *donor_prob = 0.0;
    *acceptor_prob = 0.0;
    return -1;

  } else if (bestproduct_fwd >= bestproduct_rev) {
    *found_cdna_direction = +1;
    *donor1 = donor1_fwd;
    *donor2 = donor2_fwd;
    *acceptor2 = acceptor2_fwd;
    *acceptor1 = acceptor1_fwd;
    *comp = comp_fwd;
    *donor_watsonp = donor_watsonp_fwd;
    *acceptor_watsonp = acceptor_watsonp_fwd;
    *donor_prob = donor_prob_fwd;
    *acceptor_prob = acceptor_prob_fwd;
    return exonexonpos_fwd;

  } else {
    *found_cdna_direction = -1;
    *donor1 = donor1_rev;
    *donor2 = donor2_rev;
    *acceptor2 = acceptor2_rev;
    *acceptor1 = acceptor1_rev;
    *comp = comp_rev;
    *donor_watsonp = donor_watsonp_rev;
    *acceptor_watsonp = acceptor_watsonp_rev;
    *donor_prob = donor_prob_rev;
    *acceptor_prob = acceptor_prob_rev;
    return exonexonpos_rev;
  }
}


