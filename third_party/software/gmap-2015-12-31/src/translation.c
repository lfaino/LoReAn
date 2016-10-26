static char rcsid[] = "$Id: translation.c 155282 2014-12-12 19:42:54Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "translation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For toupper */
#include "mem.h"
#include "comp.h"
#include "pairdef.h"
#include "complement.h"
#include "list.h"


#define IGNORE_MARGIN 6
#define HORIZON 99
#define MIN_NPAIRS 30


/* Finding coding sequence bounds */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Translation via reference */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Mark and assign cdna */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Translate cDNA for PMAP */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif



typedef enum {FRAME0, FRAME1, FRAME2, NOFRAME} Frame_T;
static char complCode[128] = COMPLEMENT_UC;
static char uppercaseCode[128] = UPPERCASE_U2T;

#define T Translation_T
struct T {
  int querypos;
  char aa;
  Frame_T frame;
};


static struct T *
Translation_array_new (struct Pair_T *pairs, int translationlen) {
  struct T *new;
  int i;

  new = (struct T *) CALLOC(translationlen,sizeof(struct T));
  for (i = 0; i < translationlen; i++) {
    new[i].querypos = pairs[i].querypos;
    new[i].aa = ' ';
    new[i].frame = NOFRAME;
  }

  return new;
}

static void
Translation_dump (struct Pair_T *pairs, struct T *translation, int translationlen) {
  int i;

  for (i = 0; i < translationlen; i++) {
    if (pairs[i].aaphase_g == 0 || pairs[i].aaphase_e == 0) {
      printf("=> %c %c ",pairs[i].aa_g,pairs[i].aa_e);
    } else {
      printf("       ");
    }
    printf("%d: querypos %d %d ",i,pairs[i].querypos,pairs[i].aapos);
    switch (translation[i].frame) {
    case NOFRAME: printf("%c %c %c",' ',' ',' '); break;
    case FRAME0: printf("%c %c %c",translation[i].aa,' ',' '); break;
    case FRAME1: printf("%c %c %c",' ',translation[i].aa,' '); break;
    case FRAME2: printf("%c %c %c",' ',' ',translation[i].aa); break;
    }
    printf("\n");
  }
  return;
}


/************************************************************************/


char
Translation_get_codon (char a, char b, char c) {
  switch (b) {
  case 'T':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'F';
      case 'A': case 'G': case 'R': return 'L';
      default: 	return 'X';
      }
    case 'C': return 'L';
    case 'A': 
      switch (c) {
      case 'G': return 'M';
      case 'T': case 'A': case 'C': case 'H': return 'I';
      default: return 'X';
      }
    case 'G': return 'V';
    default: return 'X';
  }
  case 'C':
    switch (a) {
    case 'T': return 'S';
    case 'C': return 'P';
    case 'A': return 'T';
    case 'G': return 'A';
    default: return 'X';
    }
  case 'A':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'Y';
      case 'A': case 'G': case 'R': return '*';
      default: return 'X';
      }
    case 'C':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'H';
      case 'A': case 'G': case 'R': return 'Q';
      default: return 'X';
      }
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'N';
      case 'A': case 'G': case 'R': return 'K';
      default: return 'X';
      }
    case 'G':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'D';
      case 'A': case 'G': case 'R': return 'E';
      default: return 'X';
      }
    default: return 'X';
    }
  case 'G':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'C';
      case 'A': return '*';
      case 'G': return 'W';
      default: return 'X';
      }
    case 'C': return 'R';
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'S';
      case 'A': case 'G': case 'R': return 'R';
      default: return 'X';
      }
    case 'G': return 'G';
    default: return 'X';
    }
  default: return 'X';
  }
  return 'X';
}


#ifndef PMAP
static void
find_bounds_forward (Frame_T *translation_frame, int *translation_starti, 
		     int *translation_endi, int *translation_length,
		     bool *endstopp, struct T *translation, 
		     int translationlen, bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = 0, start1 = 0, start2 = 0;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = 0; i < translationlen; i++) {
    debug(printf("%d %c: %d %d %d\n",i,translation[i].aa,orf0,orf1,orf2));
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    beststart0 = start0;
	    bestend0 = i;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  debug(printf("Incrementing orf0\n"));
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    beststart1 = start1;
	    bestend1 = i;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  debug(printf("Incrementing orf1\n"));
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    beststart2 = start2;
	    bestend2 = i;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  debug(printf("Incrementing orf2\n"));
	  orf2++;
	}
      } else {
	fprintf(stderr,"No frame at %d\n",i);
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    beststart0 = start0;
    bestend0 = translationlen-1;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    beststart1 = start1;
    bestend1 = translationlen-1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    beststart2 = start2;
    bestend2 = translationlen-1;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_starti = beststart2;
    *translation_endi = bestend2;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_starti = beststart1;
    *translation_endi = bestend1;
  } else if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_starti = beststart0;
    *translation_endi = bestend0;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}
#endif

#ifndef PMAP
static void
find_bounds_backward (Frame_T *translation_frame, int *translation_starti,
		      int *translation_endi, int *translation_length,
		      bool *endstopp, struct T *translation, int translationlen, bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = translationlen-1, start1 = translationlen-1, start2 = translationlen-1;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = translationlen-1; i >= 0; --i) {
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    bestend0 = i;
	    beststart0 = start0;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    bestend1 = i;
	    beststart1 = start1;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    bestend2 = i;
	    beststart2 = start2;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  orf2++;
	}
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    bestend0 = 0;
    beststart0 = start0;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    bestend1 = 0;
    beststart1 = start1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    bestend2 = 0;
    beststart2 = start2;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_starti = beststart0;
    *translation_endi = bestend0;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_starti = beststart1;
    *translation_endi = bestend1;
  } else if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_starti = beststart2;
    *translation_endi = bestend2;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}
#endif


#ifndef PMAP
static void
find_bounds_forward_fromstart (Frame_T *translation_frame, int *translation_starti, 
			       int *translation_endi, int *translation_length,
			       bool *endstopp, struct T *translation, int translationlen,
			       int cds_startpos) {
  Frame_T frame_fromstart;
  int phase_fromstart;
  int beststart0, bestend0;
  int bestorf0 = 0, orf0 = 0;
  int start0 = 0;
  bool endstop0p = false;
  char codon;
  int i;

  phase_fromstart = (cds_startpos - 1) % 3;
  if (translation[0].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[0].frame;
  } else if (translation[1].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[1].frame;
  } else if (translation[2].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[2].frame;
  } else {
    fprintf(stderr,"Something wrong with Translation_T object\n");
    abort();
  }

  for (i = 0; i < translationlen && endstop0p == false; i++) {
    if (translation[i].querypos >= cds_startpos - 1 &&
	translation[i].frame == frame_fromstart &&
	(codon = translation[i].aa) != ' ') {
      debug(printf("%d %c\n",i,translation[i].aa));
      if (orf0 == 0) {
	start0 = i;
      }
      if (codon == '*') {
	orf0++;
	if (orf0 > bestorf0) {
	  debug(printf("Frame 0: Best orf is %d\n",orf0));
	  bestorf0 = orf0;
	  beststart0 = start0;
	  bestend0 = i;
	  endstop0p = true;
	}

      } else {
	debug(printf("Incrementing orf0\n"));
	orf0++;
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    beststart0 = start0;
    bestend0 = translationlen-1;
    endstop0p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;

  debug(printf("Assigning frame %d\n",frame_fromstart));
  *translation_frame = frame_fromstart;
  *translation_starti = beststart0;
  *translation_endi = bestend0;

  return;
}
#endif


#ifndef PMAP
static void
find_bounds_backward_fromstart (Frame_T *translation_frame, int *translation_starti,
				int *translation_endi, int *translation_length,
				bool *endstopp, struct T *translation, int translationlen,
				int cds_startpos, int querylength) {
  Frame_T frame_fromstart;
  int phase_fromstart;
  int beststart0, bestend0;
  int bestorf0 = 0, orf0 = 0;
  int start0 = translationlen-1;
  bool endstop0p = false;
  char codon;
  int i;

  phase_fromstart = (translationlen - cds_startpos) % 3;
  if (translation[translationlen-1].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[translationlen-1].frame;
  } else if (translation[translationlen-2].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[translationlen-2].frame;
  } else if (translation[translationlen-3].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[translationlen-3].frame;
  } else {
    fprintf(stderr,"Something wrong with Translation_T object\n");
    abort();
  }
    
  for (i = translationlen-1; i >= 0 && endstop0p == false; --i) {
    if (translation[i].querypos >= cds_startpos - 1 &&
	translation[i].frame == frame_fromstart &&
	(codon = translation[i].aa) != ' ') {
      debug(printf("%d %c\n",i,translation[i].aa));
      if (orf0 == 0) {
	start0 = i;
      }
      if (codon == '*') {
	orf0++;
	if (orf0 > bestorf0) {
	  debug(printf("Frame 0: Best orf is %d\n",orf0));
	  bestorf0 = orf0;
	  bestend0 = i;
	  beststart0 = start0;
	  endstop0p = true;
	}
      } else {
	debug(printf("Incrementing orf0\n"));
	orf0++;
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    bestend0 = 0;
    beststart0 = start0;
    endstop0p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;

  debug(printf("Assigning frame %d\n",frame_fromstart));
  *translation_frame = frame_fromstart;
  *translation_starti = beststart0;
  *translation_endi = bestend0;

  return;
}
#endif



#ifdef PMAP
static void
translate_pairs_cdna (int *translation_starti, int *translation_endi, int *translation_length,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr) {
  struct Pair_T *ptr, *pair;
  int i;

  /* Go backward so we put aa at beginning of codon */
  /* printf("lastquerypos is %d\n",pairs[npairs-1].querypos); */
  i = npairs-1;
  while (i >= 0 && pairs[i].querypos % 3 != 2) {
    i--;
  }
  *translation_endi = i;
  if (*translation_endi < 0) {
    *translation_endi = 0;
  }

  debug3(printf("Entering translate_pairs_cdna with npairs=%d, translation_endi = %d\n",
		npairs,*translation_endi));

  *translation_starti = *translation_endi;
  *translation_length = 0;

  ptr = &(pairs[(*translation_endi)+1]);
  for (i = *translation_endi; i >= 0; --i) {
    pair = --ptr;
    pair->aapos = pair->querypos/3 + 1;

    if (pair->cdna == ' ') {
      /* pair->aa_e = ' '; */
    } else {
      if ((pair->aaphase_e = pair->querypos % 3) == 0) {
	pair->aa_e = queryaaseq_ptr[pair->querypos/3];
	*translation_starti = i;
	(*translation_length) += 1;
      }
    }
  }

  return;
}  
#endif


static struct Translation_T *
translate_pairs_forward (struct Pair_T *pairs, int npairs, bool revcompp) {
  struct T *translation;
  struct Pair_T *ptr, *pair;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';

  translation = Translation_array_new(pairs,npairs);

  /* Go backward so we put aa at beginning of codon */
  ptr = &(pairs[npairs]);
  for (i = npairs-1; i >= 0; --i) {
    pair = --ptr;
    if (pair->gapp == true) {
      /* translation[i].aa = ' '; */
    } else if (pair->genome == ' ') {
      /* translation[i].aa = ' '; */
    } else {
      nt2 = nt1;
      nt1 = nt0;
      nt0 = revcompp ? complCode[(int) pair->genome] : uppercaseCode[(int) pair->genome];

      codon = Translation_get_codon(nt0,nt1,nt2);
      if (gpos < 2 && codon == 'X') {
	/* translation[i].aa = ' '; */
      } else {
	translation[i].aa = codon;
	switch (gpos % 3) {
	case 0: translation[i].frame = FRAME0; break;
	case 1: translation[i].frame = FRAME1; break;
	case 2: translation[i].frame = FRAME2; break;
	}
      }
      gpos++;
    }
  }

  return translation;
}  

static struct Translation_T *
translate_pairs_backward (struct Pair_T *pairs, int npairs, bool revcompp) {
  struct T *translation;
  struct Pair_T *ptr, *pair;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';

  translation = Translation_array_new(pairs,npairs);

  /* Go forwards so we put aa at beginning of codon */
  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    pair = ptr++;
    if (pair->gapp == true) {
      /* translation[i].aa = ' '; */
    } else if (pair->genome == ' ') {
      /* translation[i].aa = ' '; */
    } else {
      nt2 = nt1;
      nt1 = nt0;
      nt0 = revcompp ? complCode[(int) pair->genome] : uppercaseCode[(int) pair->genome];

      codon = Translation_get_codon(nt0,nt1,nt2);
      if (gpos < 2 && codon == 'X') {
	/* translation[i].aa = ' '; */
      } else {
	translation[i].aa = codon;
	switch (gpos % 3) {
	case 0: translation[i].frame = FRAME0; break;
	case 1: translation[i].frame = FRAME1; break;
	case 2: translation[i].frame = FRAME2; break;
	}
      }
      gpos++;
    }
  }

  return translation;
}  


/************************************************************************/


#ifndef PMAP
static int
count_cdna_forward_strict (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;

  j = starti;
  while (j < npairs) {
    if (ncdna >= 3 && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j++;
  }

  *nexti = j;
  return ncdna;
}
#endif


static int
count_cdna_forward (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;

  j = starti;
  while (j <= endi) {
    if (j > starti && pairs[j].aaphase_g == 0 && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j++;
  }

  *nexti = j;
  return ncdna;
}


static int
count_cdna_forward_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;
  
  j = starti;
  while (j <= endi && ncdna <= HORIZON) {
    if (j > starti && pairs[j].aaphase_g == 0 && pairs[j].cdna != ' ' &&
	(ncdna % 3) == 0) {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j++;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}

#ifndef PMAP
static int
count_cdna_backward_strict (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;

  j = starti;
  while (j >= 0) {
    if (ncdna >= 3 && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j--;
  }

  *nexti = j;
  return ncdna;
}
#endif

static int
count_cdna_backward (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;

  j = starti;
  while (j >= endi) {
    if (j < starti && pairs[j].aaphase_g == 0 && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j--;
  }

  *nexti = j;
  return ncdna;
}


static int
count_cdna_backward_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ncdna = 0, j;
  
  j = starti;
  while (j >= endi && ncdna <= HORIZON) {
    if (j < starti && pairs[j].aaphase_g == 0 && pairs[j].cdna != ' ' &&
	(ncdna % 3) == 0) {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].cdna != ' ') {
      ncdna++;
    }
    j--;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}


static int
count_genomic_strict (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ngenomic = 0, j;

  j = starti;
  while (j < npairs) {
    /* Need to check gapp in addition to genome, because genome
       characters are provided at ends of introns */
    if (ngenomic >= 3 && pairs[j].gapp == false && pairs[j].genome != ' ') {
      *nexti = j;
      return ngenomic;
    } else if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      ngenomic++;
    }
    j++;
  }

  *nexti = j;
  return ngenomic;
}


static int
count_genomic (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ngenomic = 0, j;

  j = starti;
  while (j <= endi) {
    if (j > starti && pairs[j].aaphase_e == 0 && pairs[j].genome != ' ') {
      *nexti = j;
      return ngenomic;
    } else if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      ngenomic++;
    }
    j++;
  }

  *nexti = j;
  return ngenomic;
}


static int
count_genomic_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int endi) {
  int ngenomic = 0, j;
  
  j = starti;
  while (j <= endi && ngenomic <= HORIZON) {
    if (j > starti && pairs[j].aaphase_e == 0 && pairs[j].genome != ' ' &&
	(ngenomic % 3) == 0) {
      *nexti = j;
      return ngenomic;
    } else if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      ngenomic++;
    }
    j++;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}

/************************************************************************/

static int
get_codon_forward (int *nexti, struct Pair_T *pairs, int npairs, int starti, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int j2 = -1, j1 = -1, j0 = -1;
  int ncdna = 0, j;

  j = starti;
  while (j < npairs && ncdna < 3) {
    if (pairs[j].cdna != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[(int) pairs[j].cdna] : uppercaseCode[(int) pairs[j].cdna];
      j0 = j1;
      j1 = j2;
      j2 = j;
      ncdna++;
    }
    j++;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j <= npairs - 1 && pairs[j].cdna == ' ') {
    j++;
  }

  *nexti = j;
  
  if (j0 >= 0) {
    pairs[j0].aaphase_e = 0;
    pairs[j1].aaphase_e = 1;
    pairs[j2].aaphase_e = 2;
  }

  return Translation_get_codon(nt0,nt1,nt2);
}

static int
get_codon_backward (int *nexti, struct Pair_T *pairs, int npairs, int starti, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int j2 = -1, j1 = -1, j0 = -1;
  int ncdna = 0, j;

  j = starti;
  while (j >= 0 && ncdna < 3) {
    if (pairs[j].cdna != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[(int) pairs[j].cdna] : uppercaseCode[(int) pairs[j].cdna];
      j0 = j1;
      j1 = j2;
      j2 = j;
      ncdna++;
    }
    --j;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j >= 0 && (pairs[j].cdna == ' ')) {
    --j;
  }

  *nexti = j;

  if (j0 >= 0) {
    pairs[j0].aaphase_e = 0;
    pairs[j1].aaphase_e = 1;
    pairs[j2].aaphase_e = 2;
  }

  return Translation_get_codon(nt0,nt1,nt2);
}


#ifdef PMAP
static int
get_codon_genomic (int *nexti, struct Pair_T *pairs, int npairs, int starti) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int j0, j1, j2;
  int ngenomic = 0, j;

  j = starti;
  while (ngenomic < 3) {
    if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = uppercaseCode[(int) pairs[j].genome];
      j0 = j1;
      j1 = j2;
      j2 = j;
      ngenomic++;
    }
    j++;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j <= npairs-1 && (pairs[j].gapp == true || pairs[j].genome == ' ')) {
    j++;
  }

  *nexti = j;

  pairs[j0].aaphase_g = 0;
  pairs[j1].aaphase_g = 1;
  pairs[j2].aaphase_g = 2;

  return Translation_get_codon(nt0,nt1,nt2);
}
#endif



static int
assign_cdna_forward (int ncdna, struct Pair_T *pairs, int npairs, bool revcompp, int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0;
  int codon = ' ';

  i = starti;
  while (i < npairs && pairs[i].cdna == ' ') {
    i++;
  }
  while (j < ncdna) {
    pair = &(pairs[i]);
    codon = pair->aa_e = get_codon_forward(&nexti,pairs,npairs,i,revcompp);
    debug2(Pair_dump_one(pair,true));
    debug2(printf(" marked with amino acid %c\n",pair->aa_e));
    i = nexti;
    j += 3;
  }
  return codon;
}

static void
terminate_cdna_forward (struct Pair_T *pairs, int npairs, bool revcompp, int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  while (i < npairs && pairs[i].cdna == ' ') {
    i++;
  }
  while (i <= npairs-3 && lastcodon != '*') {
    pair = &(pairs[i]);
    lastcodon = pair->aa_e = get_codon_forward(&nexti,pairs,npairs,i,revcompp);
    debug2(Pair_dump_one(pair,true));
    debug2(printf(" marked with amino acid %c\n",pair->aa_e));
    i = nexti;
  }
  return;
}

static int
assign_cdna_backward (int ncdna, struct Pair_T *pairs, int npairs, bool revcompp, int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0;
  int codon = ' ';

  i = starti;
  while (i >= 0 && pairs[i].cdna == ' ') {
    --i;
  }
  while (j < ncdna) {
    pair = &(pairs[i]);
    codon = pair->aa_e = get_codon_backward(&nexti,pairs,npairs,i,revcompp);
    debug2(Pair_dump_one(pair,true));
    debug2(printf(" marked with amino acid %c\n",pair->aa_e));
    i = nexti;
    j += 3;
  }
  return codon;
}

static void
terminate_cdna_backward (struct Pair_T *pairs, int npairs, bool revcompp, int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  while (i >= 0 && pairs[i].cdna == ' ') {
    --i;
  }

  /* i > 1 is equivalent to i >= 2 */
  while (i > 1 && lastcodon != '*') {
    pair = &(pairs[i]);
    lastcodon = pair->aa_e = get_codon_backward(&nexti,pairs,npairs,i,revcompp);
    debug2(Pair_dump_one(pair,true));
    debug2(printf(" marked with amino acid %c\n",pair->aa_e));
    i = nexti;
  }
  return;
}

#ifdef PMAP
static int
assign_genomic (int ngenomic, struct Pair_T *pairs, int npairs, int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0, codon;

  i = starti;
  while (j < ngenomic) {
    pair = &(pairs[i]);
    codon = pair->aa_g = get_codon_genomic(&nexti,pairs,npairs,i);
    i = nexti;
    j += 3;
  }
  return codon;
}

static void
terminate_genomic (struct Pair_T *pairs, int npairs, int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  while (i <= npairs - 3 && lastcodon != '*') {
    pair = &(pairs[i]);
    lastcodon = pair->aa_g = get_codon_genomic(&nexti,pairs,npairs,i);
    i = nexti;
  }
  return;
}
#endif



#ifndef PMAP
static void
mark_cdna_forward_strict (struct Pair_T *pairs, int npairs, bool revcompp, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, ncdna, codon = ' ';

  debug2(printf("mark_cdna_forward_strict called with pairs #%d..%d\n",starti,endi));

  i = starti;

  /* Advance to same start as genomic */
  pair = &(pairs[i]);
  while (i < endi && pair->aaphase_g != 0) {
    i++;
    pair = &(pairs[i]);
  }

  while (i < npairs && codon != '*') {
    pair = &(pairs[i]);
    ncdna = count_cdna_forward_strict(&nexti,pairs,npairs,i,endi);
    if (ncdna == 3) {
      codon = assign_cdna_forward(ncdna,pairs,npairs,revcompp,i);
    }
    i = nexti;
  }

  if (codon != '*') {
    terminate_cdna_forward(pairs,npairs,revcompp,i);
  }

  return;
}  
#endif

static void
mark_cdna_forward (struct Pair_T *pairs, int npairs, bool revcompp, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ncdna, ncdna_alt;

  debug2(printf("mark_cdna_forward called with pairs #%d..%d\n",starti,endi));

  i = starti;
  while (i < endi) {
    pair = &(pairs[i]);
    if (pair->aaphase_g != 0) {
      i++;
    } else {
      ncdna = count_cdna_forward(&nexti,pairs,npairs,i,endi);
      if (ncdna == 3) {
	assign_cdna_forward(ncdna,pairs,npairs,revcompp,i);
      } else if ((ncdna % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ncdna));
	assign_cdna_forward(ncdna,pairs,npairs,revcompp,i);
      } else if (i + 2 > endi) {
	debug2(printf("At %d, saw last codon: %d+2 > %d\n",pair->aapos,i,endi));
	assign_cdna_forward(ncdna,pairs,npairs,revcompp,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ncdna));	
	ncdna_alt = count_cdna_forward_mod3(&nexti_alt,pairs,npairs,i,endi);

	debug2(printf("  Alternate search yields %d\n",ncdna_alt));
	if ((ncdna_alt % 3) == 0) {
	  debug2(printf("  Using alternate search\n"));
	  assign_cdna_forward(ncdna_alt,pairs,npairs,revcompp,i);
	  nexti = nexti_alt;

	} else if (ncdna < 3) {
	  debug2(printf("  Skipping\n"));

	} else {
	  debug2(printf("  Using original count - 3 = %d\n",ncdna-3));
	  assign_cdna_forward(ncdna-3,pairs,npairs,revcompp,i);
	}
      }
      i = nexti;
    }
  }

  debug2(printf("Calling terminate_cdna_forward\n"));
  terminate_cdna_forward(pairs,npairs,revcompp,i);

  return;
}  

#ifndef PMAP
static void
mark_cdna_backward_strict (struct Pair_T *pairs, int npairs, bool revcompp, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, ncdna, codon = ' ';

  debug2(printf("mark_cdna_backward_strict called with pairs #%d..%d\n",starti,endi));

  i = starti;

  /* Advance to same start as genomic */
  pair = &(pairs[i]);
  while (i > endi && pair->aaphase_g != 0) {
    i--;
    pair = &(pairs[i]);
  }

  while (i >= 0 && codon != '*') {
    pair = &(pairs[i]);
    ncdna = count_cdna_backward_strict(&nexti,pairs,npairs,i,endi);
    if (ncdna == 3) {
      codon = assign_cdna_backward(ncdna,pairs,npairs,revcompp,i);
    }
    i = nexti;
  }

  if (codon != '*') {
    terminate_cdna_backward(pairs,npairs,revcompp,i);
  }

  return;
}  
#endif

static void
mark_cdna_backward (struct Pair_T *pairs, int npairs, bool revcompp, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ncdna, ncdna_alt;

  debug2(printf("mark_cdna_backward called with pairs #%d..%d\n",starti,endi));

  i = starti;
  while (i > endi) {
    pair = &(pairs[i]);
    if (pair->aaphase_g != 0) {
      i--;
    } else {
      ncdna = count_cdna_backward(&nexti,pairs,npairs,i,endi);
      if (ncdna == 3) {
	assign_cdna_backward(ncdna,pairs,npairs,revcompp,i);
      } else if ((ncdna % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ncdna));
	assign_cdna_backward(ncdna,pairs,npairs,revcompp,i);
      } else if (i - 2 < endi) {
	debug2(printf("At %d, saw last codon: %d-2 < %d\n",pair->aapos,i,endi));
	assign_cdna_backward(ncdna,pairs,npairs,revcompp,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ncdna));	
	ncdna_alt = count_cdna_backward_mod3(&nexti_alt,pairs,npairs,i,endi);

	debug2(printf("  Alternate search yields %d\n",ncdna_alt));
	if ((ncdna_alt % 3) == 0) {
	  debug2(printf("  Using alternate search\n"));
	  assign_cdna_backward(ncdna_alt,pairs,npairs,revcompp,i);
	  nexti = nexti_alt;

	} else if (ncdna < 3) {
	  debug2(printf("Skipping\n"));

	} else {
	  debug2(printf("  Using original count - 3 = %d\n",ncdna-3));
	  assign_cdna_backward(ncdna-3,pairs,npairs,revcompp,i);
	}
      }
      i = nexti;
    }
  }

  debug2(printf("Calling terminate_cdna_backward\n"));
  terminate_cdna_backward(pairs,npairs,revcompp,i);

  return;
}  


#ifdef PMAP
static void
mark_genomic_strict (struct Pair_T *pairs, int npairs, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, ngenomic, codon = ' ';

  debug2(printf("mark_genomic_strict called with pairs #%d %d\n",starti,endi));

  i = starti;

  /* Advance to same start as cDNA */
  pair = &(pairs[i]);
  while (i < endi && pair->aaphase_e != 0) {
    i++;
    pair = &(pairs[i]);
  }

  while (i < npairs && codon != '*') {
    pair = &(pairs[i]);
    ngenomic = count_genomic_strict(&nexti,pairs,npairs,i,endi);
    if (ngenomic == 3) {
      codon = assign_genomic(ngenomic,pairs,npairs,i);
    }
    i = nexti;
  }

  if (codon != '*') {
    terminate_genomic(pairs,npairs,i);
  }

  return;
}  

static void
mark_genomic (struct Pair_T *pairs, int npairs, int starti, int endi) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ngenomic, ngenomic_alt;

  debug2(printf("mark_genomic called with pairs #%d %d\n",starti,endi));

  i = starti;
  while (i < endi) {
    pair = &(pairs[i]);
    if (pair->aaphase_e != 0) {
      i++;
    } else {
      ngenomic = count_genomic(&nexti,pairs,npairs,i,endi);
      if (ngenomic == 3) {
	assign_genomic(ngenomic,pairs,npairs,i);
      } else if ((ngenomic % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ngenomic));
	assign_genomic(ngenomic,pairs,npairs,i);
      } else if (i + 2 > endi) {
	debug2(printf("At %d, saw last codon: %d+2 > %d\n",pair->aapos,i,endi));
	assign_genomic(ngenomic,pairs,npairs,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ngenomic));
	ngenomic_alt = count_genomic_mod3(&nexti_alt,pairs,npairs,i,endi);

	debug2(printf("  Alternate search yields %d\n",ngenomic_alt));
	if ((ngenomic_alt % 3) == 0) {
	  debug2(printf("  Using alternate search\n"));
	  assign_genomic(ngenomic_alt,pairs,npairs,i);
	  nexti = nexti_alt;

	} else if (ngenomic < 3) {
	  debug2(printf("Skipping\n"));

	} else {
	  debug2(printf("  Using original count - 3 = %d\n",ngenomic-3));
	  assign_genomic(ngenomic-3,pairs,npairs,i);
	}
      }
      i = nexti;
    }
  }

  terminate_genomic(pairs,npairs,i);

  return;
}  
#endif



#ifdef PMAP
void
Translation_via_cdna (int *translation_leftpos, int *translation_rightpos, int *translation_length,
		      int *relaastart, int *relaaend,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr, bool strictp) {
  int translation_starti, translation_endi, i;

  /* Don't check for MIN_NPAIRS */

  for (i = 0; i < npairs; i++) {
    pairs[i].refquerypos = pairs[i].querypos;
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  translate_pairs_cdna(&translation_starti,&translation_endi,&(*translation_length),
		       pairs,npairs,queryaaseq_ptr);

  *translation_leftpos = pairs[translation_starti].querypos;
  *translation_rightpos = pairs[translation_endi].querypos;

  debug(printf("Translation start = pair #%d (querypos %d)\n",translation_starti,*translation_leftpos));
  debug(printf("Translation end = pair #%d (querypos %d)\n",translation_endi,*translation_rightpos));

  /* Take care of genomic side */
    
  *relaastart = pairs[translation_starti].aapos;
  *relaaend = pairs[translation_endi].aapos;

  if (strictp == true) {
    mark_genomic_strict(pairs,npairs,translation_starti,translation_endi);
  } else {
    mark_genomic(pairs,npairs,translation_starti,translation_endi);
  }

  return;
}
#else
void
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 int *relaastart, int *relaaend,
			 struct Pair_T *pairs, int npairs, bool backwardp, bool revcompp, bool fulllengthp,
			 int cds_startpos, int querylength, bool strictp) {
  char lastaa;
  struct T *translation;
  bool endstopp;
  int i, aapos = 0;
  Frame_T translation_frame;
  int translation_starti = 0, translation_endi = 0, phase;
  int minpos, maxpos;

  if (npairs < MIN_NPAIRS) {
    *translation_leftpos = 0;
    *translation_rightpos = 0;
    *translation_length = 0;
    *relaastart = 0;
    *relaaend = 0;
    return;
  }

  if (backwardp == false) {
    translation = translate_pairs_forward(pairs,npairs,revcompp);
    if (cds_startpos > 0) {
      find_bounds_forward_fromstart(&translation_frame,&translation_starti,
				    &translation_endi,&(*translation_length),&endstopp,
				    translation,npairs,cds_startpos);
    } else {
      find_bounds_forward(&translation_frame,&translation_starti,
			  &translation_endi,&(*translation_length),&endstopp,
			  translation,npairs,fulllengthp);
      if (fulllengthp == true && *translation_length == 0) {
	/* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
	find_bounds_forward(&translation_frame,&translation_starti,
			    &translation_endi,&(*translation_length),&endstopp,
			    translation,npairs,/*fulllengthp*/false);
      }
    }

  } else {
    translation = translate_pairs_backward(pairs,npairs,revcompp);
    if (cds_startpos > 0) {
      find_bounds_backward_fromstart(&translation_frame,&translation_starti,
				     &translation_endi,&(*translation_length),&endstopp,
				     translation,npairs,cds_startpos,querylength);
    } else {
      find_bounds_backward(&translation_frame,&translation_starti,
			   &translation_endi,&(*translation_length),&endstopp,
			   translation,npairs,fulllengthp);
      if (fulllengthp == true && *translation_length == 0) {
	/* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
	find_bounds_backward(&translation_frame,&translation_starti,
			     &translation_endi,&(*translation_length),&endstopp,
			     translation,npairs,/*fulllengthp*/false);
      }
    }
  }
  /* debug(printf("ref:\n")); */
  debug(Translation_dump(pairs,translation,npairs));

  for (i = 0; i < npairs; i++) {
    pairs[i].refquerypos = pairs[i].querypos;
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  if (translation_starti < 0 || translation_endi < 0) {
    *translation_leftpos = *translation_rightpos = -1;
    *relaastart = *relaaend = -1;
  } else {
    minpos = pairs[npairs-1].querypos;
    maxpos = pairs[0].querypos;
    if (backwardp == false) {
      /* forward */
      debug(printf("Translation is forward pairs %d..%d\n",translation_starti,translation_endi));
      for (i = translation_starti; i <= translation_endi; i++) {
	/* Avoid problems with genome position advancing prematurely */
	if (pairs[i].genome != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((pairs[i].aa_g = translation[i].aa) != ' ') {
	      if (pairs[i].querypos < minpos) {
		minpos = pairs[i].querypos;
	      }
	      if (pairs[i].querypos > maxpos) {
		maxpos = pairs[i].querypos;
	      }
	      lastaa = pairs[i].aa_g;
	      aapos++;
	      pairs[i].aaphase_g = 0;
	    }
	  } else if (translation[i].frame != 3) {
	    if ((phase = translation_frame - translation[i].frame) < 0) {
	      phase += 3;
	    }
	    pairs[i].aaphase_g = phase;
	  }
	}
	pairs[i].aapos = aapos;
      }
      *translation_leftpos = minpos;
      if ((*translation_rightpos = maxpos + 2) >= querylength) {
	*translation_rightpos = querylength - 1;
      }
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
#if 0
      if (i < npairs) {
	*translation_rightpos += 1;
	pairs[i++].aapos = aapos;
      }
      if (i < npairs) {
	*translation_rightpos += 1;
	pairs[i].aapos = aapos;
      }
#endif
    
      /* Fill in aapos to the end */
      for ( ; i < npairs; i++) {
	pairs[i].aapos = aapos;
      }

    } else {
      /* backward */
      debug(printf("Translation is backward pairs %d..%d\n",translation_starti,translation_endi));
      
      for (i = translation_starti; i >= translation_endi; --i) {
	/* Avoid problems with genome position advancing prematurely */
	if (pairs[i].genome != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((pairs[i].aa_g = translation[i].aa) != ' ') {
	      if (pairs[i].querypos < minpos) {
		minpos = pairs[i].querypos;
	      }
	      if (pairs[i].querypos > maxpos) {
		maxpos = pairs[i].querypos;
	      }
	      lastaa = pairs[i].aa_g;
	      aapos++;
	      pairs[i].aaphase_g = 0;
	    }
	  } else if (translation[i].frame != 3) {
	    if ((phase = translation_frame - translation[i].frame) < 0) {
	      phase += 3;
	    }
	    pairs[i].aaphase_g = phase;
	  }
	}
	pairs[i].aapos = aapos;
      }
      if ((*translation_leftpos = minpos - 2) < 0) {
	*translation_leftpos = 0;
      }
      *translation_rightpos = maxpos;
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
#if 0
      if (i >= 0) {
	*translation_leftpos -= 1;
	pairs[i--].aapos = aapos;
      }
      if (i >= 0) {
	*translation_leftpos -= 1;
	pairs[i].aapos = aapos;
      }
#endif

      /* Fill in aapos to the end */
      for ( ; i >= 0; --i) {
	pairs[i].aapos = aapos;
      }
    }
    
    /* Take care of cDNA side */
    
    *relaastart = pairs[translation_starti].aapos;
    *relaaend = pairs[translation_endi].aapos;

    debug(printf("Translation start = pair #%d (querypos %d, aa #%d)\n",
		 translation_starti,*translation_rightpos,*relaaend));
    debug(printf("Translation end = pair #%d (querypos %d, aa#%d)\n",
		 translation_endi,*translation_leftpos,*relaastart));
    debug(printf("Translation length = %d aa\n",*translation_length));
    
    if (strictp == true) {
      if (revcompp == false) {
	mark_cdna_forward_strict(pairs,npairs,revcompp,translation_starti,translation_endi);
      } else {
	mark_cdna_backward_strict(pairs,npairs,revcompp,translation_starti,translation_endi);
      }
    } else {
      if (revcompp == false) {
	mark_cdna_forward(pairs,npairs,revcompp,translation_starti,translation_endi);
      } else {
	mark_cdna_backward(pairs,npairs,revcompp,translation_starti,translation_endi);
      }
    }
  }

  FREE(translation);

  return;
}
#endif

/* Pairs are always ordered by ascending cDNA position and genomic
   position.  For Crick strand matches, the genomic position needs to
   be standardized. */
static void
bound_via_reference (int *start, int *end, struct Pair_T *pairs, int npairs, bool watsonp, 
		     struct Pair_T *refpairs, int nrefpairs, bool refwatsonp) {
  int i, j, aapos = 0;
  int refquerypos, genomepos, refgenomepos;

  debug(Pair_dump_array(refpairs,nrefpairs,false));

  *start = *end = -1;
  if (refwatsonp == true && watsonp == true) {
    debug1(printf("refwatsonp == true && watsonp == true\n"));
    i = 0;
    j = 0;
    while (i < nrefpairs && j < npairs) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	debug(printf("Not incrementing aapos %d\n",aapos));
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else {
	debug(printf("looking at refpairs %d, genomepos %d (%c)\n",
		     i,refgenomepos,refpairs[i].aa_e));
	if (refpairs[i].aa_e != ' ') {
	  if (*start < 0) {
	    *start = j;
	  }
	  *end = j;
	}
	pairs[j].aaphase_g = refpairs[i].aaphase_g;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      }
    }

    /*
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    */

  } else if (refwatsonp == true && watsonp == false) {
    debug1(printf("refwatsonp == true && watsonp == false\n"));
    i = 0;
    j = npairs-1;
    while (i < nrefpairs && j >= 0) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*end < 0) {
	    *end = j;
	    /*
	    aapos = refpairs[i].aapos;
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    */
	  }
	  *start = j;
	}
	pairs[j].aaphase_g = refpairs[i].aaphase_g;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      }
    }

  } else if (refwatsonp == false && watsonp == true) {
    debug1(printf("refwatsonp == false && watsonp == true\n"));
    i = nrefpairs-1;
    j = 0;
    while (i >= 0 && j < npairs) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*start < 0) {
	    *start = j;
	  }
	  *end = j;
	}
	pairs[j].aaphase_g = refpairs[i].aaphase_g;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      }
    }

    /*
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    */

  } else {
    debug1(printf("refwatsonp == false && watsonp == false\n"));
    i = nrefpairs-1;
    j = npairs-1;
    while (i >= 0 && j >= 0) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*end < 0) {
	    *end = j;
	    /*
	    aapos = refpairs[i].aapos;
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    */
	  }
	  *start = j;
	}
	pairs[j].aaphase_g = refpairs[i].aaphase_g;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      }
    }
  }    

  debug(printf("Bound by reference: start = %d, end = %d\n",*start,*end));
  debug(Pair_dump_array(pairs,npairs,false));

  return;
}


void
Translation_via_reference (int *relaastart, int *relaaend,
			   struct Pair_T *pairs, int npairs, bool watsonp, bool backwardp, bool revcompp,
			   struct Pair_T *refpairs, int nrefpairs, bool refwatsonp,
			   bool fixshiftp) {
  struct T *translation;
  int start, end, i;

  if (npairs < MIN_NPAIRS) {
    *relaastart = 0;
    *relaaend = 0;
    return;
  }

  for (i = 0; i < npairs; i++) {
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  debug2(printf("Translation_via_reference called with backwardp = %d\n",backwardp));
  bound_via_reference(&start,&end,pairs,npairs,watsonp,refpairs,nrefpairs,refwatsonp);

  if (start < 0 || end < 0) {
    *relaastart = *relaaend = -1;
  } else {
    *relaastart = pairs[start].aapos;
    *relaaend = pairs[end].aapos;

    if (backwardp == false) {
      translation = translate_pairs_forward(pairs,npairs,revcompp);
    } else {
      translation = translate_pairs_backward(pairs,npairs,revcompp);
    }

    for (i = 0; i < npairs; i++) {
      if (pairs[i].aaphase_g == 0) {
	pairs[i].aa_g = translation[i].aa;
      } else {
	pairs[i].aa_g = ' ';
      }
      pairs[i].aa_e = ' ';
    }

    debug2(printf("Bounding pairs #%d..%d to yield amino acid coordinates %d (%c) to %d (%c)\n",
		  start,end,*relaastart,pairs[start].aa_g,*relaaend,pairs[end].aa_g));

    if (backwardp == false) {
      mark_cdna_forward(pairs,npairs,revcompp,start,end);
    } else {
      /* Note that we need to flip start and end here */
      mark_cdna_backward(pairs,npairs,revcompp,end,start);
    }

    debug1(Translation_dump(pairs,translation,npairs));
    FREE(translation);

  }

  return;
}

#ifndef PMAP
static char
lookup_aa (int aapos, struct Pair_T *refpairs, int nrefpairs) {
  int i;

  for (i = 0; i < nrefpairs; i++) {
    if (refpairs[i].aapos == aapos && refpairs[i].aaphase_g == 0) {
      return refpairs[i].aa_g;
    }
  }
  return 'X';
}
#endif

static int
next_aapos_fwd (struct Pair_T *pairs, int i, int npairs, int aapos) {
  struct Pair_T *this;

  this = &(pairs[i]);
  while (i < npairs && this->aapos == aapos) {
    this++;
    i++;
  }

  /* Advance to next amino acid on query side, since aapos assures
     us we have reached it on genome side */
  while (i < npairs && this->aa_e == ' ') {
    this++;
    i++;
  }

  return i;
}

static int
next_aapos_rev (struct Pair_T *pairs, int i, int npairs, int aapos) {
  struct Pair_T *this;

  this = &(pairs[i]);
  while (i >= 0 && this->aapos == aapos) {
    --this;
    --i;
  }

  /* Advance to next amino acid on query side, since aapos assures
     us we have reached it on genome side */
  while (i >= 0 && this->aa_e == ' ') {
    --this;
    --i;
  }

  return i;
}


#define MAXMUT 100

static void
fill_aa_fwd (int *strlen_g, int *strlen_c, int *netchars, char *aa_genomicseg, char *aa_queryseq,
	     char *nt_genomicseg, char *nt_queryseq, struct Pair_T *start, struct Pair_T *end) {
  int k1 = 0, k2 = 0, i1 = 0, i2 = 0;
  struct Pair_T *this;

  *netchars = 0;
  for (this = start; this <= end && i1 < MAXMUT && k1 < MAXMUT; this++) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	nt_genomicseg[i1++] = uppercaseCode[(int) this->genome];
      } else {
	*netchars += 1;
      }
      if (this->aa_g != ' ') {
	aa_genomicseg[k1++] = this->aa_g;
      }
    }
  }

  for (this = start; this <= end && i2 < MAXMUT && k2 < MAXMUT; this++) {
    if (this->gapp == false) {
      if (this->cdna != ' ') {
	nt_queryseq[i2++] = uppercaseCode[(int) this->cdna];
      } else {
	*netchars -= 1;
      }
      if (this->aa_e != ' ') {
	aa_queryseq[k2++] = this->aa_e;
      }
    }
  }

  if (i1 >= MAXMUT || k1 >= MAXMUT || i2 >= MAXMUT || k2 >= MAXMUT) {
    i1 = i2 = k1 = k2 = 0;
  }
  nt_genomicseg[i1] = '\0';
  aa_genomicseg[k1] = '\0';
  nt_queryseq[i2] = '\0';
  aa_queryseq[k2] = '\0';

  *strlen_g = k1;
  *strlen_c = k2;

  return;
}

static void
fill_aa_rev (int *strlen_g, int *strlen_c, int *netchars, char *aa_genomicseg, char *aa_queryseq,
	     char *nt_genomicseg, char *nt_queryseq, struct Pair_T *start, struct Pair_T *end) {
  int k1 = 0, k2 = 0, i1 = 0, i2 = 0;
  struct Pair_T *this;

  *netchars = 0;
  for (this = end; this >= start && i1 < MAXMUT && k1 < MAXMUT; --this) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	nt_genomicseg[i1++] = uppercaseCode[(int) this->genome];
      } else {
	*netchars += 1;
      }
      if (this->aa_g != ' ') {
	aa_genomicseg[k1++] = this->aa_g;
      }
    }
  }

  for (this = end; this >= start && i2 < MAXMUT && k2 < MAXMUT; --this) {
    if (this->gapp == false) {
      if (this->cdna != ' ') {
	nt_queryseq[i2++] = uppercaseCode[(int) this->cdna];
      } else {
	*netchars -= 1;
      }
      if (this->aa_e != ' ') {
	aa_queryseq[k2++] = this->aa_e;
      }
    }
  }

  if (i1 >= MAXMUT || k1 >= MAXMUT || i2 >= MAXMUT || k2 >= MAXMUT) {
    i1 = i2 = k1 = k2 = 0;
  }
  nt_genomicseg[i1] = '\0';
  aa_genomicseg[k1] = '\0';
  nt_queryseq[i2] = '\0';
  aa_queryseq[k2] = '\0';

  *strlen_g = k1;
  *strlen_c = k2;

  return;
}


static void
print_mutation (Filestring_T fp, bool *printp, int aapos, int strlen_g, int strlen_c, int refquerypos,
		char *aa_genomicseg, char *aa_queryseq, char *nt_genomicseg, char *nt_queryseq) {
  bool print_refquerypos_p = true;
#if 0
  bool print_nt_p = true;
#endif

  if (strlen_g > strlen_c) {
    if (*printp == true) FPRINTF(fp,", "); else *printp = true;
    if (aa_genomicseg[0] == aa_queryseq[0]) {
      FPRINTF(fp,"del%s%d%s ",&(aa_genomicseg[1]),aapos+1,&(aa_queryseq[1]));
      refquerypos += 3;
    } else {
      FPRINTF(fp,"del%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
    }
  } else if (strlen_g < strlen_c) {
    if (*printp == true) FPRINTF(fp,", "); else *printp = true;
    if (strlen_c - strlen_g > 4) {
      FPRINTF(fp,"ins%d+%daa ",aapos,strlen_c-strlen_g);
#if 0
      print_nt_p = false;
#endif
    } else if (aa_genomicseg[0] == aa_queryseq[0]) {
      FPRINTF(fp,"ins%s%d%s ",&(aa_genomicseg[1]),aapos,&(aa_queryseq[1]));
    } else {
      FPRINTF(fp,"ins%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
    }
  } else if (aa_genomicseg[0] == 'X' || aa_queryseq[0] == 'X') {
#if 0
    print_nt_p = false;
#endif
    print_refquerypos_p = false;
  } else {
    if (*printp == true) FPRINTF(fp,", "); else *printp = true;
    FPRINTF(fp,"%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
  }
#if 0
  if (print_nt_p == true) {
    FPRINTF(fp,"(%s>%s) ",nt_genomicseg,nt_queryseq);
  }
#endif
  if (print_refquerypos_p == true) {
#ifdef PMAP
    FPRINTF(fp,"[%d]",refquerypos+2);
#else    
    FPRINTF(fp,"[%d]",refquerypos);
#endif
  }
  
  return;
}

static void
print_large_deletion (Filestring_T fp, bool *printp, int lastaapos, int nextaapos, int refquerypos) {

  if (*printp == true) FPRINTF(fp,", "); else *printp = true;
  FPRINTF(fp,"del%d-%daa ",lastaapos+1,nextaapos-lastaapos-1);
  FPRINTF(fp,"[%d]",refquerypos+3);
  
  return;
}


void
Translation_print_comparison (Filestring_T fp, struct Pair_T *pairs, int npairs, struct Pair_T *refpairs, int nrefpairs,
			      int cdna_direction, int relaastart, int relaaend, int maxmutations) {
  int i, j;
  int aapos, strlen_g, strlen_c, netchars;
  bool printp = false;
  char nt_genomicseg[MAXMUT], aa_genomicseg[MAXMUT], nt_queryseq[MAXMUT], aa_queryseq[MAXMUT];

  FPRINTF(fp,"    Amino acid changes: ");

  if (relaastart < relaaend) {
    if ((aapos = pairs[i=0].aapos) == 0) {
      i = next_aapos_fwd(pairs,0,npairs,0);
    }
    while (i < npairs) {
      aapos = pairs[i].aapos;
      j = next_aapos_fwd(pairs,i,npairs,aapos);
      if (pairs[i].aa_g != ' ' && pairs[i].aa_e != ' ') {
	/* not a frameshift */
	fill_aa_fwd(&strlen_g,&strlen_c,&netchars,aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq,
		    &(pairs[i]),&(pairs[j-1]));
	if (strcmp(aa_genomicseg,aa_queryseq) && strcmp(nt_genomicseg,nt_queryseq)) {
	  if (netchars % 3 == 0 || netchars > 12 || netchars < -12) {
	    print_mutation(fp,&printp,aapos,strlen_g,strlen_c,pairs[i].refquerypos,
			   aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq);
	  }
	} else if (j < npairs && pairs[j].aapos - aapos > 4) {
	  print_large_deletion(fp,&printp,aapos,pairs[j].aapos,pairs[i].refquerypos);
	}
      }
      i = j;
    }

  } else {
    if ((aapos = pairs[i=npairs-1].aapos) == 0) {
      i = next_aapos_rev(pairs,0,npairs,0);
    }
    while (i >= 0) {
      aapos = pairs[i].aapos;
      j = next_aapos_rev(pairs,i,npairs,aapos);
      if (pairs[i].aa_g != ' ' && pairs[i].aa_e != ' ') {
	/* not a frameshift */
	fill_aa_rev(&strlen_g,&strlen_c,&netchars,aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq,
		    &(pairs[i]),&(pairs[j+1]));
	if (strcmp(aa_genomicseg,aa_queryseq) && strcmp(nt_genomicseg,nt_queryseq)) {
	  if (netchars % 3 == 0 || netchars > 12 || netchars < -12) {
	    print_mutation(fp,&printp,aapos,strlen_g,strlen_c,pairs[i].refquerypos,
			   aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq);
	  }
	} else if (j >= 0 && pairs[j].aapos - aapos > 4) {
	  print_large_deletion(fp,&printp,aapos,pairs[j].aapos,pairs[i].refquerypos);
	}
      }
      i = j;
    }
  }

  FPRINTF(fp,"\n");

  return;
}


