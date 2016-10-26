static char rcsid[] = "$Id: intron.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intron.h"
#include <stdlib.h>		/* For abort() */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


int
Intron_type (char left1, char left2, char right2, char right1,
	     char left1_alt, char left2_alt, char right2_alt, char right1_alt,
	     int cdna_direction
#ifdef INTRON_HELP
	     , IIT_T splicesites_iit, int *splicesites_divint_crosstable,
	     int donor_typeint, int acceptor_typeint, Chrnum_T chrnum,
	     Chrpos_T leftgenomepos, Chrpos_T rightgenomepos,
	     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
	     ) {
  int introntype, leftdi, rightdi;

  if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
    leftdi = LEFT_GT;
  } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
    leftdi = LEFT_GC;
  } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
    leftdi = LEFT_AT;
#ifndef PMAP
  } else if ((left1 == 'C' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
    leftdi = LEFT_CT;
#endif

#ifdef INTRON_HELP
    /* Not tested */
  } else if (splicesites_iit == NULL) {
    debug(printf("splicesites_iit is NULL\n"));
    return NONINTRON;
  } else if (cdna_direction > 0) {
    if (watsonp) {
      debug(printf("cdna_direction %d, watsonp %d, looking for donor at %u..%u, sign +1\n",
		   cdna_direction,watsonp,chroffset+leftgenomepos,chroffset+leftgenomepos+1U));
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chroffset+leftgenomepos,chroffset+leftgenomepos+1U,
					     donor_typeint,/*sign*/+1) == true) {
	leftdi = LEFT_GT;
      } else {
	return NONINTRON;
      }
    } else {
      debug(printf("cdna_direction %d, watsonp %d, looking for donor at %u..%u, sign -1\n",
		   cdna_direction,watsonp,chrhigh-leftgenomepos,chrhigh-leftgenomepos+1U));
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chrhigh-leftgenomepos,chrhigh-leftgenomepos+1U,
					     donor_typeint,/*sign*/-1) == true) {
	leftdi = LEFT_GT;
      } else {
	return NONINTRON;
      }
    }
  } else if (cdna_direction < 0) {
    if (watsonp) {
      debug(printf("cdna_direction %d, watsonp %d, looking for acceptor at %u..%u, sign -1\n",
		   cdna_direction,watsonp,chroffset+leftgenomepos,chroffset+leftgenomepos+1U));
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chroffset+leftgenomepos,chroffset+leftgenomepos+1U,
					     acceptor_typeint,/*sign*/-1) == true) {
	leftdi = LEFT_CT;
      } else {
	return NONINTRON;
      }
    } else {
      debug(printf("cdna_direction %d, watsonp %d, looking for acceptor at %u..%u, sign +1\n",
		   cdna_direction,watsonp,chrhigh-leftgenomepos,chrhigh-leftgenomepos+1U));
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chrhigh-leftgenomepos,chrhigh-leftgenomepos+1U,
					     acceptor_typeint,/*sign*/+1) == true) {
	leftdi = LEFT_CT;
      } else {
	return NONINTRON;
      }
    }
#endif

  } else {
    return NONINTRON;
  }

  if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
    rightdi = RIGHT_AG;
  } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
    rightdi = RIGHT_AC;
#ifndef PMAP
  } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
    rightdi = RIGHT_GC;
  } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
    rightdi = RIGHT_AT;
#endif

#ifdef INTRON_HELP
    /* Not tested */
  } else if (splicesites_iit == NULL) {
    return NONINTRON;
  } else if (cdna_direction > 0) {
    if (watsonp) {
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chroffset+rightgenomepos,chroffset+rightgenomepos+1U,
					     acceptor_typeint,/*sign*/+1) == true) {
	rightdi = RIGHT_AG;
      } else {
	return NONINTRON;
      }
    } else {
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chrhigh-rightgenomepos,chrhigh-rightgenomepos+1U,
					     acceptor_typeint,/*sign*/-1) == true) {
	rightdi = RIGHT_AG;
      } else {
	return NONINTRON;
      }
    }
  } else if (cdna_direction < 0) {
    if (watsonp) {
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chroffset+rightgenomepos,chroffset+rightgenomepos+1U,
					     donor_typeint,/*sign*/-1) == true) {
	rightdi = RIGHT_AC;
      } else {
	return NONINTRON;
      }
    } else {
      if (IIT_exists_with_divno_typed_signed(splicesites_iit,splicesites_divint_crosstable[chrnum],
					     chrhigh-rightgenomepos,chrhigh-rightgenomepos+1U,
					     donor_typeint,/*sign*/-1) == true) {
	rightdi = RIGHT_AC;
      } else {
	return NONINTRON;
      }

    }
#endif

  } else {
    return NONINTRON;
  }


  if ((introntype = leftdi & rightdi) == 0x00) {
    return NONINTRON;
  } else if (cdna_direction > 0) {
    if (introntype < 0x08) {
      return NONINTRON;
    } else {
      return introntype;
    }
  } else if (cdna_direction < 0) {
    if (introntype > 0x04) {
      return NONINTRON;
    } else {
      return introntype;
    }
  } else {
    /* Should happen only from Stage3_merge_local_splice */
    /* return NONINTRON; */
    return introntype;		/* Needed for guess */
  }
}


char *
Intron_type_string (int introntype) {
  switch (introntype) {
  case GTAG_FWD: return "GT-AG, fwd";
  case GCAG_FWD: return "GC-AG, fwd";
  case ATAC_FWD: return "AT-AC, fwd";
#ifndef PMAP
  case GTAG_REV: return "GT-AG, rev";
  case GCAG_REV: return "GC-AG, rev";
  case ATAC_REV: return "AT-AC, rev";
#endif
  default: return "nonintron";
  }
}    

char *
Intron_left_dinucl_string (int dinucl) {
  switch (dinucl) {
  case LEFT_GT: return "GT-";
  case LEFT_GC: return "GC-";
  case LEFT_AT: return "AT-";
#ifndef PMAP
  case LEFT_CT: return "CT-";
#endif
  default: return "XX-";
  }
}

char *
Intron_right_dinucl_string (int dinucl) {
  switch (dinucl) {
  case RIGHT_AG: return "-AG";
  case RIGHT_AC: return "-AC";
#ifndef PMAP
  case RIGHT_GC: return "-GC";
  case RIGHT_AT: return "-AT";
#endif
  default: return "-XX";
  }
}


bool
Intron_canonical_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'G') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_canonical_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'C' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_gcag_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'C' &&
      acceptor2 == 'A' && acceptor1 == 'G') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_atac_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'A' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_gcag_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'C' && donor2 == 'T' &&
      acceptor2 == 'G' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_atac_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'T') {
    return true;
  } else {
    return false;
  }
}

