#ifndef MAXENT_INCLUDED
#define MAXENT_INCLUDED

#define DONOR_MODEL_LEFT_MARGIN 3 /* Amount in exon.  Does not include GT */
#define DONOR_MODEL_RIGHT_MARGIN 6 /* Amount in intron */

#define ACCEPTOR_MODEL_LEFT_MARGIN 20 /* Amount in intron.  Includes AG */
#define ACCEPTOR_MODEL_RIGHT_MARGIN 3 /* Amount in exon */

#define MAXENT_MAXLENGTH 24	/* larger of donor left+right and acceptor left+right, plus 1 */

extern double
Maxent_donor_prob (char *sequence);
extern double
Maxent_donor_prob_revcomp (char *sequence);
extern double
Maxent_donor_prob_nucleotides (unsigned char *nucleotides);
extern double
Maxent_donor_logodds (char *sequence);
extern double
Maxent_donor_logodds_nucleotides (unsigned char *nucleotides);

extern double
Maxent_acceptor_prob (char *sequence);
extern double
Maxent_acceptor_prob_revcomp (char *sequence);
extern double
Maxent_acceptor_prob_nucleotides (unsigned char *nucleotides);
extern double
Maxent_acceptor_logodds (char *sequence);
extern double
Maxent_acceptor_logodds_nucleotides (unsigned char *nucleotides);

#endif

