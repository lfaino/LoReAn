package HelminthGenomeAnalysis::AvrilAlignUtils;

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Moose;
use Math::Round; # HAS THE nearest() FUNCTION
use Carp::Assert; # HAS THE assert() FUNCTION 
use Scalar::Util qw(looks_like_number);

use base 'Exporter';
our @EXPORT_OK = qw( run_blast format_blastdb record_blast_hits merge_overlapping_blast_hits );

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: format_blastdb(): GIVEN A FASTA FILE OF SEQUENCES, FORMATS IT AS A BLAST DATABASE:

sub format_blastdb
{
   my $blast_path            = $_[0]; # PATH TO THE BLAST EXECUTABLE, eg. /software/pubseq/bin/ncbi_blast+/
   my $fasta_file            = $_[1]; # FASTA FILE TO MAKE A BLAST DATABASE FOR
   my $fasta_file_type       = $_[2]; # TYPE OF SEQUENCE IN $fasta_file: prot/nucl
   my $cmd;                           # COMMAND TO RUN
   my $systemcall;                    # RETURN VALUE FROM THE SYSTEM CALL 

   # THROW AN EXCEPTION IF $fasta_file_type IS NOT prot/nucl
   throw Error::Simple("ERRORCODE=1: format_blastdb: fasta_file_type is $fasta_file_type, should be prot/nucl") if ($fasta_file_type ne "prot" && $fasta_file_type ne "nucl");

   # FORMAT THE BLAST DATABASE:   
   if ($fasta_file_type eq 'nucl')
   {
      $cmd                   = $blast_path."/makeblastdb -in $fasta_file -input_type fasta -dbtype $fasta_file_type -parse_seqids -out $fasta_file";
   }
   elsif ($fasta_file_type eq 'prot')
   {
      $cmd                   = $blast_path."/makeblastdb -in $fasta_file -input_type fasta -dbtype $fasta_file_type -parse_seqids -out $fasta_file";
   }

   $systemcall               = system "$cmd";
   # THROW AN EXCEPTION IF $systemcall IS NOT 0:
   throw Error::Simple("ERRORCODE=2: format_blastdb: return value $systemcall when tried to run: $cmd") if ($systemcall != 0);
   sleep(1); # WAIT A SECOND TO CHECK THE OUTPUT FILE IS FINISHED

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: run_blast(): GIVEN A QUERY SEQUENCE FILE, AND DATABASE FILE, RUNS BLAST:

sub run_blast
{
   my $blast_path            = $_[0]; # PATH TO THE BLAST EXECUTABLE, eg. /software/pubseq/bin/ncbi_blast+/
   my $blast_output          = $_[1]; # BLAST OUTPUT FILE NAME
   my $blast_type            = $_[2]; # TYPE OF BLAST TO RUN, eg. tblastn, blastp, blastn, blastx 
   my $db_file               = $_[3]; # DATABASE TO SEARCH AGAINST
   my $query_file            = $_[4]; # QUERY FILE 
   my $cmd;                           # COMMAND TO RUN
   my $systemcall;                    # RETURN VALUE FROM THE SYSTEM CALL 

   # THROW AN EXCEPTION IF $blast_type IS NOT blastp/blastn/blastx/tblastn:
   throw Error::Simple("ERRORCODE=1: run_blast: blast_type is $blast_type, should be blastp/blastn/tblastn/blastx") if ($blast_type ne 'blastp' && $blast_type ne 'blastn' && $blast_type ne 'blastx' && $blast_type ne 'tblastn');
 
   # RUN BLAST:
   # -outfmt 7 GIVES TABULAR OUTPUT, SEE http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.Obtaining_Sample_data
   if    ($blast_type eq 'blastp')
   {
      $cmd                   = $blast_path."/blastp -db $db_file -query $query_file -out $blast_output -outfmt 7";
   }
   elsif ($blast_type eq 'blastn')
   {
      $cmd                   = $blast_path."/blastn -db $db_file -query $query_file -out $blast_output -outfmt 7";
   }
   elsif ($blast_type eq 'blastx')
   {
      $cmd                   = $blast_path."/blastx -db $db_file -query $query_file -out $blast_output -outfmt 7";
   }
   elsif ($blast_type eq 'tblastn')
   {
      $cmd                   = $blast_path."/tblastn -db $db_file -query $query_file -out $blast_output -outfmt 7";
   }
   $systemcall               = system "$cmd";
   # THROW AN EXCEPTION IF $systemcall IS NOT 0:
   throw Error::Simple("ERRORCODE=2: run_blast: return value $systemcall when tried to run: $cmd") if ($systemcall != 0);
   sleep(1); # WAIT A SECOND TO CHECK THE OUTPUT FILE IS FINISHED

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: run_exonerate(): RUN EXONERATE BETWEEN A SCAFFOLD AND A PROTEIN:

sub run_exonerate
{
   my $exonerate_output      = $_[0]; # EXONERATE OUTPUT FILE
   my $scaffold_seq_file     = $_[1]; # FASTA FILE OF SCAFFOLDS
   my $input_pep_file        = $_[2]; # PROTEIN SEQUENCE FILE
   my $outputdir             = $_[3]; # DIRECTORY TO WRITE OUTPUT FILES IN 
   my $type                  = $_[4]; # 'est' OR 'prot' - SAYS WHETHER THE QUERY SEQUENCE IS EST OR PROTEIN 
   my $cmd;                           # COMMAND TO RUN
   my $systemcall;                    # RETURN VALUE FROM THE SYSTEM CALL 
   my $temp_exonerate_output;         # A TEMPORARY EXONERATE OUTPUT FILE 
   my $line;                          # 
   my @temp;                          # 
   my $gff_start             = 0;     # SAYS WHETHER WE FOUND THE START OF THE GFF FILE 

   # THROW AN EXCEPTION IF $type IS NOT 'est' OR 'prot':
   throw Error::Simple("ERRORCODE=2: run_exonerate: type is $type, but should be est/prot") if ($type ne 'est' && $type ne 'prot'); # TESTED FOR

   $temp_exonerate_output    = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);  
   open(TEMP_EXONERATE_OUTPUT,">$temp_exonerate_output") || die "ERROR: run_exonerate: cannot open $temp_exonerate_output\n";
   close(TEMP_EXONERATE_OUTPUT);
   if ($type eq 'prot')
   {
      $cmd                   = "exonerate -t $scaffold_seq_file -q $input_pep_file --model protein2genome --bestn 1 --showtargetgff > $temp_exonerate_output";
   }
   elsif ($type eq 'est')
   {
      $cmd                   = "exonerate -t $scaffold_seq_file -q $input_pep_file --model coding2genome --bestn 1 --showtargetgff > $temp_exonerate_output";
   }
   $systemcall               = system "$cmd";
   # PRINT A WARNING MESSAGE IF $systemcall IS NOT 0:
   # NOTE THAT SOMETIMES EXONERATE DOES FAIL FOR SOME INPUT SEQUENCES (eg. WITH AN ERROR MESSAGE ABOUT 'Bad seed'), SO WE PROBABLY JUST WANT TO PRINT A WARNING AND CARRY ON TO THE NEXT 
   # INPUT SEQUENCES, RATHER THAN HAVE THE SCRIPT DIE:
   if ($systemcall != 0)
   {
      print STDERR "WARNING: run_exonerate: failed to run exonerate for this query and target!\n"; 
   } 
   sleep(1); # WAIT A SECOND TO CHECK THE OUTPUT FILE IS FINISHED

   # TAKE THE GFF LINES FROM $temp_exonerate_output, AND PUT THEM IN $exonerate_output:
   open(OUTPUT,">$exonerate_output") || die "ERROR: run_exonerate: cannot open $exonerate_output\n";
   open(EXONERATE,"$temp_exonerate_output") || die "ERROR: run_exonerate: cannot open $temp_exonerate_output\n";
   while(<EXONERATE>)
   {
      $line                  = $_;
      chomp $line;
      @temp                  = split(/\t+/,$line);
      if ($line eq '##gff-version 2')           { $gff_start = 1;}
      if ($line eq '# --- END OF GFF DUMP ---') { $gff_start = 0;}
      if ($gff_start == 1 && substr($line,0,1) ne '#')
      {
         # THE EXONERATE GFF SHOULD HAVE 8 OR 9 COLUMNS:
         assert($#temp == 7 || $#temp == 8);
         if    ($#temp == 7) # cds LINES DO NOT HAVE 9 COLUMNS
         {
            print OUTPUT "$line";
            print OUTPUT "cds\n";
         }
         elsif ($#temp == 8)
         {
            print OUTPUT "$line\n";
         }
      }
   }
   close(EXONERATE);
   close(OUTPUT);

   # DELETE TEMPORARY FILES:
   system "rm -f $temp_exonerate_output";   
    
   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: record_blast_hits(): GIVEN A BLAST OUTPUT FILE, RECORDS THE BLAST HITS IN A HASH:

sub record_blast_hits
{
   my $blast_output          = $_[0]; # BLAST OUTPUT FILE
   my $evalue_cutoff         = $_[1]; # EVALUE CUTOFF
   my $flank_length          = $_[2]; # FLANK SIZE TO TAKE ON EITHER SIDE OF EACH BLAST HIT IN A SUBJECT
   my %BLAST                 = ();    # HASH TABLE OF BLAST HITS
   my $line;                          # 
   my @temp;                          # 
   my $subject;                       # BLAST HIT
   my $subject_start;                 # START OF THE HIT IN $subject  
   my $subject_end;                   # END OF THE HIT IN $subject
   my $tmp;                           #  
   my $evalue;                        # EVALUE FOR THE HIT
   my $pos;                           # POSITION OF THE BLAST HIT IN $subject
   my %SEEN                  = ();    # HASH TABLE OF BLAST HITS THAT WE HAVE SEEN ALREADY
   my $hits;                          # HITS FOR A PARTICULAR SUBJECT
   my $new_hits;                      # NON-OVERLAPPING HITS FOR A PARTICULAR SUBJECT

   open(BLAST,"$blast_output") || die "ERROR: cannot open $blast_output\n";
   while(<BLAST>)
   {
      $line                  = $_;
      chomp $line;
      if (substr($line,0,1) ne '#')
      {
         @temp                  = split(/\t+/,$line);
         # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
         $subject            = $temp[1];
         $subject_start      = $temp[8];
         $subject_end        = $temp[9];
         $evalue             = $temp[10];
         if ($subject_end < $subject_start)
         {
            $tmp             = $subject_end;
            $subject_end     = $subject_start;
            $subject_start   = $tmp;
         }
         if ($evalue <= $evalue_cutoff)
         {
            # TAKE $flank_length ON EITHER SIDE OF EACH BLAST HIT:
            $subject_start   = $subject_start - $flank_length;
            $subject_end     = $subject_end + $flank_length;
            $pos             = $subject_start."=".$subject_end;  
            if (!($SEEN{$subject."=".$pos}))
            {
               $SEEN{$subject."=".$pos} = 1;
               if (!($BLAST{$subject})) { $BLAST{$subject} = $pos;}
               else {$BLAST{$subject} = $BLAST{$subject}.",".$pos;}
            }
         }
      }
   }
   close(BLAST);

   # MERGE OVERLAPPING HITS IN ANY SUBJECT:
   foreach $subject (keys %BLAST)
   {
      $hits                  = $BLAST{$subject};
      $new_hits              = HelminthGenomeAnalysis::AvrilAlignUtils::merge_overlapping_blast_hits($hits);
      $BLAST{$subject}       = $new_hits;
   }

   return(\%BLAST);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS merge_overlapping_blast_hits(): GIVEN A LIST OF START-END COORDINATES OF BLAST HITS, MERGES OVERLAPPING HITS:

sub merge_overlapping_blast_hits
{
   my $hits                  = $_[0]; # POSITIONS OF BLAST HITS, eg. 50-100,99-105,107-109
   my @hits;                          # ARRAY OF HITS 
   my $i;                             #    
   my @temp;                          # 
   my $hit_start;                     # START OF A HIT
   my $hit_end;                       # END OF A HIT
   my %START                 = ();    # HASH TABLE OF STARTS OF HITS
   my $hit;                           # A HIT
   my $prev_hit              = "none";# THE PREVIOUS HIT
   my $prev_hit_start;                # START OF THE PREVIOUS HIT
   my $prev_hit_end;                  # END OF THE PREVIOUS HIT
   my $new_hits              = "";    # THE FINAL MERGED HITS

   # SORT THE HITS BY START POINT:
   @hits                     = split(/\,/,$hits);
   for ($i = 1; $i <= $#hits + 1; $i++)
   {
      $hit                   = $hits[($i-1)];
      @temp                  = split(/=/,$hit);
      $hit_start             = $temp[0];
      # THROW AN EXCEPTION IF WE ALREADY HAVE STORED THE START OF $hit:
      throw Error::Simple("ERRORCODE=1: merge_overlapping_blast_hits: already have stored start for $hit") if (defined($START{$hit}));
      $START{$hit}           = $hit_start;
   }
   @hits                     = reverse sort { $START{$b} <=> $START{$a} } keys %START;

   # MERGE OVERLAPPING HITS: 
   for ($i = 1; $i <= $#hits + 1; $i++)
   {
      $hit                   = $hits[($i-1)];
      @temp                  = split(/=/,$hit);
      $hit_start             = $temp[0];
      $hit_end               = $temp[1];
      # CHECK IF THIS HIT OVERLAPS WITH THE PREVIOUS HIT:
      if ($i > 1)
      {
         if (($hit_start >= $prev_hit_start && $hit_start <= $prev_hit_end) ||
             ($hit_end   >= $prev_hit_start && $hit_end   <= $prev_hit_end) ||
             ($prev_hit_start >= $hit_start && $prev_hit_start <= $hit_end) ||
             ($prev_hit_end   >= $hit_start && $prev_hit_end   <= $hit_end))
         {
            # THIS HIT OVERLAPS WITH THE PREVIOUS HIT, SO MERGE THEM:
            assert($hit_start >= $prev_hit_start);
            if ($hit_end > $prev_hit_end) { $prev_hit_end = $hit_end;}
         }
         else # THIS HIT DOES NOT OVERLAP WITH THE PREVIOUS HIT, SO RECORD THE PREVIOUS HIT:
         {
            $prev_hit        = $prev_hit_start."=".$prev_hit_end;
            $new_hits        = $new_hits.",".$prev_hit; 
            $prev_hit_start  = $hit_start;
            $prev_hit_end    = $hit_end;
         }
      }
      elsif ($i == 1)
      {
         $prev_hit_start     = $hit_start;
         $prev_hit_end       = $hit_end;
      }
   }
   $prev_hit                 = $prev_hit_start."=".$prev_hit_end;
   $new_hits                 = $new_hits.",".$prev_hit; 
   assert($new_hits ne '');
   $new_hits               = substr($new_hits,1,length($new_hits)-1);
 
   return($new_hits);
}

#------------------------------------------------------------------#

1;
