package HelminthGenomeAnalysis::AvrilFastaUtils;

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Moose;
use Math::Round; # HAS THE nearest() FUNCTION
use Carp::Assert; # HAS THE assert() FUNCTION 
use Scalar::Util qw(looks_like_number);

has 'fasta_file'   => (is => 'ro', isa => 'Str', required => 1);
has '_contigs2seq' => (is => 'ro', isa => 'HashRef[Str]', lazy => 1, builder => "_build__contigs2seq");

use base 'Exporter';
our @EXPORT_OK = qw( build_contigs2len count_internal_stops print_seq_to_fasta count_number_seqs get_fasta_stats );

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: _build__contigs2seq: SUBROUTINE TO READ IN SEQUENCES OF SCAFFOLDS INTO A HASH %contig2seq

sub _build__contigs2seq
{
    my ($self)               = @_;

    # THROW AN EXCEPTION IF THE FASTA FILE DOES NOT EXIST:
    throw Error::Simple("ERRORCODE=1: _build__contig2seq: fasta file does not exist") if !(-e $self->fasta_file);

    my $fasta_file           = $self->fasta_file; # FASTA FILE
    my %contigs2seq;                              # HASH TABLE OF CONTIG/SCAFFOLD SEQUENCES
    my %SEEN                 = ();                # HASH TABLE OF SEQUENCES WE HAVE SEEN ALREADY
    my $name;                                     # SEQUENCE NAME
    my $num_seqs             = 0;                 # NUMBER OF SEQUENCES IN THE FASTA FILE

    # READ IN EACH SEQUENCE FROM THE FASTA FILE AT A TIME:
    my $seqin                = Bio::SeqIO->new( '-format' => 'Fasta', '-file' => $fasta_file);
    while ( $_ = $seqin->next_seq() ) 
    {
        $contigs2seq{$_->display_id()} = $_->seq();
        $name                = $_->display_id();
        # THROW AN EXCEPTION IF A SEQUENCE NAME APPEARS TWICE IN THE FASTA FILE (bioperl DIDN'T PICK THIS UP):
        throw Error::Simple("ERRORCODE=2: _build__contig2seq: sequence $name appears twice in $fasta_file") if ($SEEN{$name});
        $SEEN{$name}         = 1;
        $num_seqs++;
    }
    # THROW AN EXCEPTION IF NO SEQUENCES WERE READ IN FROM THE FASTA FILE:
    throw Error::Simple("ERRORCODE=3: _build__contig2seq: no sequences in $fasta_file") if $num_seqs == 0;

    return(\%contigs2seq);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: build_contigs2len: MAKE A HASH TABLE WITH THE LENGTHS OF SEQUENCES IN A FASTA FILE.

sub build_contigs2len
{
   my ($self)                = @_;
   my %contigs2len;                                 # HASH TABLE OF CONTIG/SCAFFOLD LENGTHS 
   my %SEEN                  = ();                  # HASH TABLE TO RECORD NAMES OF SCAFFOLDS WE HAVE SEEN 
   my $contigs2seq           = $self->_contigs2seq; # HASH TABLE OF SEQUENCES OF SCAFFOLDS/CONTIGS
  
   # GET THE LENGTH OF EACH SEQUENCE:
   foreach my $name (keys %$contigs2seq)
   {
      assert(!(defined($SEEN{$name}))); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $SEEN{$name}           = 1;
      my $seq                = $contigs2seq->{$name};
      my $seqlength          = length($seq);
      assert(looks_like_number($seqlength)); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      assert($seqlength >= 1); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      assert(!($seqlength =~ /\D/)); # $seqlength IS NOT AN INTEGER; THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $contigs2len{$name}    = $seqlength;
        
   }
 
   return(\%contigs2len);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: count_internal_stops(): CHECK WHETHER ANY OF THE PROTEIN SEQUENCES IN A FASTA FILE HAVE INTERNAL STOP CODONS. IT ASSUMES STOP CODONS ARE REPRESENTED BY X OR *. IT RETURNS A HASH TABLE WITH THE NUMBER OF INTERNAL STOP CODONS IN SEQUENCES THAT HAVE AT LEAST ONE INTERNAL STOP CODON.

sub count_internal_stops
{
   my ($self)                = @_;
   my %SEEN                  = ();                  # HASH TABLE TO RECORD NAMES OF SCAFFOLDS WE HAVE SEEN 
   my $contigs2seq           = $self->_contigs2seq; # HASH TABLE OF SEQUENCES OF SCAFFOLDS/CONTIGS
   my $internal_stops;                              # NUMBER OF INTERNAL STOP CODONS IN A SEQUENCE 
   my %protname2internalstops    = ();              # HASH TABLE OF THE NUMBER INTERNAL STOP CODONS IN SEQUENCES. 

   # LOOK AT EACH SEQUENCE IN TURN:
   foreach my $name (keys %$contigs2seq)
   {
      assert(!(defined($SEEN{$name}))); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $SEEN{$name}           = 1;
      my $seq                = $contigs2seq->{$name};
      $seq                   =~ tr/[a-z]/[A-Z]/;
      # IGNORE THE STOP CODON AT THE END OF THE SEQUENCE (IF THERE IS ONE):
      if (substr($seq,length($seq)-1,1) eq '*') { chop($seq);}
      # CHECK WHETHER THERE IS AN INTERNAL STOP CODON IN THE SEQUENCE:
      $internal_stops        = ($seq =~ tr/\*//); 
      assert(looks_like_number($internal_stops)); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      assert($internal_stops >= 0); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      assert(!($internal_stops =~ /\D/)); # $internal_stops IS NOT AN INTEGER; THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.      
      if ($internal_stops >= 1)
      {
         assert(!(defined($protname2internalstops{$name}))); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
         $protname2internalstops{$name} = $internal_stops;
      }
   }
 
   return(\%protname2internalstops);

}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: print_seq_to_fasta(): PRINTS A SEQUENCE TO A FASTA OUTPUT FILE, IN LINES OF LENGTH 60. THIS CAN EITHER APPEND TO AN EXISTING OUTPUT FILE, OR WRITE A NEW FILE.

sub print_seq_to_fasta
{
   my $fasta                 = $_[0]; # FASTA FILE
   my $seq                   = $_[1]; # SEQUENCE
   my $name                  = $_[2]; # SEQUENCE NAME
   my $append                = $_[3]; # SAYS WHETHER TO APPEND TO AN EXISTING FILE OR NOT (yes/no)
   my $length;                        # LENGTH OF THE SEQUENCE  
   my $offset;                        # OFFSET USED WHILE PRINTING OUT THE SEQUENC
   my $a_line;                        # A LINE OF THE SEQUENC
 
   # THROW AN EXCEPTION IF $append IS NOT yes/no:
   throw Error::Simple("ERRORCODE=1: print_seq_to_file: append is $append, should be yes/no") if ($append ne 'yes' && $append ne 'no');

   # OPEN THE OUTPUT FASTA FILE:
   if    ($append eq 'yes') { open(FASTA,">>$fasta") || die "ERROR: print_seq_to_file: cannot open $fasta\n"; }
   elsif ($append eq 'no')  { open(FASTA,">$fasta") || die "ERROR: print_seq_to_file: cannot open $fasta\n";  }
   
   print FASTA ">$name\n";
   $length                 = length($seq);
   $offset                 = 0;
   while ($offset < $length)
   {
      $a_line              = substr($seq,$offset,60);
      print FASTA "$a_line\n";
      $offset              = $offset + 60;
   }
   # CLOSE THE OUTPUT FILE:
   close(FASTA);   

   return(1); 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: count_number_seqs(): COUNT THE NUMBER OF SEQUENCES IN A FASTA FILE:

sub count_number_seqs
{
   my ($self)                = @_;
   my $contigs2seq           = $self->_contigs2seq; # HASH TABLE OF SEQUENCES OF SCAFFOLDS/CONTIGS
   my $number_seqs           = 0;                   # NUMBER OF SEQUENCES IN THE FASTA FILE 
   my %SEEN                  = ();                  # HASH TABLE TO RECORD NAMES OF SEQUENCES WE HAVE SEEN 

   # LOOK AT EACH SEQUENCE IN TURN:
   foreach my $name (keys %$contigs2seq)
   {
      assert(!(defined($SEEN{$name}))); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $SEEN{$name}           = 1;
      $number_seqs++;
   }
 
   return($number_seqs);
  
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: get_fasta_stats(): SUBROUTINE TO CALCULATE THE NUMBER OF SEQUENCES, AVERAGE SEQUENCE LENGTH AND MAXIMUM LENGTH FOR A FASTA FILE:

sub get_fasta_stats
{
   my $file                  = $_[0]; # INPUT FASTA FILE
   my $cmd;                           # COMMAND TO RUN
   my $num_seqs              = -1;    # NUMBER OF SEQUENCES
   my $av_seq_len            = -1;    # AVERAGE SEQUENCE LENGTH
   my $max_seq_len           = -1;    # MAXIMUM SEQUENCE LENGTH
   my $tot_seq_len           = -1;    # TOTAL SEQUENCE LENGTH 
   my $line;                          # 
   my @temp;                          # 

   $cmd                      = "~mh12/bin/stats $file";
   open(TEMP,"$cmd |");
   while(<TEMP>)
   {
      $line                  = $_;
      chomp $line;
      if ($line =~ /sum/)
      {
         # sum = 4920469, n = 27189, ave = 180.97, largest = 6633
         @temp               = split(/\s+/,$line);
         $tot_seq_len        = $temp[2]; chop($tot_seq_len);# REMOVE COMMA
         $num_seqs           = $temp[5]; chop($num_seqs);   # REMOVE COMMA
         $av_seq_len         = $temp[8]; chop($av_seq_len); # REMOVE COMMA
         $max_seq_len        = $temp[11];
      } 
   }
   close(TEMP);
   assert(looks_like_number($tot_seq_len)); assert($tot_seq_len != -1);
   assert(looks_like_number($num_seqs)); assert($num_seqs != -1);
   assert(looks_like_number($av_seq_len)); assert($av_seq_len != -1);
   assert(looks_like_number($max_seq_len)); assert($max_seq_len != -1);

   return($num_seqs, $av_seq_len, $max_seq_len, $tot_seq_len);
}

#------------------------------------------------------------------#

1;
