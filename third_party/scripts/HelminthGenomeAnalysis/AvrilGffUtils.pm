package HelminthGenomeAnalysis::AvrilGffUtils;

use strict;
use warnings;
use Math::Round; # HAS THE nearest() FUNCTION
use Carp::Assert; # HAS THE assert() FUNCTION 
use Scalar::Util qw(looks_like_number);

use base 'Exporter';
our @EXPORT_OK = qw( read_gene_ids replace_gene_ids_in_gff get_gff_lines_for_feature_types add_flanking_region_to_gff_features make_fasta_for_gff make_gff_for_features_in_subsequences read_genenames_for_transcripts get_gene_name add_gene_features_to_gff_for_subsequences make_gff_for_features_in_sequences sort_gff merge_overlapping_cds sort_gff_lines_for_genes count_features_in_region convert_exonerate_gff_to_standard_gff );

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: make_fasta_for_gff(): MAKE A FASTA FILE OF SEQUENCES SPECIFIED BY A GFF FILE.

sub make_fasta_for_gff
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $input_fasta         = $_[1]; # INPUT FASTA FILE
   my $output_fasta        = $_[2]; # OUTPUT FASTA FILE
   my $outputdir           = $_[3]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $cmd;                         # COMMAND TO RUN
   my $systemcall;                  # RETURN VALUE FROM SYSTEM CALL
   my $temp_output_fasta;           # TEMPORARY OUTPUT FASTA FILE 
   my $line;                        # 
   my @temp;                        # 
   my $region;                      # NAME OF REGION IN FASTA FILE
   my $scaffold;                    # SCAFFOLD OF REGION IN FASTA FILE
   my $start;                       # START OF REGION IN FASTA FILE IN $scaffold
   my $end;                         # END OF REGION IN FASTA FILE IN $scaffold

   $temp_output_fasta      = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); # MAKE A TEMPORARY FILE IN THE CURRENT DIRECTORY

   # MAKE A FASTA FILE WITH THE SEQUENCES SPECIFIED BY $input_gff:
   $cmd                    = "bedtools getfasta -fi $input_fasta -bed $input_gff -fo $temp_output_fasta";
   $systemcall             = system "$cmd";
   # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
   throw Error::Simple("ERRORCODE=1: add_flanking_region_to_gff_features: return value $systemcall when tried to run $cmd") if ($systemcall != 0); # TESTED FOR THIS
   system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED

   # NOTE: THE OUTPUT FASTA FILE $output_fasta HAS THE COORDINATES LIKE IN BED FILES: 
   # BED FILES ARE ZERO-BASED, HALF-OPEN: http://code.google.com/p/bedtools/wiki/FAQ#What_does_zero-based,_half-open_mean?
   # WE WANT TO RE-WRITE THESE COORDINATES SO THEY ARE 1-BASED:
   open(OUTPUT_FASTA,">$output_fasta") || die "ERROR: make_fasta_for_gff: cannot open $output_fasta\n";
   open(TEMP_OUTPUT_FASTA,"$temp_output_fasta") || die "ERROR: make_fasta_for_gff: cannot open $temp_output_fasta\n";
   while(<TEMP_OUTPUT_FASTA>)
   {
      $line                = $_;
      chomp $line;
      if (substr($line,0,1) eq '>')
      {
         @temp             = split(/\s+/,$line);
         $region           = $temp[0];
         $region           = substr($region,1,length($region)-1); # eg. scaff1:0-7
         @temp             = split(/:/,$region);
         $scaffold         = $temp[0]; # eg. scaff1
         $start            = $temp[1]; # eg. 0-7
         @temp             = split(/-/,$start);
         $start            = $temp[0]; # eg. 0
         $end              = $temp[1]; # eg. 7
         # CHANGE 0-BASED START COORDINATE TO A 1-BASED START COORDINATE:
         $start            = $start + 1;
         print OUTPUT_FASTA ">$scaffold:$start-$end\n";
      }
      else { print OUTPUT_FASTA "$line\n";}
   }
   close(TEMP_OUTPUT_FASTA);

   return(1);
 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: make_gff_for_features_in_subsequences(): GIVEN AN INPUT GFF OF FEATURES (eg. GENES) IN SCAFFOLDS, AND A SECOND INPUT GFF OF SUBSEQUENCES OF THE SCAFFOLDS, MAKES AN OUTPUT GFF OF FEATURES (eg. GENES) IN THE SUBSEQUENCES. 

sub make_gff_for_features_in_subsequences
{
   my $input_gff           = $_[0]; # INPUT GFF FILE OF FEATURES IN SCAFFOLDS
   my $input_gff2          = $_[1]; # INPUT GFF FILE OF SUBSEQUENCES OF SCAFFOLDS
   my $output_gff          = $_[2]; # OUTPUT GFF FILE OF FEATURES IN SUBSEQUENCES
   my $line;                        # 
   my @temp;                        # 
   my $feature;                     # FEATURE NAME
   my %START               = ();    # START POSITIONS OF FEATURES IN SCAFFOLDS
   my %END                 = ();    # END POSITIONS OF FEATURES IN SCAFFOLDS 
   my $feature_start_in_scaff;      # START OF FEATURE $feature IN THE SCAFFOLD 
   my $feature_end_in_scaff;        # END OF FEATURE $feature IN THE SCAFFOLD 
   my $subsequence_start_in_scaff;  # START OF THE SUBSEQUENCE IN THE SCAFFOLD
   my $subsequence_end_in_scaff;    # END OF THE SUBSEQUENCE IN THE SCAFFOLD 
   my $feature_start_in_subsequence;# START OF THE FEATURE $feature IN THE SUBSEQUENCE
   my $feature_end_in_subsequence;  # END OF THE FEATURE $feature IN THE SUBSEQUENCE
   my $scaffold;                    # SCAFFOLD NAME
   my $feature_length;              # FEATURE LENGTH 
   my $diff_in_start;               # DIFFERENCE IN THE START OF A FEATURE BETWEEN $input_gff AND $output_gff
   my %DIFF_IN_START       = ();    # HASH TABLE TO STORE $diff_in_start 
 
   # READ IN THE GFF FILE OF FEATURES IN SCAFFOLDS:
   open(INPUT_GFF,"$input_gff") || die "ERROR: make_gff_for_features_in_subsequences: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      # eg. scaff1  source	feature	20	30	.	+	.	feature1 
      @temp                = split(/\t+/,$line);
      $feature_start_in_scaff = $temp[3];
      $feature_end_in_scaff= $temp[4];
      $feature             = $temp[8];
      # THROW AN ERROR IF WE ALREADY HAVE STORED THE START POSITION OF $feature:
      throw Error::Simple("ERRORCODE=1: make_gff_for_features_in_subsequences: already stored start for $feature") if (defined($START{$feature})); # TESTED FOR
      $START{$feature}     = $feature_start_in_scaff;
      # THROW AN ERROR IF WE ALREADY HAVE STORED THE END POSITION OF $feature:
      throw Error::Simple("ERRORCODE=2: make_gff_for_features_in_subsequences: already stored end for $feature") if (defined($END{$feature})); # TESTED FOR
      $END{$feature}       = $feature_end_in_scaff; 
   }
   close(INPUT_GFF);   

   # OPEN THE OUTPUT GFF FILE:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: make_gff_for_features_in_subsequences: cannot open $output_gff\n";

   # READ IN THE GFF FILE OF SUBSEQUENCES IN SCAFFOLDS:
   open(INPUT_GFF2,"$input_gff2") || die "ERROR: make_gff_for_features_in_subsequences: cannot open $input_gff2\n";
   while(<INPUT_GFF2>)
   {
      $line                = $_;
      chomp $line;
      # eg. scaff1	source	feature	15	35	.	+	.	feature1 
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $subsequence_start_in_scaff = $temp[3];
      $subsequence_end_in_scaff = $temp[4]; 
      $feature             = $temp[8];
      # THROW AN ERROR IF WE DO NOT KNOW THE START OF THIS FEATURE IN THE SCAFFOLD:
      assert(defined($START{$feature})); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $feature_start_in_scaff = $START{$feature};
      # THROW AN ERROR IF WE DO NOT KNOW THE END OF THIS FEATURE IN THE SCAFFOLD:
      assert(defined($END{$feature})); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $feature_end_in_scaff = $END{$feature};
      $feature_length      = $feature_end_in_scaff - $feature_start_in_scaff;
      # FIND THE START AND END OF THE FEATURE WITH RESPECT TO THE SUBSEQUENCE:
      $feature_start_in_subsequence = $feature_start_in_scaff - $subsequence_start_in_scaff + 1; 
      # eg. feature_start_in_scaff=20, subsequence_start_in_scaff=15 ---> 20-15+1 = 6
      $feature_end_in_subsequence = $feature_start_in_subsequence + $feature_length;  
      # eg. feature_end_in_scaff = 30, feature_length = 10, feature_end_in_subsequence = 6 + 10 = 16
      print OUTPUT_GFF "$scaffold:$subsequence_start_in_scaff-$subsequence_end_in_scaff\t$temp[1]\t$temp[2]\t$feature_start_in_subsequence\t$feature_end_in_subsequence\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\n"; 
      # RECORD THE DIFFERENCE IN THE START OF THE FEATURE BETWEEN $input_gff AND $output_gff:
      $diff_in_start       = $feature_start_in_scaff - $feature_start_in_subsequence;   
      assert(!(defined($DIFF_IN_START{$feature}))); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.
      $DIFF_IN_START{$feature} = $diff_in_start;
   }
   close(INPUT_GFF2);
   close(OUTPUT_GFF);
 
   return(\%DIFF_IN_START); 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: add_flanking_region_to_gff_features(): ADD flank_size ON EITHER SIDE OF EACH GFF FEATURE IN AN INPUT GFF FILE, TO MAKE A NEW OUTPUT GFF.

sub add_flanking_region_to_gff_features
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $flank_size          = $_[1]; # FLANKING REGION SIZE
   my $scaffold_length_file= $_[2]; # FILE WITH LENGTHS OF SCAFFOLDS
   my $output_gff          = $_[3]; # OUTPUT GFF FILE
   my $cmd;                         # COMMAND TO RUN
   my $systemcall;                  # RETURN VALUE FROM SYSTEM CALL
   my $line;                        # 
   my @temp;                        # 

   # MAKE AN OUTPUT GFF BY ADDING A FLANKING REGION OF $flank_size bp TO EACH SIDE OF EACH FEATURE IN $input_gff:
   $cmd                    = "bedtools slop -i $input_gff -g $scaffold_length_file -b $flank_size > $output_gff";
   $systemcall             = system "$cmd";
   # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
   throw Error::Simple("ERRORCODE=1: add_flanking_region_to_gff_features: return value $systemcall when tried to run $cmd") if ($systemcall != 0); # TESTED FOR THIS
   system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: get_gff_lines_for_feature_types(): EXTRACT THE GFF LINES FOR PARTICULAR FEATURE TYPES FROM A GFF FILE, AND PUT IN A NEW GFF FILE. THIS DOES NOT INCLUDE FASTA SEQUENCE IN THE OUTPUT GFF.

sub get_gff_lines_for_feature_types
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $feature_types       = $_[1]; # ARRAY OF FEATURE TYPES TO TAKE
   my $output_gff          = $_[2]; # OUTPUT GFF FILE 
   my $new_gff;                     # NEW GFF FILE  
   my $line;                        # 
   my @temp;                        # 
   my %TAKE                = ();    # HASH OF FEATURE TYPES TO TAKE
   my $i;                           # 
   my $feature_type;                # A FEATURE TYPE
   my $end_of_gff          = 0;     # SAYS WHETHER WE HAVE REACHED THE END OF THE INPUT GFF FILE
   
   # RECORD THE FEATURE TYPES TO TAKE IN A HASH: 
   for ($i = 0; $i < @$feature_types; $i++)
   {
      $feature_type        = $feature_types->[$i];
      $TAKE{$feature_type} = 1;
   }

   # OPEN THE OUTPUT GFF:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: get_gff_lines_for_feature_types: cannot open output_gff $output_gff\n";

   # READ THROUGH THE INPUT GFF:
   open(INPUT_GFF,"$input_gff") || die "ERROR: get_gff_lines_for_feature_types: cannot open input_gff $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      if ($line eq '##FASTA') { $end_of_gff = 1;}
      if (substr($line,0,1) ne '#' && $end_of_gff == 0)
      {
         @temp             = split(/\t+/,$line);
         # THROW AN EXCEPTION IF THERE ARE NOT 9 COLUMNS IN THE GFF LINE:
         throw Error::Simple("ERRORCODE=1: get_gff_lines_for_feature_types: do not have 9 columns in $input_gff: line $line") if ($#temp != 8); # TESTED FOR
         $feature_type     = $temp[2];
         if ($TAKE{$feature_type}) { print OUTPUT_GFF "$line\n";}
      }
   }
   close(INPUT_GFF);
   close(OUTPUT_GFF); 
 
   return(1);
}
 
#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: read_gene_ids(): READ IN A FILE OF NEW VERSUS OLD GENE IDs, AND RETURN A HASH:

sub read_gene_ids
{
   my $new_gene_ids        = $_[0]; # FILE WITH NEW VERSUS OLD GENE IDS
   my %OLD2NEW             = ();    # HASH TABLE TO STORE NEW GENE IDS   
   my $line;                        # 
   my @temp;                        # 
   my $new;                         # NEW GENE ID
   my $old;                         # OLD GENE ID 
 
   # READ IN THE FILE OF OLD VERSUS NEW GENE IDS:
   # eg. 
   # SRAE_0036900 (old) = SRAE_0000000100 (new)
   # SRAE_0037000 (old) = SRAE_0000000200 (new)
   open(GENE_IDS,"$new_gene_ids") || die "ERROR: read_gene_ids: cannot open $new_gene_ids\n";
   while(<GENE_IDS>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\s+/,$line);
      $new                 = $temp[3];
      $old                 = $temp[0];
      # THROW AN EXCEPTION IF WE ALREADY HAVE A NEW GENE ID. FOR $old:
      throw Error::Simple("ERRORCODE=1: read_gene_ids: already have new gene id. for old gene id. $old") if (defined($OLD2NEW{$old})); # TESTED FOR
      $OLD2NEW{$old}       = $new;
   }
   close(GENE_IDS);

   return(\%OLD2NEW);   
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: replace_gene_ids_in_gff(): READ IN AN OLD GFF FILE, AND REPLACE THE GENE NAMES, WRITE OUT A NEW GFF FILE:

sub replace_gene_ids_in_gff
{
   my $input_gff           = $_[0]; # INPUT GFF FILE NAME
   my $output_gff          = $_[1]; # OUTPUT GFF FILE NAME
   my $OLD2NEW             = $_[2]; # HASH TABLE OF NEW GENE IDS.
   my $line;                        #  
   my @temp;                        # 
   my @temp2;                       # 
   my $feature;                     # GFF FEATURE
   my $gene;                        # GENE IN GFF LINE 
   my $new_gene;                    # NEW GENE NAME
   my $name;                        # LAST COLUMN OF GFF LINE 
   my $transcript;                  # TRANSCRIPT FOR GENE
   my $cds;                         # CDS FOR GENE 

   # OPEN THE OUTPUT GFF FILE:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: replace_gene_ids_in_gff: cannot open $output_gff\n"; 

   # READ IN THE INPUT GFF FILE:
   open(INPUT_GFF,"$input_gff") || die "ERROR: replace_gene_ids_in_gff: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $feature             = $temp[2];
      if ($feature eq 'gene')
      {
         $gene             = $temp[8];  # eg. ID=SRAE_2424300
         @temp2            = split(/ID=/,$gene);
         $gene             = $temp2[1]; # eg. SRAE_2424300
         # THROW AN EXCEPTION IF WE DO NOT HAVE THE NEW GENE ID. FOR $gene:
         throw Error::Simple("ERRORCODE=2: replace_gene_ids_in_gff: do not know new gene id. for $gene") if !(defined($OLD2NEW->{$gene})); # TESTED FOR
         $new_gene         = $OLD2NEW->{$gene};
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$new_gene\n";
      }
      elsif ($feature eq 'mRNA')
      {
         $name             = $temp[8]; # eg. ID=SRAE_2424300.t1:mRNA;Parent=SRAE_2424300
         @temp2            = split(/ID=/,$name); 
         $transcript       = $temp2[1]; # eg. SRAE_2424300.t1:mRNA;Parent=SRAE_2424300
         @temp2            = split(/\;/,$transcript);
         $transcript       = $temp2[0]; # eg. SRAE_2424300.t1:mRNA  
         @temp2            = split(/\./,$transcript);
         $gene             = $temp2[0]; # eg. SRAE_2424300
         $transcript       = $temp2[1]; # .t1:mRNA     
         # THROW AN EXCEPTION IF WE DO NOT HAVE THE NEW GENE ID. FOR $gene:
         throw Error::Simple("ERRORCODE=2: replace_gene_ids_in_gff: do not know new gene id. for $gene") if !(defined($OLD2NEW->{$gene})); # TESTED FOR
         $new_gene         = $OLD2NEW->{$gene};
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$new_gene.$transcript;Parent=$new_gene\n";   
      }
      elsif ($feature eq 'CDS')
      {
         $name             = $temp[8]; # ID=SRAE_2424300.t1:exon:1;Parent=SRAE_2424300.t1:mRNA
         @temp2            = split(/ID=/,$name);
         $cds              = $temp2[1]; # SRAE_2424300.t1:exon:1;Parent=SRAE_2424300.t1:mRNA
         @temp2            = split(/\;/,$cds);
         $cds              = $temp2[0]; # eg. SRAE_X125400.t1:exon:1
         @temp2            = split(/\./,$cds);
         $gene             = $temp2[0]; # eg. SRAE_X125400 
         $cds              = $temp2[1]; # eg. .t1:exon:1
         @temp2            = split(/Parent=/,$name);
         $transcript       = $temp2[1]; # eg. SRAE_2424300.t1:mRNA
         @temp2            = split(/\./,$transcript);
         $gene             = $temp2[0]; # eg. SRAE_2424300 
         $transcript       = $temp2[1]; # eg. t1:mRNA
         # THROW AN EXCEPTION IF WE DO NOT HAVE THE NEW GENE ID. FOR $gene:
         throw Error::Simple("ERRORCODE=2: replace_gene_ids_in_gff: do not know new gene id. for $gene") if !(defined($OLD2NEW->{$gene})); # TESTED FOR
         $new_gene         = $OLD2NEW->{$gene};
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$new_gene.$cds;Parent=$new_gene.$transcript\n";       
      }
      else 
      {
         # THROW AN EXCEPTION IF THE FEATURE IS NOT OF TYPE 'gene'/'mRNA'/'CDS':
         throw Error::Simple("ERRORCODE=3: replace_gene_ids_in_gff: feature $feature"); # TESTED FOR
      }
   }
   close(INPUT_GFF);
   close(OUTPUT_GFF);
  
   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: add_gene_features_to_gff_for_subsequences(): GIVEN A GFF OF GENES, AND A SECOND GFF OF THE GENES WITH RESPECT TO SUBSEQUENCES OF THE SCAFFOLDS, AND A HASH SAYING HOW MUCH THE START OF THE GENES DIFFER BETWEEN THE FIRST AND SECOND GFFs, MAKE AN OUTPUT GFF WITH ALL THE GENES' FEATURES (eg. EXONS, etc.) WITH RESPECT TO THE SUBSEQUENCES OF SCAFFOLDS:

sub add_gene_features_to_gff_for_subsequences
{
   my $input_gff           = $_[0]; # INPUT GFF FILE OF GENE FEATURES IN SCAFFOLDS
   my $gene_gff            = $_[1]; # SECOND GFF FILE, OF THE GENES WITH RESPECT TO SUBSEQUENCES OF THE SCAFFOLDS
   my $DIFF_IN_START       = $_[2]; # HASH TABLE SAYING HOW MUCH THE START OF THE GENES DIFFER BETWEEN THE FIRST AND SECOND GFFs
   my $output_gff          = $_[3]; # OUTPUT GFF FILE
   my $TRANSCRIPT2GENE     = $_[4]; # HASH TABLE OF GENE NAMES FOR TRANSCRIPT NAMES 
   my $line;                        # 
   my @temp;                        # 
   my $feature_type;                # FEATURE TYPE
   my $feature_name;                # FEATURE NAME 
   my $feature_start;               # FEATURE START POSITION
   my $feature_end;                 # FEATURE END POSITION
   my $gene;                        # GENE NAME
   my $diff_in_start;               # DIFFERENCE IN START OF A GENE BETWEEN THE FIRST AND SECOND GFFs 
   my $scaffold;                    # SCAFFOLD NAME
   my %SCAFFOLD            = ();    # HASH TABLE TO RECORD SCAFFOLD NAME FOR A GENE
   my $feature_length;              # LENGTH OF A FEATURE 

   # OPEN THE OUTPUT GFF FILE:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: add_gene_features_to_gff_for_subsequences: cannot open output_gff $output_gff\n";
  
   # PUT ALL THE LINES FROM $gene_gff IN $output_gff:
   open(GENE_GFF,"$gene_gff") || die "ERROR: add_gene_features_to_gff_for_subsequences: cannot open $gene_gff\n";
   while(<GENE_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $gene                = $temp[8];
      # THROW AN EXCEPTION IF WE ALREADY HAVE STORED A SCAFFOLD NAME FOR THIS GENE:
      throw Error::Simple("ERRORCODE=3: add_gene_features_to_gff_for_subsequences: already stored a scaffold name $SCAFFOLD{$gene} for gene $gene, scaffold $scaffold") if (defined($SCAFFOLD{$gene})); # TESTED FOR
      $SCAFFOLD{$gene}     = $scaffold; 
      print OUTPUT_GFF "$line\n";
   } 
   close(GENE_GFF);

   # OPEN THE INPUT GFF FILE:
   open(INPUT_GFF,"$input_gff") || die "ERROR: add_gene_features_to_gff_for_subsequences: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\s+/,$line);
      # eg.
      # SSTP.contig.00005.803316	source	CDS	8340	9119	.	+	.	T1.SSTP.contig.00005.803316_1_1;Parent=T1.SSTP.contig.00005.803316_1_mRNA
      # SSTP.contig.00005.803316	source	gene	8171	9119	.	+	.	T1.SSTP.contig.00005.803316_1
      # SSTP.contig.00005.803316	source	mRNA	8171	9119	.	+	.	T1.SSTP.contig.00005.803316_1_mRNA;Parent=T1.SSTP.contig.00005.803316_1
      $scaffold            = $temp[0];
      $feature_type        = $temp[2];    
      $feature_start       = $temp[3];
      $feature_end         = $temp[4];
      $feature_length      = $feature_end - $feature_start + 1;
      $feature_name        = $temp[8];
      if ($feature_type eq 'CDS' || $feature_type eq 'mRNA' || $feature_type eq 'gene')
      {
         if ($feature_type eq 'CDS' || $feature_type eq 'mRNA')
         {
            # GET THE GENE NAME:
            $gene          = HelminthGenomeAnalysis::AvrilGffUtils::get_gene_name($feature_type,$feature_name,$TRANSCRIPT2GENE,$scaffold);
            # THROW AN EXCEPTION IF WE DO NOT KNOW THE DIFFERENCE IN START FOR THIS GENE:
            throw Error::Simple("ERRORCODE=2: add_gene_features_to_gff_for_subsequences: do not know diff_in_start for $gene") if !(defined($DIFF_IN_START->{$gene})); # TESTED FOR
            $diff_in_start = $DIFF_IN_START->{$gene};
            $feature_start = $feature_start - $diff_in_start;
            $feature_end   = $feature_start + $feature_length - 1;
            # GET THE SCAFFOLD NAME:
            throw Error::Simple("ERRORCODE=4: add_gene_features_to_gff_for_subsequences: do not know scaffold for gene $gene") if !(defined($SCAFFOLD{$gene})); # TESTED FOR
            $scaffold      = $SCAFFOLD{$gene};
            print OUTPUT_GFF "$scaffold\t$temp[1]\t$temp[2]\t$feature_start\t$feature_end\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\n";
         }
      }
      else
      {
         # THROW AN EXCEPTION IF THE FEATURE IS NOT OF TYPE 'gene'/'mRNA'/'CDS':
         throw Error::Simple("ERRORCODE=1: add_gene_features_to_gff_for_subsequences: feature_type $feature_type"); # TESTED FOR
      }
   }
   close(INPUT_GFF);   
   close(OUTPUT_GFF);

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS get_gene_name(): GIVEN A FEATURE NAME IN A GFF FILE, RETURN THE GENE NAME:

sub get_gene_name
{
   my $feature_type        = $_[0]; # FEATURE TYPE
   my $feature_name        = $_[1]; # FEATURE NAME
   my $TRANSCRIPT2GENE     = $_[2]; # HASH TABLE OF TRANSCRIPT NAMES TO GENE NAMES 
   my $scaffold            = $_[3]; # THE SCAFFOLD NAME
   my @temp;                        # 
   my $transcript;                  # TRANSCRIPT NAME
   my $gene;                        # GENE NAME 
 
   if ($feature_type eq 'CDS')
   {
      # eg. T1.SSTP.contig.00005.803316_1_1;Parent=T1.SSTP.contig.00005.803316_1_mRNA
      # eg. ID=augustus_masked-scaff_000010-processed-gene-0.3-mRNA-1:cds;Parent=augustus_masked-scaff_000010-processed-gene-0.3-mRNA-1
      if ($feature_name =~ /ID/)
      {
         @temp             = split(/ID=/,$feature_name);
         $feature_name     = $temp[1]; 
      } 
      @temp                = split(/Parent=/,$feature_name);
      $transcript          = $temp[1]; # eg. T1.SSTP.contig.00005.803316_1_mRNA
      @temp                = split(/\;/,$transcript);
      $transcript          = $temp[0];
      $transcript          = $scaffold."=".$transcript;
      # THROW AN EXCEPTION IF THE GENE NAME IS NOT KNOWN FOR THIS TRANSCRIPT:
      throw Error::Simple("ERRORCODE=1: get_gene_name: gene name is not known for transcript $transcript") if !(defined($TRANSCRIPT2GENE->{$transcript})); # TESTED FOR
      $gene                = $TRANSCRIPT2GENE->{$transcript};
   } 
   elsif ($feature_type eq 'exon')
   {
      # eg. ID=CGOC_0001473701-mRNA-1:exon:1;Parent=CGOC_0001473701-mRNA-1
      if ($feature_name =~ /ID/)
      {
         @temp             = split(/ID=/,$feature_name);
         $feature_name     = $temp[1];
      } 
      @temp                = split(/Parent=/,$feature_name);
      $transcript          = $temp[1]; # eg. T1.SSTP.contig.00005.803316_1_mRNA
      @temp                = split(/\;/,$transcript);
      $transcript          = $temp[0];
      $transcript          = $scaffold."=".$transcript;
      # THROW AN EXCEPTION IF THE GENE NAME IS NOT KNOWN FOR THIS TRANSCRIPT:
      throw Error::Simple("ERRORCODE=2: get_gene_name: gene name is not known for transcript $transcript") if !(defined($TRANSCRIPT2GENE->{$transcript})); # TESTED FOR
      $gene                = $TRANSCRIPT2GENE->{$transcript};
   }
   elsif ($feature_type eq 'mRNA')
   {
      # eg. T1.SSTP.contig.00005.803316_1_mRNA;Parent=T1.SSTP.contig.00005.803316_1
      # eg. ID=augustus_masked-scaff_000010-processed-gene-0.3-mRNA-1;Parent=augustus_masked-scaff_000010-processed-gene-0.3;Name=augustus_masked-scaff_000010-processed-gene-0.3-mRNA-1;_AED=0.00;_eAED=0.00;_QI=0|-1|0|1|-1|1|1|0|392
      if ($feature_name =~ /ID/)
      {
         @temp             = split(/ID=/,$feature_name);
         $feature_name     = $temp[1]; 
      } 
      @temp                = split(/Parent=/,$feature_name);
      $gene                = $temp[1]; # eg. augustus_masked-scaff_000010-processed-gene-0.3;Name=augustus_masked-scaff_000010-processed-gene-0.3-mRNA-1;_AED=0.00;_eAED=0.00;_QI=0|-1|0|1|-1|1|1|0|392
      @temp                = split(/\;/,$gene);
      $gene                = $temp[0]; # eg. augustus_masked-scaff_000010-processed-gene-0.3
   }
   elsif ($feature_type eq 'three_prime_UTR' || $feature_type eq 'five_prime_UTR')
   {
      # eg. ID=CGOC_0000014201-mRNA-1:five_prime_utr;Parent=CGOC_0000014201-mRNA-1 
      if ($feature_name =~ /ID/)
      {
         @temp             = split(/ID=/,$feature_name);
         $feature_name     = $temp[1]; 
      } 
      @temp                = split(/Parent=/,$feature_name);
      $transcript          = $temp[1]; # eg. CGOC_0000014201-mRNA-1 
      @temp                = split(/\;/,$transcript);
      $transcript          = $temp[0];
      $transcript          = $scaffold."=".$transcript;
      # THROW AN EXCEPTION IF THE GENE NAME IS NOT KNOWN FOR THIS TRANSCRIPT:
      throw Error::Simple("ERRORCODE=3: get_gene_name: gene name is not known for transcript $transcript") if !(defined($TRANSCRIPT2GENE->{$transcript})); 
      $gene                = $TRANSCRIPT2GENE->{$transcript};
   }
   elsif ($feature_type eq 'gene')
   {
      # eg. ID=CGOC_0001473701;Name=CGOC_0001473701
      if ($feature_name =~ /ID/)
      {
         @temp             = split(/ID=/,$feature_name);
         $feature_name     = $temp[1];
      } 
      @temp                = split(/\;/,$feature_name);
      $gene                = $temp[0]; # eg. CGOC_0001473701 
   }
   else
   {
      throw Error::Simple("ERRORCODE=2: get_gene_name: feature_type $feature_type"); # TESTED FOR
   }

   return($gene);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS read_genenames_for_transcripts(): READS IN THE GENE NAME FOR EACH TRANSCRIPT IN AN INPUT GFF FILE, AND RETURNS IT IN A HASH TABLE (WITH KEYS AS 'scaffold=transcript'):

sub read_genenames_for_transcripts
{
   my $gff                 = $_[0]; # INPUT GFF FILE
   my %TRANSCRIPT2GENE     = ();    # HASH TABLE OF GENE NAMES FOR TRANSCRIPT NAMES
   my $line;                        # 
   my @temp;                        #   
   my $end_of_gff          = 0;     # SAYS WHETHER WE HAVE REACHED THE END OF THE GFF FILE
   my $scaffold;                    # NAME OF THE SCAFFOLD
   my $gene;                        # GENE NAME
   my $transcript;                  # TRANSCRIPT NAME

   open(GFF,"$gff") || die "ERROR: read_genenames_for_transcripts: cannot open gff $gff\n";
   while(<GFF>)
   { 
      $line                = $_;
      chomp $line;
      if ($line eq '##FASTA') { $end_of_gff = 1;}
      if (substr($line,0,1) ne '#' && $end_of_gff == 0)
      {
         @temp             = split(/\t+/,$line);
         # THROW AN EXCEPTION IF THE GFF LINE DOES NOT HAVE 9 COLUMNS:
         throw Error::Simple("ERRORCODE=1 read_genenames_for_transcripts: line $line of gff $gff does not have 9 columns") if ($#temp != 8); # TESTED FOR
         # IF THIS IS A LINE WITH A TRANSCRIPT ON IT:
         if ($temp[2] eq 'mRNA') 
         {
            # eg. scaff_000010    maker   mRNA    88793   89828   .       +       .       ID=maker-scaff_000010-augustus-gene-0.21-mRNA-1;Parent=maker-scaff_000010-augustus-gene-0.21;... 
            # eg. SSTP.contig.00005.803316	source	mRNA	8171	9119	.	+	.	T1.SSTP.contig.00005.803316_1_mRNA;Parent=T1.SSTP.contig.00005.803316_1
            $scaffold      = $temp[0];
            $gene          = $temp[8];
            # GET THE NAME OF THE TRANSCRIPT:
            $transcript    = $temp[8];
            if ($transcript =~ /ID=/)
            {
               @temp       = split(/ID=/,$transcript);
               $transcript = $temp[1]; # eg. maker-scaff_000010-augustus-gene-0.21-mRNA-1;...
            }
            @temp          = split(/\;/,$transcript);
            $transcript    = $temp[0]; # eg. maker-scaff_000010-augustus-gene-0.21-mRNA-1
            $transcript    = $scaffold."=".$transcript; 
            # GET THE NAME OF THE GENE:
            @temp          = split(/Parent=/,$gene);
            $gene          = $temp[1]; # eg. maker-scaff_000010-augustus-gene-0.21;...
            @temp          = split(/\;/,$gene);
            $gene          = $temp[0]; # eg. maker-scaff_000010-augustus-gene-0.21  
            # RECORD THE GENE NAME FOR THE mRNA:
            throw Error::Simple("ERRORCODE=2 read_genenames_for_transcripts: already have gene name for transcript $transcript") if (defined($TRANSCRIPT2GENE{$transcript})); # TESTED FOR
            $TRANSCRIPT2GENE{$transcript} = $gene;
         }
      }
   }

   return(\%TRANSCRIPT2GENE);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: make_gff_for_features_in_sequence(): GIVEN AN INPUT GFF OF FEATURES (eg. GENES) IN SUBSEQUENCES OF A SCAFFOLD, AND THE START POSITION OF THE SUBSEQUENCES IN THE SCAFFOLDS, MAKES AN OUTPUT GFF OF FEATURES (eg. GENES) IN THE SEQUENCES.

sub make_gff_for_features_in_sequences
{
   my $input_gff           = $_[0]; # INPUT GFF OF FEATURES IN SUBSEQUENCES OF SCAFFOLDS
   my $input_gff2          = $_[1]; # INPUT GFF OF SUBSEQUENCES IN SCAFFOLDS
   my $output_gff          = $_[2]; # OUTPUT GFF OF FEATURES IN SCAFFOLD SEQUENCES
   my $line;                        # 
   my @temp;                        # 
   my $feature_start_in_subsequence;# START OF THE FEATURE IN THE SUBSEQUENCE OF THE SCAFFOLD
   my $feature_end_in_subsequence;  # END OF THE FEATURE IN THE SUBSEQUENCE OF THE SCAFFOLD
   my $feature;                     # FEATURE NAME 
   my $subsequence_start_in_scaff;  # START OF A SUBSEQUENCE IN A SCAFFOLD
   my $subsequence;                 # NAME OF A SUBSEQUECE
   my %START               = ();    # START POSITIONS OF SUBSEQUENCES IN SCAFFOLDS
   my $scaffold;                    # SCAFFOLD NAME 
   my $feature_start;               # FEATURE START IN SCAFFOLD
   my $feature_end;                 # FEATURE END IN SCAFFOLD 
   my @temp2;                       #  

   # OPEN THE OUTPUT GFF FILE:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: make_gff_for_features_in_sequences: cannot open $output_gff\n";

   # READ IN THE GFF FILE OF SUBSEQUENCES IN SCAFFOLDS:
   open(INPUT_GFF2,"$input_gff2") || die "ERROR: make_gff_for_features_in_sequences: cannot open $input_gff2\n";
   while(<INPUT_GFF2>)
   {
       $line                = $_;
       chomp $line;
       @temp                = split(/\t+/,$line);
       $subsequence_start_in_scaff = $temp[3];
       $subsequence         = $temp[8]; # eg. scaffold1=10=20
       # THROW AN ERROR IF WE HAVE ALREADY STORED THE START POSITION OF $subsequence:
       throw Error::Simple("ERRORCODE=1: make_gff_for_features_in_sequences: already stored start for $subsequence") if (defined($START{$subsequence})); # TESTED FOR
       # THROW AN ERROR IF WE ALREADY HAVE STORED THE END POSITION OF $subsequence:
       $START{$subsequence} = $subsequence_start_in_scaff;
   }
   close(INPUT_GFF2);
  
   # READ IN THE GFF FILE OF FEATURES IN SUBSEQUENCES OF SCAFFOLDS:
   open(INPUT_GFF, "$input_gff") || die "ERROR: make_gff_for_features_in_sequences: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                 = $_;
      chomp $line; 
      @temp                 = split(/\t+/,$line);
      $subsequence          = $temp[0]; # eg. scaffold1=10=20
      $feature_start_in_subsequence = $temp[3];
      $feature_end_in_subsequence = $temp[4];
      $feature              = $temp[8];
      assert(defined($START{$subsequence}));
      $subsequence_start_in_scaff = $START{$subsequence};
      # FIND THE START AND END OF THE FEATURE WITH RESPECT TO THE SCAFFOLD:
      $feature_start        = $feature_start_in_subsequence + $subsequence_start_in_scaff - 1;
      $feature_end         = $feature_end_in_subsequence + $subsequence_start_in_scaff - 1; 
      # eg. feature_start_in_subsequence = 2, subsequence_start_in_scaff = 10 ---> 10+2-1 = 11
      @temp2               = split(/=/,$subsequence);
      $scaffold            = $temp2[0]; # eg. scaffold1
      print OUTPUT_GFF "$scaffold\t$temp[1]\t$temp[2]\t$feature_start\t$feature_end\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\n";
   }
   close(INPUT_GFF);
   close(OUTPUT_GFF);

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: sort_gff: SORTS A GFF FILE BY CHROMOSOME, THEN BY THE LAST COLUMN (eg. GENE NAME), THEN BY START POINT (4th COLUMN):

sub sort_gff
{
   my $gff                 = $_[0]; # INPUT GFF FILE
   my $sorted_gff          = $_[1]; # NAME OF SORTED GFF FILE
   my $gff_type            = $_[2]; # GFF TYPE: CAN BE 'exon' OR 'CDS' 
   my $cmd;                         # COMMAND TO RUN
   my $systemcall;                  # RETURN VALUE FROM SYSTEM CALL TO RUN $cmd 
   my $temp_gff;                    # TEMPORARY GFF FILE 
   my $temp_gff2;                   # TEMPORARY GFF FILE
   my $line;                        # 
   my @temp;                        #  

   # THROW AN ERROR IF $gff_type IS NOT 'exon' OR 'CDS':
   throw Error::Simple("ERRORCODE=1: sort_gff: gff_type is $gff_type (should be exon/CDS") if ($gff_type ne 'exon' && $gff_type ne 'CDS'); # TESTED FOR
 
   if ($gff_type eq 'CDS')
   {
      $cmd                 = "sort $gff -k1,1 -k9,9 -k4,4n > $sorted_gff";
      $systemcall          = system "$cmd";
      # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
      throw Error::Simple("ERRORCODE=2: sort_gff: return value $systemcall when tried to run $cmd") if ($systemcall != 0); 
      system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED
   }
   elsif ($gff_type eq 'exon')
   {
      # eg. ID=CGOC_0000249501-mRNA-1:exon:1;Parent=CGOC_0000249501-mRNA-1 IN LAST COLUMN.
      # BREAK LAST COLUMN INTO TWO COLUMNS:
      $temp_gff            = HelminthGenomeAnalysis::AvrilFileUtils::make_filename('.'); 
      $cmd                 = "cat $gff | tr '\;' '\t' > $temp_gff";
      $systemcall          = system "$cmd";
      # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
      throw Error::Simple("ERRORCODE=3: sort_gff: return value $systemcall when tried to run $cmd") if ($systemcall != 0); 
      system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED

      # SORT $temp_gff BY CHROMOSOME, THEN GENE (TRANSCRIPT) NAME, THEN START POINT: 
      $temp_gff2           = HelminthGenomeAnalysis::AvrilFileUtils::make_filename('.'); 
      $cmd                 = "sort $temp_gff -k1,1 -k10,10 -k4,4n > $temp_gff2";
      $systemcall          = system "$cmd";
      # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
      throw Error::Simple("ERRORCODE=4: sort_gff: return value $systemcall when tried to run $cmd") if ($systemcall != 0); 
      system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED

      # STICK THE LAST TWO COLUMNS OF $temp_gff2 TOGETHER IN $output:
      open(OUTPUT,">$sorted_gff") || die "ERROR: sort_gff: cannot open $sorted_gff\n";
      open(TEMP_GFF2,"$temp_gff2") || die "ERROR: sort_gff: cannot open $temp_gff2\n";
      while(<TEMP_GFF2>)
      {
         $line             = $_;
         chomp $line;
         @temp             = split(/\t+/,$line);
         print OUTPUT "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\;$temp[9]\n";
      }
      close(TEMP_GFF2);
      close(OUTPUT);

      # DELETE TEMPORARY FILES:
      system "rm -f $temp_gff";
      system "rm -f $temp_gff2"; 
   }

   return(1);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: merge_overlapping_cds: MERGE OVERLAPPING CDS (EXON) FEATURES, OR CDS FEATURES THAT HAVE A 0-bp INTRON BETWEEN THEM: 

sub merge_overlapping_cds
{
   my $sorted_cds_gff      = $_[0]; # SORTED CDS GFF FILE
   my $sorted_exon_gff     = $_[1]; # SORTED EXON GFF FILE 
   my $input_gff           = $_[2]; # THE INPUT GFF FILE
   my $output_gff          = $_[3]; # THE OUTPUT GFF FILE 
   my $outputdir           = $_[4]; # DIRECTORY TO WRITE OUTPUT FILES IN 
   my $scaffold_length_file= $_[5]; # FILE WITH LENGTHS OF SCAFFOLDS  
   my $line;                        # 
   my $line2;                       # GFF LINE FOR CDS 2 
   my @temp;                        # 
   my @temp2;                       # 
   my $scaffold;                    # SCAFFOLD NAME
   my $scaffold2;                   # SCAFFOLD 2
   my $cds;                         # CDS NAME
   my $cds2;                        # CDS 2 NAME
   my $cds_start;                   # CDS START 
   my $cds_start2;                  # CDS 2 START  
   my $cds_end;                     # CDS END
   my $cds_end2;                    # CDS 2 END
   my $transcript;                  # TRANSCRIPT NAME
   my $transcript2;                 # TRANSCRIPT 2 NAME 
   my $end_of_gff          = 0;     # SAYS WHETHER WE'VE REACHED THE END OF THE GFF FILE 
   my $feature_type;                # FEATURE TYPE IN GFF FILE 
   my $feature_name;                # FEATURE NAME IN GFF FILE 
   my %PRINTED             = ();    # HASH TABLE TO RECORD WHICH LINES WERE PRINTED TO THE OUTPUT GFF ALREADY  
   my $exon;                        # EXON NAME
   my $exon2;                       # EXON 2 NAME
   my $exon_start;                  # EXON START
   my $exon_start2;                 # EXON 2 START
   my $exon_end;                    # EXON END 
   my $exon_end2;                   # EXON 2 END
   my $cmd;                         # COMMAND TO RUN 
   my $overlaps;                    # FILE WITH OVERLAPS  
   my $systemcall;                  # RETURN VALUE FROM SYSTEM CALL 
   my %SEEN                = ();    # SAYS WHETHER WE HAVE ALREADY SEEN A PARTICULAR PAIR OF EXONS/CDSs 
   my $sorted_cds_gff2;             # COPY OF $sorted_cds_gff, WITH EACH CDS EXTENDED BY 1 bp IN EACH DIRECTION. 
   my $sorted_exon_gff2;            # COPY OF $sorted_exon_gff, WITH EACH EXON EXTENDED BY 1 IN EACH DIRECTION.
   my $returnvalue;                 # RETURN VALUE FROM A FUNCTION 
   my %REPLACE             = ();    # HASH TABLE TO RECORD WHICH NUMBER TO REPLACE A COORDINATE IN A GENE WITH.
   my $feature_start;               # START OF FEATURE
   my $feature_end;                 # END OF FEATURE

   # MAKE A COPY OF $sorted_cds_gff THAT HAS EACH CDS EXTENDED BY 1 bp IN EACH DIRECTION:
   $sorted_cds_gff2        = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::add_flanking_region_to_gff_features($sorted_cds_gff,1,$scaffold_length_file,$sorted_cds_gff2);

   # MAKE A COPY OF $sorted_exon_gff THAT HAS EACH EXON EXTENDED BY 1 bp IN EACH DIRECTION:
   $sorted_exon_gff2        = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::add_flanking_region_to_gff_features($sorted_exon_gff,1,$scaffold_length_file,$sorted_exon_gff2);

   # FIND CASES OF OVERLAPPING CDS IN THE $sorted_cds_gff FILE:
   $overlaps               = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $cmd                    = "bedtools intersect -loj -a $sorted_cds_gff2 -b $sorted_cds_gff2 > $overlaps"; 
   $systemcall             = system "$cmd";
   # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
   throw Error::Simple("ERRORCODE=2: merge_overlapping_cds: return value $systemcall when tried to run $cmd") if ($systemcall != 0); # TESTED FOR THIS
   system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED
   open(FILE,"$overlaps") || die "ERROR: merge_overlapping_cds: cannot open $overlaps\n";
   while(<FILE>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $cds_start           = $temp[3] + 1; # ADJUSTED TO BE LIKE IN $sorted_cds_gff
      $cds_end             = $temp[4] - 1; # ADJUSTED TO BE LIKE IN $sorted_cds_gff
      $cds                 = $temp[8];
      $cds                 = $scaffold."_".$cds_start."_".$cds_end."_".$cds;
      $scaffold2           = $temp[9];
      $cds_start2          = $temp[12] + 1; # ADJUSTED TO BE LIKE IN $sorted_cds_gff
      $cds_end2            = $temp[13] - 1; # ADJUSTED TO BE LIKE IN $sorted_cds_gff
      $cds2                = $temp[17];
      $cds2                = $scaffold2."_".$cds_start2."_".$cds_end2."_".$cds2;
      # DEFINE GFF LINE FOR CDS 1:
      $line                = $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$cds_start."\t".$cds_end."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8];
      # DEFINE GFF LINE FOR CDS 2:
      $line2               = $temp[9]."\t".$temp[10]."\t".$temp[11]."\t".$cds_start2."\t".$cds_end2."\t".$temp[14]."\t".$temp[15]."\t".$temp[16]."\t".$temp[17];
      # GET THE TRANSCRIPT NAME:
      @temp2               = split(/Parent=/,$cds);
      $transcript          = $temp2[1]; # eg. CGOC_0000249501-mRNA-1
      @temp2               = split(/\;/,$transcript);
      $transcript          = $temp2[0]; 
      # GET THE TRANSCRIPT 2 NAME:
      @temp2               = split(/Parent=/,$cds2);
      $transcript2         = $temp2[1];
      @temp2               = split(/\;/,$transcript);
      $transcript          = $temp2[0];
      # CHECK IF THE TWO TRANSCRIPTS ARE THE SAME:
      if ($transcript eq $transcript2 && $scaffold eq $scaffold2 && $cds ne $cds2)
      {    
         if ($cds_start2 <= $cds_start) 
         { 
            if (!(defined($REPLACE{$transcript."_CDSstart_".$cds_start}))) 
            { 
               $REPLACE{$transcript."_CDSstart_".$cds_start} = $cds_start2;
            }
            else
            {
               if ($cds_start2 < $REPLACE{$transcript."_CDSstart_".$cds_start})
               {
                  $REPLACE{$transcript."_CDSstart_".$cds_start} = $cds_start2;
               }
            } 
         } 
         if ($cds_end2 >= $cds_end)     
         { 
            if (!(defined($REPLACE{$transcript."_CDSend_".$cds_end})))
            {
               $REPLACE{$transcript."_CDSend_".$cds_end} = $cds_end2;
            }
            else
            {
               if ($cds_end2 > $REPLACE{$transcript."_CDSend_".$cds_end})
               {
                  $REPLACE{$transcript."_CDSend_".$cds_end} = $cds_end2;
               }
            }
         }
      }
   }
   close(FILE);
   system "rm -f $overlaps";

   # FIND CASES OF OVERLAPPING EXONS IN THE $sorted_exon_gff FILE:
   $overlaps               = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $cmd                    = "bedtools intersect -loj -a $sorted_exon_gff2 -b $sorted_exon_gff2 > $overlaps"; 
   $systemcall             = system "$cmd";
   # THROW AN ERROR IF THE RETURN VALUE FROM THE SYSTEM CALL IS NOT 0:
   throw Error::Simple("ERRORCODE=2: merge_overlapping_cds: return value $systemcall when tried to run $cmd") if ($systemcall != 0); # TESTED FOR THIS
   system "sleep 1"; # SLEEP FOR 1 SECOND TO MAKE SURE THE OUTPUT IS COMPLETED
   open(FILE,"$overlaps") || die "ERROR: merge_overlapping_exon: cannot open $overlaps\n";
   while(<FILE>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $exon_start          = $temp[3] + 1; # ADJUSTED TO BE LIKE IN $sorted_exon_gff
      $exon_end            = $temp[4] - 1; # ADJUSTED TO BE LIKE IN $sorted_exon_gff
      $exon                = $temp[8];
      $scaffold2           = $temp[9];
      $exon_start2         = $temp[12] + 1; # ADJUSTED TO BE LIKE IN $sorted_exon_gff
      $exon_end2           = $temp[13] - 1; # ADJUSTED TO BE LIKE IN $sorted_exon_gff
      $exon2               = $temp[17];
      # DEFINE GFF LINE FOR EXON 1:
      $line                = $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$exon_start."\t".$exon_end."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8];
      # DEFINE GFF LINE FOR EXON 2:
      $line2               = $temp[9]."\t".$temp[10]."\t".$temp[11]."\t".$exon_start2."\t".$exon_end2."\t".$temp[14]."\t".$temp[15]."\t".$temp[16]."\t".$temp[17];
      # GET THE TRANSCRIPT NAME:
      @temp2               = split(/Parent=/,$exon);
      $transcript          = $temp2[1]; # eg. CGOC_0000249501-mRNA-1
      @temp2               = split(/\;/,$transcript);
      $transcript          = $temp2[0]; 
      # GET THE TRANSCRIPT 2 NAME:
      @temp2               = split(/Parent=/,$exon2);
      $transcript2         = $temp2[1];
      @temp2               = split(/\;/,$transcript);
      $transcript          = $temp2[0];
      # CHECK IF THE TWO TRANSCRIPTS ARE THE SAME:
      if ($transcript eq $transcript2 && $scaffold eq $scaffold2 && $exon ne $exon2)
      {    
         if ($exon_start2 <= $exon_start) 
         { 
            if (!(defined($REPLACE{$transcript."_exonstart_".$exon_start}))) 
            { 
               $REPLACE{$transcript."_exonstart_".$exon_start} = $exon_start2;
            }
            else
            {
               if ($exon_start2 < $REPLACE{$transcript."_exonstart_".$exon_start})
               {
                  $REPLACE{$transcript."_exonstart_".$exon_start} = $exon_start2;
               }
            } 
         } 
         if ($exon_end2 >= $exon_end)     
         { 
            if (!(defined($REPLACE{$transcript."_exonend_".$exon_end})))
            {
               $REPLACE{$transcript."_exonend_".$exon_end} = $exon_end2;
            }
            else
            {
               if ($exon_end2 > $REPLACE{$transcript."_exonend_".$exon_end})
               {
                  $REPLACE{$transcript."_exonend_".$exon_end} = $exon_end2;
               }
            }
         }
      }
   }
   close(FILE);
   system "rm -f $overlaps";

   # CLOSE THE OUTPUT GFF FILE:
   open(OUTPUTGFF,">$output_gff") || die "ERROR: merge_overlapping_cds: cannot open $output_gff\n";

   # READ IN THE INPUT GFF AND REPLACE CDS LINES:
   open(INPUTGFF,"$input_gff") || die "ERROR: merge_overlapping_cds: cannot open $input_gff\n";
   while(<INPUTGFF>)
   {
      $line                = $_;
      chomp $line;
      if ($line eq '##FASTA') { $end_of_gff = 1;}
      if (substr($line,0,1) ne '#' && $end_of_gff == 0)
      {
         @temp             = split(/\t+/,$line);
         $feature_type     = $temp[2]; 
         $feature_start    = $temp[3];
         $feature_end      = $temp[4];
         $feature_name     = $temp[8];

         # THROW AN EXCEPTION IF THERE ARE NOT 9 COLUMNS IN THE GFF LINE:
         throw Error::Simple("ERRORCODE=1: merge_overlapping_cds: do not have 9 columns in $input_gff: line $line") if ($#temp != 8); # TESTED FOR 
         if ($feature_type eq 'CDS' || $feature_type eq 'exon')
         {
            # GET THE TRANSCRIPT NAME:
            @temp2         = split(/Parent=/,$feature_name);
            $transcript    = $temp2[1]; # eg. CGOC_0000249501-mRNA-1
            @temp2         = split(/\;/,$transcript);
            $transcript    = $temp2[0]; 
    
            # SEE IF WE WANT TO REPLACE THE START OR END POINT:
            if (defined($REPLACE{$transcript."_".$feature_type."start_".$feature_start})) { $feature_start = $REPLACE{$transcript."_".$feature_type."start_".$feature_start};}
            if (defined($REPLACE{$transcript."_".$feature_type."end_".$feature_end}))   { $feature_end   = $REPLACE{$transcript."_".$feature_type."end_".$feature_end};  } 

            # PRINT OUT THE LINE:
            $line          = $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$feature_start."\t".$feature_end."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8];

            if (!($PRINTED{$transcript."_".$feature_type."_".$feature_start."_".$feature_end}))
            {
               # PRINT OUT THE CDS LINE:
                print OUTPUTGFF "$line\n";
                $PRINTED{$transcript."_".$feature_type."_".$feature_start."_".$feature_end} = 1;
            }
         }
         else
         {
            print OUTPUTGFF "$line\n";
         }
      }
      else
      {
         print OUTPUTGFF "$line\n";
      }
   }
   close(INPUTGFF);
   close(OUTPUTGFF);

   # DELETE TEMPORARY FILES:
   system "rm -f $sorted_cds_gff2";
   system "rm -f $sorted_exon_gff2";

   return(1); 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: sort_gff_lines_for_genes(): ORDERS THE GFF LINES FOR EACH GENE IN ORDER 'gene', 'mRNA', 'exon', 'CDS':

sub sort_gff_lines_for_genes
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $output_gff          = $_[1]; # NAME OF THE OUTPUT GFF FILE
   my $line;                        # 
   my @temp;                        # 
   my $end_of_gff          = 0;     # SAYS WHETHER WE HAVE REACHED THE END OF THE GFF FILE
   my $feature_type;                # TYPE OF FEATURE 
   my $feature_name;                # FEATURE NAME
   my $scaffold;                    # SCAFFOLD NAME
   my $TRANSCRIPT2GENE;             # HASH TABLE OF GENE NAMES FOR TRANSCRIPT NAMES
   my %CDSLINES            = ();    # HASH TABLE TO RECORD CDS LINES FOR GENES
   my %MRNALINES           = ();    # HASH TABLE TO RECORD mRNA LINES FOR GENES
   my %EXONLINES           = ();    # HASH TABLE TO RECORD exon LINES FOR GENES   
   my %THREEPRIMEUTRLINES  = ();    # HASH TABLE TO RECORD three_prime_UTR LINES FOR GENES
   my %FIVEPRIMEUTRLINES   = ();    # HASH TABLE TO RECORD five_prime_UTR LINES FOR GENES 
   my $gene;                        # GENE NAME 

   # READ IN THE GENE NAMES FOR TRANSCRIPT NAMES:
   $TRANSCRIPT2GENE        = HelminthGenomeAnalysis::AvrilGffUtils::read_genenames_for_transcripts($input_gff);

   # READ IN THE INPUT GFF FILE, AND RECORD THE LINES FOR EACH GENE:
   open(INPUT_GFF,"$input_gff") || die "ERROR: sort_gff_lines_for_genes: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      if ($line eq '##FASTA') { $end_of_gff = 1;}
      if (substr($line,0,1) ne '#' && $end_of_gff == 0)
      {
         $scaffold         = $temp[0];
         $feature_type     = $temp[2]; 
         if ($feature_type eq 'CDS' || $feature_type eq 'mRNA' || $feature_type eq 'exon' || $feature_type eq 'three_prime_UTR' || $feature_type eq 'five_prime_UTR')
         {
            $feature_name  = $temp[8];
            # GET THE GENE NAME:
            $gene          = HelminthGenomeAnalysis::AvrilGffUtils::get_gene_name($feature_type,$feature_name,$TRANSCRIPT2GENE,$scaffold);
            # RECORD THE CDS/mRNA/exon LINES FOR $gene:
            if    ($feature_type eq 'CDS')
            {
               if (!(defined($CDSLINES{$gene}))) { $CDSLINES{$gene} = $line;}
               else { $CDSLINES{$gene} = $CDSLINES{$gene}."\n".$line;       }
            } 
            elsif ($feature_type eq 'exon')
            {
               if (!(defined($EXONLINES{$gene}))) { $EXONLINES{$gene} = $line;}
               else { $EXONLINES{$gene} = $EXONLINES{$gene}."\n".$line;       }
            }
            elsif ($feature_type eq 'mRNA')
            {
               if (!(defined($MRNALINES{$gene}))) { $MRNALINES{$gene} = $line;}
               else { $MRNALINES{$gene} = $MRNALINES{$gene}."\n".$line;       }
            }
            elsif ($feature_type eq 'three_prime_UTR')
            {
               if (!(defined($THREEPRIMEUTRLINES{$gene}))) { $THREEPRIMEUTRLINES{$gene} = $line;}
               else { $THREEPRIMEUTRLINES{$gene} = $THREEPRIMEUTRLINES{$gene}."\n".$line;       }
            }
            elsif ($feature_type eq 'five_prime_UTR')
            {
               if (!(defined($FIVEPRIMEUTRLINES{$gene}))) { $FIVEPRIMEUTRLINES{$gene} = $line;}
               else { $FIVEPRIMEUTRLINES{$gene} = $FIVEPRIMEUTRLINES{$gene}."\n".$line;       }
            }
         }
      }
   }
   close(INPUT_GFF);

   # OPEN THE OUTPUT GFF FILE:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: sort_gff_lines_for_genes: cannot open $output_gff\n";
   
   # NOW WRITE OUT TO THE OUTPUT GFF FILE:
   $end_of_gff             = 0; 
   open(INPUT_GFF,"$input_gff") || die "ERROR: sort_gff_lines_for_genes: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      if ($line eq '##FASTA') { $end_of_gff = 1;}
      if (substr($line,0,1) ne '#' && $end_of_gff == 0)
      {
         $scaffold         = $temp[0];
         $feature_type     = $temp[2]; 
         if ($feature_type eq 'gene')
         {
            $feature_name  = $temp[8];
            # GET THE GENE NAME:
            $gene          = HelminthGenomeAnalysis::AvrilGffUtils::get_gene_name($feature_type,$feature_name,$TRANSCRIPT2GENE,$scaffold);
            # PRINT THE LINE FOR THE GENE:
            print OUTPUT_GFF "$line\n";
            # THERE SHOULD BE mRNA LINES DEFINED FOR $gene, AS THEY WERE USED IN read_genenames_for_transcripts():
            assert(defined($MRNALINES{$gene})); 
            print OUTPUT_GFF "$MRNALINES{$gene}\n"; 
            # THROW AN ERROR IF exon LINES ARE NOT DEFINED FOR $gene:
            throw Error::Simple("ERRORCODE=2: sort_gff_lines_for_genes: no exon lines recorded for gene $gene") if (!(defined($EXONLINES{$gene}))); # TESTED FOR
            print OUTPUT_GFF "$EXONLINES{$gene}\n";
            # CHECK IF THERE ARE three_prime_UTR LINES FOR $gene:
            if (defined($THREEPRIMEUTRLINES{$gene})) { print OUTPUT_GFF "$THREEPRIMEUTRLINES{$gene}\n"; }  
            # CHECK IF THERE ARE five_prime_UTR LINES FOR $gene:
            if (defined($FIVEPRIMEUTRLINES{$gene})) { print OUTPUT_GFF "$FIVEPRIMEUTRLINES{$gene}\n"; }
            # THROW AN ERROR IF CDS LINES ARE NOT DEFINED FOR $gene:
            throw Error::Simple("ERRORCODE=3: sort_gff_lines_for_genes: no CDS lines recorded for gene $gene") if (!(defined($CDSLINES{$gene}))); # TESTED FOR
            print OUTPUT_GFF "$CDSLINES{$gene}\n"; 
         }
         else
         {
            # IGNORE 'mRNA'/'exon'/'CDS' LINES:
            if ($feature_type ne 'mRNA' && $feature_type ne 'CDS' && $feature_type ne 'exon' && $feature_type ne 'three_prime_UTR' && $feature_type ne 'five_prime_UTR')
            {
               print OUTPUT_GFF "$line\n";
            }
         }
      }
      else
      {
         print OUTPUT_GFF "$line\n";
      }
   }
   close(INPUT_GFF);

   return(1); 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: count_features_in_region: SUBROUTINE TO COUNT THE NUMBER OF FEATURES IN A GFF FILE BETWEEN $the_start AND $the_end ON SCAFFOLD $the_scaffold:

sub count_features_in_region
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $the_scaffold        = $_[1]; # THE SCAFFOLD OF INTEREST
   my $the_start           = $_[2]; # THE START OF THE REGION OF INTEREST
   my $the_end             = $_[3]; # THE END OF THE REGION OF INTEREST 
   my $line;                        # 
   my @temp;                        #  
   my $num_features        = 0;     # NUMBER OF FEATURES BETWEEN $the_start AND $the_end
   my $scaffold;                    # SCAFFOLD
   my $start;                       # FEATURE START
   my $end;                         # FEATURE END 

   open(INPUT_GFF,"$input_gff") || die "ERROR: count_features_in_region: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $start               = $temp[3];
      $end                 = $temp[4];
      if ($scaffold eq $the_scaffold && $start >= $the_start && $end >= $the_start && $start <= $the_end && $end <= $the_end)
      {
         $num_features++;
      }
   }
   close(INPUT_GFF);

   return($num_features);
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: convert_exonerate_gff_to_standard_gff: SUBROUTINE TO CONVERT EXONERATE GFF OUTPUT TO A MORE STANDARD GFF FORMAT.

sub convert_exonerate_gff_to_standard_gff
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $output_gff          = $_[1]; # OUTPUT GFF FILE
   my $line;                        # 
   my @temp;                        # 
   my $feature_type;                # FEATURE TYPE
   my $scaffold;                    # SCAFFOLD START
   my $feature_start;               # FEATURE START
   my $feature_end;                 # FEATURE END
   my $strand;                      # FEATURE STRAND 
   my $gene_num            = 0;     # GENE NUMBER IN THE EXONERATE FILE
   my $gene;                        # GENE NAME 
   my $exon_num            = 0;     # EXON NUMBER IN A GENE 

   # OPEN THE OUTPUT GFF:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: convert_exonerate_gff_to_standard_gff: cannot open $output_gff\n";

   # READ IN THE EXONERATE GFF:
   open(INPUT_GFF,"$input_gff") || die "ERROR: convert_exonerate_gff_to_standard_gff: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      # THROW AN EXCEPTION IF THERE ARE NOT 9 COLUMNS IN THE GFF LINE:
      throw Error::Simple("ERRORCODE=1: convert_exonerate_gff_to_standard_gff: do not have 9 columns in $input_gff: line $line") if ($#temp != 8); # TESTED FOR
      $scaffold            = $temp[0];
      $feature_type        = $temp[2];
      $feature_start       = $temp[3];
      $feature_end         = $temp[4];
      $strand              = $temp[6];
      if    ($feature_type eq 'gene')
      {
         # eg. scaffold_000039 exonerate:protein2genome:local  gene    252337  253077  685     +       .       gene_id 0 ; sequence SRAE_2000311000.t1:mRNA ; gene_orientation .
         $gene_num++;
         $exon_num         = 0;
         $gene             = "gene".$gene_num;
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene;Parent=$gene\n";
      } 
      elsif ($feature_type eq 'cds')
      {
         # scaffold_000039 exonerate:protein2genome:local  cds     252337  253077  .       +       .       cds
         print OUTPUT_GFF "$temp[0]\t$temp[1]\tCDS\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:cds;Parent=$gene\n";     
      }
      elsif ($feature_type eq 'exon')
      {
         # scaffold_000039 exonerate:protein2genome:local  exon    252337  253077  .       +       .       insertions 0 ; deletions 2
         $exon_num++;
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:exon:$exon_num;Parent=$gene\n";
      } 
      elsif ($feature_type eq 'similarity')
      {
         # scaffold_000039 exonerate:protein2genome:local  similarity      252337  253077  685     +       .       alignment_id 0 ; Query SRAE_2000311000.t1:mRNA ; Align 25175 81 564 ; Align 25739 271 177
         print OUTPUT_GFF "$line\n";
      }
      elsif ($feature_type eq 'splice5')
      {
         # scaffold_000043 exonerate:protein2genome:local  splice5 45991   45992   .       +       .       intron_id 1 ; splice_site "GT"
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:splice5;Parent=$gene\n";
      }
      elsif ($feature_type eq 'splice3')
      {
         # scaffold_000043 exonerate:protein2genome:local  splice3 46138   46139   .       +       .       intron_id 0 ; splice_site "AT"
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:splice3;Parent=$gene\n";
      }
      elsif ($feature_type eq 'intron')
      {
         # scaffold_000043 exonerate:protein2genome:local  intron  45991   46139   .       +       .       intron_id 1
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:intron;Parent=$gene\n";
      }
      else
      {
         throw Error::Simple("ERRORCODE=1: convert_exonerate_gff_to_standard_gff: unknown feature type $feature_type"); # TESTED FOR 
      }
   } 
   close(INPUT_GFF);
   close(OUTPUT_GFF);

   return(1);
}

#------------------------------------------------------------------#

1;
