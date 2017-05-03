#!/usr/bin/env perl

=head1 NAME

    run_exonerate_afterblast.pl

=head1 SYNOPSIS
 
    run_exonerate_afterblast.pl input_fasta input_pep output outputdir eval_cutoff flank_length blast_path type
        where input_fasta is the input fasta file of scaffolds in the assembly,
              input_pep is the input file of proteins or ESTs to be aligned
              output is the exonerate output file, 
              outputdir is the output directory for writing output files,
              eval_cutoff is the evalue cutoff to use for blast matches,
              flank_length is the length of DNA to take on either side of a blast match,
              blast_path is the path to the blast and makeblastdb executables,
              type (prot/est) says whether the queries in input_pep are protein or EST.

=head1 DESCRIPTION

    This script takes an input fasta file of scaffolds (<input_fasta>) and an input
    file of proteins (<input_pep>). For each protein and scaffold, we first check is there
    a blast match between the protein and scaffold. If so, we run exonerate on the scaffold,
    using that protein. The output is then written in the output file <output>.
    
=head1 VERSION
  
    Perl script last edited 15-Jul-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script run_exonerate_afterblast.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 15-Jul-13.
# Last edited 15-Jul-2013.
# SCRIPT SYNOPSIS: run_exonerate_afterblast.pl: run exonerate using input proteins, in scaffold regions that have blast matches to those proteins.
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

# xxx
# BEGIN {
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules');
# }

use HelminthGenomeAnalysis::AvrilFileUtils;
use HelminthGenomeAnalysis::AvrilFastaUtils;
use HelminthGenomeAnalysis::AvrilAlignUtils;
use HelminthGenomeAnalysis::AvrilGffUtils;

my $num_args               = $#ARGV + 1;
if ($num_args != 8)
{
    print "Usage of run_exonerate_afterblast.pl\n\n";
    print "perl run_exonerate_afterblast.pl <input_fasta> <input_pep> <output> <outputdir> <eval_cutoff> <flank_length> <blast_path> <type>\n";
    print "where <input_fasta> is the input fasta file of scaffolds in the assembly,\n";
    print "      <input_pep> is the input file of proteins or ESTs to be aligned,\n";
    print "      <output> is the genewise output file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <eval_cutoff> is the evalue cutoff to use for blast matches,\n";
    print "      <flank_length> is the length of DNA to take on either side of a blast match,\n";
    print "      <blast_path> is the path to the blast and makeblastdb executables,\n";
    print "      <type> (prot/EST) says whether the queries in input_pep are protein or EST\n";
    print "For example, >perl run_exonerate_afterblast.pl \n";
    print "assembly.fa proteins.fa exonerate.out . 0.05 25000 /software/pubseq/bin/ncbi_blast+/ prot\n";
    exit;
}

# FIND THE PATH TO THE INPUT FASTA FILE:                     

my $input_fasta            = $ARGV[0];

# FIND THE INPUT FILE OF PROTEINS:

my $input_pep              = $ARGV[1];

# FIND THE EXONERATE OUTPUT FILE:

my $output                 = $ARGV[2];

# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      

my $outputdir              = $ARGV[3];

# FIND THE EVALUE CUTOFF TO USE FOR BLAST MATCHES:

my $evalue_cutoff          = $ARGV[4];

# FIND THE LENGTH OF DNA TO TAKE ON EITHER SIDE OF A BLAST MATCH:

my $flank_length           = $ARGV[5];

# FIND THE PATH TO THE FORMATDB AND BLASTALL EXECUTABLES:

my $blast_path             = $ARGV[6];

# FOUND OUT WHETHER THE QUERY SEQUENCES IN $input_pep ARE PROTEIN/EST:

my $type                   = $ARGV[7];
if ($type ne 'est' && $type ne 'prot') { print STDERR "ERROR: type is $type but should be prot or EST\n"; exit;}

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($outputdir,$input_fasta,$input_pep,$output,$evalue_cutoff,$flank_length,$blast_path,$type);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_fasta         = $_[1]; # THE INPUT FASTA FILE
   my $input_pep           = $_[2]; # THE INPUT FILE OF PROTEINS
   my $output              = $_[3]; # THE OUTPUT FILE NAME     
   my $evalue_cutoff       = $_[4]; # EVALUE CUTOFF TO USE FOR BLAST MATCHES
   my $flank_length        = $_[5]; # LENGTH OF SEQUENCE TO TAKE ON EITHER SIDE OF A BLAST MATCH
   my $blast_path          = $_[6]; # PATH TO THE FORMATDB AND BLASTALL EXECUTABLES
   my $type                = $_[7]; # SAYS WHETHER THE QUERY SEQUENCES ARE EST/PROTEIN
   my $input_pep_fasta_obj;         # OBJECT FOR THE INPUT FASTA FILE $input_pep 
   my $input_pep_contigs2seq;       # HASH TABLE OF SEQUENCES IN $input_pep 
   my $input_fasta_obj;             # OBJECT FOR THE INPUT FASTA FILE $input_fasta
   my $input_fasta_contigs2seq;     # HASH TABLE OF SEQUENCES IN $input_fasta 
   my $protein;                     # AN INPUT PROTEIN IN $input_pep
   my $input_pep_seq;               # SEQUENCE OF $protein
   my $input_pep_file;              # FILE CONTAINING $input_pep_seq 
   my $returnvalue;                 # RETURNVALUE FROM FUNCTIONS
   my $blast_output;                # BLAST OUTPUT FILE FOR $protein 
   my $no_scaffolds_with_hits;      # NUMBER OF SCAFFOLDS WITH BLAST HITS TO $protein  
   my $scaffold;                    # SCAFFOLD WITH A BLAST MATCH
   my $scaffold_seq;                # SEQUENCE OF $scaffold
   my $scaffold_seq_file;           # FILE FOR STORING $scaffold_seq 
   my $hits;                        # BLAST HITS IN A SEQUENCE
   my @hits;                        # BLAST HITS IN A SEQUENCE
   my $i;                           # 
   my $hit;                         # BLAST HIT IN A SEQUENCE 
   my @temp;                        # 
   my $hit_start;                   # HIT START IN A SEQUENCE
   my $hit_end;                     # HIT END IN A SEQUENCE
   my $scaffold_subseq;             # A SUBSEQUENCE IN A SCAFFOLD, WITH A BLAST HIT
   my $scaffold_subseq_gff;         # GFF FILE OF THE POSITIONS OF SUBSEQUENCES WITH HITS IN THE SCAFFOLDS
   my $scaffold_exonerate_outputs;  # GFF FILE OF EXONERATE OUTPUTS TO $protein IN SCAFFOLDS, WITH COORDINATES WRT SUBSEQUENCES OF SCAFFOLDS 
   my $scaffold_exonerate_outputs2; # GFF FILE OF EXONERATE OUTPUTS IN A SCAFFOLD, WITH COORDINATES WRT SCAFFOLDS
   my $BLAST;                       # HASH TABLE OF BLAST HITS 
   my $exonerate_output;            # EXONERATE OUTPUT FOR ONE BLAST HIT REGION 
   my $temp_output;                 # TEMPORARY OUTPUT FILE

   # MAKE A TEMPORARY OUTPUT FILE:
   $temp_output            = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   open(TEMP_OUTPUT,">$temp_output") || die "ERROR: run_main_program: cannot open $temp_output\n"; 
   close(TEMP_OUTPUT);

   # READ IN THE SEQUENCES IN THE INPUT SCAFFOLD FILE:
   $input_fasta_obj        = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $input_fasta);
   $input_fasta_contigs2seq= $input_fasta_obj->_contigs2seq;

   # READ IN THE SEQUENCES IN THE INPUT FASTA FILE OF PROTEINS:
   $input_pep_fasta_obj    = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $input_pep);
   $input_pep_contigs2seq  = $input_pep_fasta_obj->_contigs2seq; 

   # FORMAT THE DATABASE OF SCAFFOLDS AS A BLAST DATABASE:
   $returnvalue            = HelminthGenomeAnalysis::AvrilAlignUtils::format_blastdb($blast_path,$input_fasta,'nucl');
   
   # FOR EACH OF THE INPUT PROTEIN SEQUENCES, FIND ITS BLAST MATCH IN THE SCAFFOLDS:
   foreach $protein (keys %{$input_pep_contigs2seq})
   {
      $input_pep_seq       = $input_pep_contigs2seq->{$protein}; 
      # WRITE THE PROTEIN SEQUENCE TO A FILE:   
      $input_pep_file      = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
      $returnvalue         = HelminthGenomeAnalysis::AvrilFastaUtils::print_seq_to_fasta($input_pep_file,$input_pep_seq,$protein,'no');
      # FIND THE BLAST MATCH OF $protein IN THE SCAFFOLDS:
      print STDERR "Checking for blast matches for $protein in the scaffolds...\n";
      $blast_output        = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);  
      if ($type eq 'prot')
      {
         $returnvalue      = HelminthGenomeAnalysis::AvrilAlignUtils::run_blast($blast_path,$blast_output,'tblastn',$input_fasta,$input_pep_file);
      }
      elsif ($type eq 'est')
      {
         $returnvalue      = HelminthGenomeAnalysis::AvrilAlignUtils::run_blast($blast_path,$blast_output,'blastn',$input_fasta,$input_pep_file);
      }
      # READ IN THE BLAST HITS:
      $BLAST               = HelminthGenomeAnalysis::AvrilAlignUtils::record_blast_hits($blast_output,$evalue_cutoff,$flank_length);
      # CHECK IF THERE ARE ANY SCAFFOLDS WITH HITS:
      $no_scaffolds_with_hits = keys %{$BLAST};
      if ($no_scaffolds_with_hits >= 1) 
      {
         # MAKE A GFF FILE OF THE SUBSEQUENCES WITH HITS TO $protein IN THE SCAFFOLDS: 
         $scaffold_subseq_gff = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);  
         open(SCAFFOLD_SUBSEQ_GFF,">$scaffold_subseq_gff") || die "ERROR: run_main_program: cannot open $scaffold_subseq_gff\n";
         # MAKE A GFF FILE OF EXONERATE OUTPUTS FOR $protein FOR ALL THE SCAFFOLDS:
         $scaffold_exonerate_outputs = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);  
         open(SCAFFOLD_EXONERATE,">$scaffold_exonerate_outputs") || die "ERROR: run_main_program: cannot open $scaffold_exonerate_outputs\n";
         close(SCAFFOLD_EXONERATE);
         # AS THERE ARE BLAST MATCHES BETWEEN THIS PROTEIN AND THE SCAFFOLDS, RUN EXONERATE FOR THE PROTEIN AGAINST THE REGIONS OF THOSE SCAFFOLDS:
         foreach $scaffold (keys %{$BLAST})
         { 
            $scaffold_seq  = $input_fasta_contigs2seq->{$scaffold};
            $hits          = $BLAST->{$scaffold};
            @hits          = split(/\,/,$hits);
            for ($i = 0; $i <= $#hits; $i++)
            { 
               $hit        = $hits[$i];
               @temp       = split(/=/,$hit);
               $hit_start  = $temp[0];
               $hit_end    = $temp[1];
               # $hit_start AND $hit_end SHOULD BE WITHIN $scaffold_seq:
               if ($hit_start < 1)                   { $hit_start = 1;                  }
               if ($hit_end > length($scaffold_seq)) { $hit_end = length($scaffold_seq);}
               $hit        = $scaffold."=".$hit_start."=".$hit_end;

               print SCAFFOLD_SUBSEQ_GFF "$scaffold\tsource\tsubseq\t$hit_start\t$hit_end\t.\t+\t.\t$hit\n";
               $scaffold_subseq = substr($scaffold_seq,$hit_start-1,$hit_end-$hit_start+1);
               # WRITE THE SCAFFOLD SUBSEQUENCE TO A FILE:
               $scaffold_seq_file = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);
               $returnvalue = HelminthGenomeAnalysis::AvrilFastaUtils::print_seq_to_fasta($scaffold_seq_file,$scaffold_subseq,$hit,'no');
               # RUN EXONERATE BETWEEN THE PROTEIN AND THE SCAFFOLD:
               print STDERR "___ running exonerate between $protein and $scaffold ($hit_start-$hit_end)...\n"; 
               $exonerate_output = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);  
               $returnvalue = HelminthGenomeAnalysis::AvrilAlignUtils::run_exonerate($exonerate_output,$scaffold_seq_file,$input_pep_file,$outputdir,$type);
               # CONCATENATE THE OUTPUT FILE TO $scaffold_exonerate_outputs: 
               system "cat $exonerate_output >> $scaffold_exonerate_outputs"; 
               # DELETE TEMPORARY FILE:
               system "rm -f $scaffold_seq_file";
               system "rm -f $exonerate_output";
            }
         } 
         close(SCAFFOLD_SUBSEQ_GFF);
         # THE $scaffold_exonerate_outputs FILE HAS EXONERATE COORDINATES WITH RESPECT TO SUBSEQUENCES (THE BLAST MATCH REGIONS) OF SCAFFOLDS. 
         # THE $scaffold_subseq_gff HAS THE COORDINATES OF THE SUBSEQUENCES WITH RESPECT TO THE SCAFFOLDS. 
         # CONVERT THE EXONERATE GFF COORDINATES SO THAT THEY ARE WITH RESPECT TO THE SCAFFOLDS:
         $scaffold_exonerate_outputs2 = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir);   
         $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::make_gff_for_features_in_sequences($scaffold_exonerate_outputs,$scaffold_subseq_gff,$scaffold_exonerate_outputs2);
         system "rm -f $scaffold_exonerate_outputs"; 
         system "rm -f $scaffold_subseq_gff";
         # CONCATENATE $scaffold_exonerate_outputs2 TO $temp_output:
         system "cat $scaffold_exonerate_outputs2 >> $temp_output"; 
         system "rm -f $scaffold_exonerate_outputs2";
      } 
      # DELETE TEMPORARY FILES:
      system "rm -f $input_pep_file";
      system "rm -f $blast_output";
   }
  
   # CONVERT THE $temp_output FILE TO MORE STANDARD GFF FORMAT:
   $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::convert_exonerate_gff_to_standard_gff($temp_output,$output);
   system "rm -f $temp_output";
}

#------------------------------------------------------------------#

