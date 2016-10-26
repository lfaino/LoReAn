#!/usr/local/bin/perl

=head1 NAME

    gff_to_genbank.pl

=head1 SYNOPSIS
 
    gff_to_genbank.pl input_gff assembly outputdir 
        where input_gff is the input gff file,
              assembly is the fasta file for the assembly,
              outputdir is the output directory for writing output files.

=head1 DESCRIPTION

    This script takes an input gff file (<input_gff>), and converts it to 
    genbank format, and writes the output genbank files for each scaffold
    (called scaffold.gb, where scaffold is the scaffold name) in directory
    <outputdir>.

=head1 VERSION
  
    Perl script last edited 4-Apr-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script gff_to_genbank.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 4-Apr-13.
# Last edited 4-Apr-2013.
# SCRIPT SYNOPSIS: gff_to_genbank.pl: convert a gff file to genbank files for the scaffolds.
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 3)
{
    print "Usage of gff_to_genbank.pl\n\n";
    print "perl gff_to_genbank.pl <input_gff> <assembly> <outputdir>\n";
    print "where <input_gff> is the input gff file,\n";
    print "      <assembly> is the fasta file for the assembly,\n";
    print "      <outputdir> is the output directory for writing output files.\n";
    print "For example, >perl gff_to_genbank.pl my.gff PTRK.v1.fa\n";
    print "/lustre/scratch108/parasites/alc/StrongyloidesAugustus/PTRK100blast\n";
    exit;
}

# FIND THE PATH TO THE INPUT GFF FILE:                     

my $input_gff              = $ARGV[0];

# FIND THE FASTA FILE FOR THE ASSEMBLY:

my $assembly               = $ARGV[1];

# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      

my $outputdir              = $ARGV[2];

#------------------------------------------------------------------#

# TEST SUBROUTINES: 

my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_convert_gff_to_genbank($outputdir);
&test_count_bases;
&test_read_gene_positions($outputdir); 
&test_sort_input_gff($outputdir);
&test_read_scaffold_lengths($outputdir);  
&test_read_assembly($outputdir);
&test_print_error;
print STDERR "Finished tests, now running main code...\n";
 
#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($outputdir,$input_gff,$assembly);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_gff           = $_[1]; # THE INPUT GFF FILE  
   my $assembly            = $_[2]; # THE FASTA FILE FOR THE ASSEMBLY
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $LEN;                         # HASH TABLE TO STORE THE LENGTH OF SCAFFOLDS
   my $SEQ;                         # HASH TABLE TO STORE SCAFFOLD SEQUENCES 
   my $GENES;                       # HASH TABLE TO STORE THE GENES ON SCAFFOLDS
   my $EXONS;                       # HASH TABLE TO STORE THE EXONS IN GENES 
   my $sorted_input_gff;            # SORTED VERSION OF THE INPUT GFF FILE 

   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # SORT THE GFF FILE SO THAT WE HAVE gene,mRNA,CDS FOR EACH GENE:
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 

   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($sorted_input_gff); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:
   ($errorcode,$errormsg)  = &convert_gff_to_genbank($outputdir,$GENES,$EXONS,$LEN,$SEQ,0);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # DELETE TEMPORARY FILES:
   system "rm -f $sorted_input_gff";

}

#------------------------------------------------------------------#

# TEST &sort_input_gff

sub test_sort_input_gff
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $input_gff;                   # INPUT GFF FILE
   my $sorted_input_gff;            # SORTED VERSION OF $input_gff 
   my $expected_sorted_input_gff;   # FILE WITH THE EXPECTED CONTENTS OF $sorted_input_gff
   my $differences;                 # DIFFERENCES BETWEEN $sorted_input_gff AND $expected_sorted_input_gff
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_sort_input_gff: cannot open $input_gff\n"; 
   print INPUT_GFF "PTRK.scaffold.00008.1453643  source	CDS	115993	116307	.	+	.	T1.PTRK.scaffold.00008.1453643_1_0;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	CDS	116360	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_1;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	gene	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	mRNA	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	CDS	535391	535554	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_0;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	CDS	535599	537591	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_1;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	CDS	537642	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_2;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	gene	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1\n";
   print INPUT_GFF "PTRK.scaffold.00008.1453643	source	mRNA	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_mRNA;Parent=T2.PTRK.scaffold.00008.1453643.embl_1\n";
   close(INPUT_GFF);
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 0) { print STDERR "ERROR: test_sort_input_gff: failed test1\n"; exit;}
   ($expected_sorted_input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_sorted_input_gff") || die "ERROR: test_sort_input_gff: cannot open $expected_sorted_input_gff\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	gene	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	mRNA	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	CDS	115993	116307	.	+	.	T1.PTRK.scaffold.00008.1453643_1_0;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	CDS	116360	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_1;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	gene	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	mRNA	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_mRNA;Parent=T2.PTRK.scaffold.00008.1453643.embl_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	CDS	535391	535554	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_0;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	CDS	535599	537591	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_1;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643	source	CDS	537642	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_2;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $sorted_input_gff $expected_sorted_input_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_sort_input_gff: failed test1 (sorted_input_gff $sorted_input_gff expected_sorted_input_gff $expected_sorted_input_gff)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $sorted_input_gff";
   system "rm -f $expected_sorted_input_gff"; 

   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_sort_input_gff: cannot open $input_gff\n"; 
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(INPUT_GFF);
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 0) { print STDERR "ERROR: test_sort_input_gff: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_sorted_input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_sorted_input_gff") || die "ERROR: test_sort_input_gff: cannot open $expected_sorted_input_gff\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print EXPECTED "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $sorted_input_gff $expected_sorted_input_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_sort_input_gff: failed test3 (sorted_input_gff $sorted_input_gff expected_sorted_input_gff $expected_sorted_input_gff)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $sorted_input_gff";
   system "rm -f $expected_sorted_input_gff"; 

   # TEST ERRORCODE=17 (mRNA LINE APPEARS TWICE IN GFF):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_sort_input_gff: cannot open $input_gff\n"; 
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(INPUT_GFF);
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 17) { print STDERR "ERROR: test_sort_input_gff: failed test3\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $sorted_input_gff";

   # TEST ERRORCODE=19 (gene LINE APPEARS MORE THAN ONCE IN GFF):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_sort_input_gff: cannot open $input_gff\n"; 
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(INPUT_GFF);
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 19) { print STDERR "ERROR: test_sort_input_gff: failed test4\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $sorted_input_gff";

   # TEST ERRORCODE=20 (NO CDS LINES FOR A GENE):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_sort_input_gff: cannot open $input_gff\n"; 
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print INPUT_GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   close(INPUT_GFF);
   ($sorted_input_gff,$errorcode,$errormsg) = &sort_input_gff($input_gff);
   if ($errorcode != 20) { print STDERR "ERROR: test_sort_input_gff: failed test5\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $sorted_input_gff";

}

#------------------------------------------------------------------#

# SORT THE INPUT GFF FILE, SO THAT FOR EACH GENE WE HAVE gene FEATURES,
# THEN mRNA, THEN CDS:

sub sort_input_gff
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = 'none';# RETURNED AS 'none' IF THERE IS NO ERROR
   my $sorted_input_gff;            # SORTED VERSION OF THE INPUT GFF FILE
   my $line;                        # 
   my @temp;                        # 
   my $scaffold;                    # SCAFFOLD NAME
   my $feature;                     # GFF FEATURE TYPE 
   my $name;                        # NAME OF THE FEATURE
   my $mRNA;                        # NAME OF THE mRNA
   my $gene;                        # NAME OF THE GENE 
   my %GENE                = ();    # HASH TABLE TO STORE GENES FOR mRNAs 
   my %GENES               = ();    # HASH TABLE TO STORE GENES ON EACH SCAFFOLD 
   my %CDSLINES            = ();    # HASH TABLE TO RECORD THE CDS LINES FOR A GENE
   my %MRNALINES           = ();    # HASH TABLE TO RECORD THE mRNA LINES FOR A GENE
   my %GENELINES           = ();    # HASH TABLE TO RECORD THE GENE LINES FOR A GENE
   my $genes;                       # LIST OF GENES ON A SCAFFOLD
   my @genes;                       # LIST OF GENES ON A SCAFFOLD
   my $i;                           # 
   my $cdslines;                    # CDS LINES FOR A GENE
   my $mrnalines;                   # mRNA LINES FOR A GENE
   my $genelines;                   # GENE LINES FOR A GENE

   # OPEN THE OUTPUT SORTED GFF:
   ($sorted_input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(SORTED_INPUT_GFF,">$sorted_input_gff") || die "ERROR: sort_input_gff: cannot open $sorted_input_gff\n";

   # READ IN THE INPUT GFF FILE TO RECORD GENES NAME FOR EACH mRNA:
   open(INPUT_GFF,"$input_gff") || die "ERROR: sort_input_gff: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_; 
      chomp $line;
      @temp                = split(/\t+/,$line); 
      $scaffold            = $temp[0];
      $feature             = $temp[2];
      $name                = $temp[8]; 
      if ($feature eq 'mRNA')
      {
         # eg. PTRK.scaffold.00008.1453643	source	mRNA	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1
         @temp             = split(/\;/,$name);
         $mRNA             = $temp[0]; # eg. T1.PTRK.scaffold.00008.1453643_1_mRNA 
         if ($mRNA =~ /ID=/)
         {
            @temp          = split(/ID=/,$mRNA);
            $mRNA          = $temp[1];
         }
         @temp             = split(/Parent=/,$name);
         $gene             = $temp[1]; # eg. T1.PTRK.scaffold.00008.1453643_1
         # RECORD THE GENE NAME FOR THIS mRNA
         if ($GENE{$mRNA}) 
         {
            $errormsg      = "ERROR: sort_input_gff: already know gene for $mRNA\n";
            $errorcode     = 17; # ERRORCODE=17 (TESTED FOR)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg);
         }
         $GENE{$mRNA}      = $gene;
      }
      elsif ($feature eq 'gene')
      {
         # eg. PTRK.scaffold.00008.1453643	source	gene	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1
         $gene             = $name; # eg. T1.PTRK.scaffold.00008.1453643_1 
         if ($gene =~ /ID=/)
         {
            @temp          = split(/ID=/,$gene);
            $gene          = $temp[1];
         }
         if (!($GENES{$scaffold})) { $GENES{$scaffold} = $gene;}
         else {$GENES{$scaffold} = $GENES{$scaffold}.",".$gene;}
      }
   }
   close(INPUT_GFF); 

   # READ IN THE INPUT GFF FILE TO RECORD THE LINES FOR EACH GENE: 
   open(INPUT_GFF,"$input_gff") || die "ERROR: sort_input_gff: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_; 
      chomp $line;
      @temp                = split(/\t+/,$line); 
      $feature             = $temp[2];
      $name                = $temp[8]; 
      if ($feature eq 'CDS')
      {
         # eg. PTRK.scaffold.00008.1453643	source	CDS	115993	116307	.	+	.	T1.PTRK.scaffold.00008.1453643_1_0;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA
         @temp             = split(/Parent=/,$name);
         $mRNA             = $temp[1]; # eg. T1.PTRK.scaffold.00008.1453643_1_mRNA
         # FIND THE GENE FOR THIS mRNA:
         if (!($GENE{$mRNA}))
         {
            $errormsg      = "ERROR: sort_input_gff: do not know gene for mRNA $mRNA\n";
            $errorcode     = 18; # ERRORCODE=18 (SHOULDN'T HAPPEN, SO CAN'T TEST)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg);
         }
         $gene             = $GENE{$mRNA};
         # RECORD THE CDS FEATURES FOR THIS GENE: 
         if (!($CDSLINES{$gene})) { $CDSLINES{$gene} = $line; }
         else {$CDSLINES{$gene} = $CDSLINES{$gene}."\n".$line;}
      }
      elsif ($feature eq 'mRNA')
      {
         # eg. PTRK.scaffold.00008.1453643	source	mRNA	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1
         @temp             = split(/Parent=/,$name);
         $gene             = $temp[1]; # eg. T1.PTRK.scaffold.00008.1453643_1
         # RECORD THE mRNA FEATURES FOR THIS GENE:
         if (!($MRNALINES{$gene})) { $MRNALINES{$gene} = $line; }
         else {$MRNALINES{$gene} = $MRNALINES{$gene}."\n".$line;} 
      }
      elsif ($feature eq 'gene')
      {
         # eg. PTRK.scaffold.00008.1453643	source	gene	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1
         $gene             = $name;
         if ($gene =~ /ID=/)
         {
            @temp          = split(/ID=/,$gene);
            $gene          = $temp[1];
         }
         # RECORD THE GENE LINE FOR THIS GENE:
         if ($GENELINES{$gene})
         {
            $errormsg      = "ERROR: sort_input_gff: gene $gene appears more than once in the input gff\n";
            $errorcode     = 19; # ERRORCODE=19 (TESTED FOR)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg);
         }
         $GENELINES{$gene} = $line;
      }
   }
   close(INPUT_GFF); 

   # PRINT OUT THE GFF LINES FOR EACH GENE:
   foreach $scaffold (sort keys %GENES) # SORTING THE KEYS MAKES IT EASIER TO TEST
   {
      $genes               = $GENES{$scaffold};
      @genes               = split(/\,/,$genes);
      for ($i = 0; $i <= $#genes; $i++)
      {
         $gene             = $genes[$i];
         if (!($CDSLINES{$gene}))
         {
            $errormsg      = "ERROR: sort_input_gff: do not find any CDS lines for gene $gene\n";
            $errorcode     = 20; # ERRORCODE=20 (TESTED FOR)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg); 
         }   
         $cdslines         = $CDSLINES{$gene};
         if (!($MRNALINES{$gene}))
         {
            $errormsg      = "ERROR: sort_input_gff: do not find any mRNA lines for gene $gene\n";
            $errorcode     = 21; # ERRORCODE=21 (SHOULDN'T HAPPEN, ACTUALLY LEADS TO ERRORCODE=18 IF YOU DELETE mRNA LINES FOR A GENE)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg);
         }
         $mrnalines        = $MRNALINES{$gene};
         if (!($GENELINES{$gene}))
         {
            $errormsg      = "ERROR: sort_input_gff: do not find any gene lines for gene $gene\n";
            $errorcode     = 22; # ERRORCODE=22 (SHOULDN'T HAPPEN, SO CAN'T TEST FOR)
            system "rm -f $sorted_input_gff";
            return($sorted_input_gff,$errorcode,$errormsg);
         }
         $genelines        = $GENELINES{$gene};
         # PRINT OUT THE gene LINE, THEN mRNA LINES, then CDS LINES:
         print SORTED_INPUT_GFF "$genelines\n"; 
         print SORTED_INPUT_GFF "$mrnalines\n";
         print SORTED_INPUT_GFF "$cdslines\n";
      }
   }
   close(SORTED_INPUT_GFF);

   return($sorted_input_gff,$errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &count_bases

sub test_count_bases
{
   my $seq;                         # SEQUENCE
   my $num_as;                      # NUMBER OF As IN THE SEQUENCE
   my $num_cs;                      # NUMBER OF Cs IN THE SEQUENCE
   my $num_gs;                      # NUMBER OF Gs IN THE SEQUENCE
   my $num_ts;                      # NUMBER OF Ts IN THE SEQUENCE
   my $num_other;                   # NUMBER OF OTHER LETTERS IN THE SEQUENCE
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
 
   $seq                    = "AAACCGGTYYAACCGT";
   ($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg) = &count_bases($seq); 
   if ($errorcode != 0 || $num_as != 5 || $num_cs != 4 || $num_gs != 3 || $num_ts != 2 || $num_other != 2) { print STDERR "ERROR: test_count_bases: failed test1\n"; exit;}
   
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: count_bases(): COUNTS THE NUMBER OF As, Cs, Gs AND Ts IN A SEQUENCE

sub count_bases
{
   my $seq                 = $_[0]; # SEQUENCE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $num_as              = 0;     # NUMBER OF As IN THE SEQUENCE
   my $num_cs              = 0;     # NUMBER OF Cs IN THE SEQUENCE
   my $num_gs              = 0;     # NUMBER OF Gs IN THE SEQUENCE
   my $num_ts              = 0;     # NUMBER OF Ts IN THE SEQUENCE
   my $num_other           = 0;     # NUMBER OF OTHER LETTERS IN THE SEQUENCE
   my $length;                      # LENGTH OF THE SEQUENC
   my $i;                           # 
   my $base;                        # A BASE IN THE SEQUENCE

   $length                 = length($seq);
   $seq                    =~ tr/[a-z]/[A-Z]/;
   for ($i = 1; $i <= $length; $i++)
   {
      $base                = substr($seq,$i-1,1);
      if    ($base eq 'A') { $num_as++;   }
      elsif ($base eq 'C') { $num_cs++;   }
      elsif ($base eq 'G') { $num_gs++;   }
      elsif ($base eq 'T') { $num_ts++;   }
      else                 { $num_other++;} 
   }

   return($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &convert_gff_to_genbank

sub test_convert_gff_to_genbank
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $input_gff;                   # INPUT GFF FILE
   my $output_genbank;              # OUTPUT GENBANK FILE 
   my $expected_output_genbank;     # FILE WITH EXPECTED CONTENTS OF $output_genbank 
   my $differences;                 # DIFFERENCE BETWEEN $output_genbank AND $expected_output_genbank
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #
   my $assembly;                    # ASSEMBLY FASTA FILE
   my $LEN;                         # HASH TABLE OF LENGTHS OF SCAFFOLDS
   my $SEQ;                         # HASH TABLE OF SEQUENCES OF SCAFFOLDS 
   my $GENES;                       # HASH TABLE OF GENES ON SCAFFOLDS
   my $EXONS;                       # HASH TABLE OF EXONS IN GENES 
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_convert_gff_to_genbank: cannot open $input_gff\n";
   print INPUT_GFF "HS04636\tsource\tgene\t966\t2995\t.\t+\t.\tmygene\n";
   print INPUT_GFF "HS04636\tsource\tmRNA\t966\t2995\t.\t+\t.\tID=mymrna;Parent=mygene\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t966\t1017\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t1818\t1934\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2055\t2198\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2852\t2995\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   close(INPUT_GFF);
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_convert_gff_to_genbank: cannot open $assembly\n";
   print ASSEMBLY ">HS04636\n";
   print ASSEMBLY "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n";
   close(ASSEMBLY);
   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($input_gff); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:
   ($errorcode,$errormsg)  = &convert_gff_to_genbank($outputdir,$GENES,$EXONS,$LEN,$SEQ,1);
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_gff_to_genbank: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_output_genbank,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_genbank") || die "ERROR: test_convert_gff_to_genbank: $expected_output_genbank\n";
   print EXPECTED "LOCUS       HS04636   80 bp  DNA\n";
   print EXPECTED "FEATURES             Location/Qualifiers\n";
   print EXPECTED "     source          1..80\n";
   print EXPECTED "     CDS             join(966..1017,1818..1934,2055..2198,2852..2995)\n";
   print EXPECTED "BASE COUNT     20 a   20 c  20 g  20 t\n";
   print EXPECTED "ORIGIN\n";
   print EXPECTED "        1 aaaaaaaaaa cccccccccc gggggggggg tttttttttt aaaaaaaaaa cccccccccc\n";
   print EXPECTED "       61 gggggggggg tttttttttt\n";
   print EXPECTED "//\n";
   close(EXPECTED);
   $differences            = "";
   $output_genbank         = $outputdir."/HS04636.gb";
   open(TEMP,"diff $output_genbank $expected_output_genbank |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_gff_to_genbank: failed test1 (output_genbank $output_genbank expected_output_genbank $expected_output_genbank)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $output_genbank";
   system "rm -f $expected_output_genbank"; 
   system "rm -f $assembly";

   # TEST ERRORCODE=3 (DO NOT KNOW EXONS IN GENE):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_convert_gff_to_genbank: cannot open $input_gff\n";
   print INPUT_GFF "HS04636\tsource\tgene\t966\t2995\t.\t+\t.\tmygene\n";
   print INPUT_GFF "HS04636\tsource\tmRNA\t966\t2995\t.\t+\t.\tID=mymrna;Parent=mygene\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t966\t1017\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t1818\t1934\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2055\t2198\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2852\t2995\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   close(INPUT_GFF);
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_convert_gff_to_genbank: cannot open $assembly\n";
   print ASSEMBLY ">HS04636\n";
   print ASSEMBLY "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n";
   close(ASSEMBLY);
   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($input_gff); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   delete($EXONS->{'mygene'});
   # READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:
   ($errorcode,$errormsg)  = &convert_gff_to_genbank($outputdir,$GENES,$EXONS,$LEN,$SEQ,1);
   if ($errorcode != 3) { print STDERR "ERROR: test_convert_gff_to_genbank: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $output_genbank";
   system "rm -f $assembly";

   # TEST ERRORCODE=7 (DO NOT KNOW LENGTH OF SCAFFOLD):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_convert_gff_to_genbank: cannot open $input_gff\n";
   print INPUT_GFF "HS04636\tsource\tgene\t966\t2995\t.\t+\t.\tmygene\n";
   print INPUT_GFF "HS04636\tsource\tmRNA\t966\t2995\t.\t+\t.\tID=mymrna;Parent=mygene\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t966\t1017\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t1818\t1934\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2055\t2198\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2852\t2995\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   close(INPUT_GFF);
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_convert_gff_to_genbank: cannot open $assembly\n";
   print ASSEMBLY ">HS04636\n";
   print ASSEMBLY "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n";
   close(ASSEMBLY);
   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($input_gff); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:
   delete($LEN->{"HS04636"});
   ($errorcode,$errormsg)  = &convert_gff_to_genbank($outputdir,$GENES,$EXONS,$LEN,$SEQ,1);
   if ($errorcode != 7) { print STDERR "ERROR: test_convert_gff_to_genbank: failed test3 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $output_genbank";
   system "rm -f $assembly";

   # TEST ERRORCODE=9 (DO NOT KNOW SEQUENCE OF SCAFFOLD):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_convert_gff_to_genbank: cannot open $input_gff\n";
   print INPUT_GFF "HS04636\tsource\tgene\t966\t2995\t.\t+\t.\tmygene\n";
   print INPUT_GFF "HS04636\tsource\tmRNA\t966\t2995\t.\t+\t.\tID=mymrna;Parent=mygene\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t966\t1017\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t1818\t1934\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2055\t2198\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   print INPUT_GFF "HS04636\tsource\tCDS\t2852\t2995\t.\t+\t.\tID=cds1;Parent=mymrna\n";
   close(INPUT_GFF);
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_convert_gff_to_genbank: cannot open $assembly\n";
   print ASSEMBLY ">HS04636\n";
   print ASSEMBLY "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT\n";
   close(ASSEMBLY);
   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($input_gff); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:
   delete($SEQ->{"HS04636"});
   ($errorcode,$errormsg)  = &convert_gff_to_genbank($outputdir,$GENES,$EXONS,$LEN,$SEQ,1);
   if ($errorcode != 9) { print STDERR "ERROR: test_convert_gff_to_genbank: failed test4 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $input_gff";
   system "rm -f $output_genbank";
   system "rm -f $assembly";

 
}

#------------------------------------------------------------------#

# READ IN THE INPUT GFF FILE, AND MAKE THE OUTPUT GENBANK FILE:

sub convert_gff_to_genbank
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $GENES               = $_[1]; # HASH TABLE OF GENES ON SCAFFOLDS
   my $EXONS               = $_[2]; # HASH TABLE OF EXONS IN GENES
   my $LEN                 = $_[3]; # HASH TABLE WITH THE LENGTHS OF SCAFFOLDS
   my $SEQ                 = $_[4]; # HASH TABLE WITH SCAFFOLD SEQUENCES
   my $testing             = $_[5]; # SAYS WHETHER THIS IS CALLED FROM A TESTING SUBROUTINE 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $scaffold;                    # NAME OF A SCAFFOLD
   my $scaffold_len;                # LENGTH OF SCAFFOLD 
   my $output;                      # OUTPUT FILE FOR A SCAFFOLD
   my $genes;                       # GENES ON A SCAFFOLD
   my @genes;                       # GENES ON A SCAFFOLD
   my $gene;                        # A GENE
   my $i;                           # 
   my $exons;                       # EXONS IN A GENE
   my @exons;                       # EXONS IN A GENE
   my $exon;                        # AN EXON IN A GENE
   my $exon_strand;                 # STRAND OF EXON 
   my $j;                           # 
   my $exon_start;                  # START OF AN EXON
   my $exon_end;                    # END OF AN EXON
   my $new_exons;                   # EXONS IN A GENE, FORMATTED TO BE PRINTED OUT
   my @temp;                        # 
   my $seq;                         # SCAFFOLD SEQUENCE
   my $num_as;                      # NUMBER OF As IN SCAFFOLD
   my $num_cs;                      # NUMBER OF Cs IN SCAFFOLD
   my $num_gs;                      # NUMBER OF Gs IN SCAFFOLD
   my $num_ts;                      # NUMBER OF Ts IN SCAFFOLD
   my $num_other;                   # NUMBER OF NON-ACGT LETTERS IN SCAFFOLD
   my $a_line;                      # A LINE OF SEQUENCE TO PRINT OUT
   my $a_chunk;                     # A CHUNK OF SEQUENCE TO PRINT OUT
   my $offset;                      # OFFSET TO USE FOR PRINTING OUT SEQUENCE
   my $mrna_start;                  # START OF A mRNA
   my $mrna_end;                    # END OF A mRNA

   # LOOP THROUGH ALL THE SCAFFOLDS WITH GENES ON THEM: 
   foreach $scaffold (keys %{$GENES})
   {
      # FIND THE LENGTH OF THE SCAFFOLD:
      if (!($LEN->{$scaffold}))
      {
         $errormsg         = "ERROR: convert_gff_to_genbank: do not know length of scaffold $scaffold\n";
         $errorcode        = 7; # ERRORCODE=7 (TESTED FOR) 
         return($errorcode,$errormsg);
      }
      $scaffold_len        = $LEN->{$scaffold};

      # MAKE AN OUTPUT FILE FOR THIS SCAFFOLD:
      $output              = $scaffold.".gb";
      $output              = $outputdir."/".$output;
      system "rm -f $output"; # DELETE THE OUTPUT FILE IF IT EXISTS ALREADY

      # SEE FOR GENBANK FORMAT DEFINITION: http://www.bioperl.org/wiki/GenBank_sequence_format
      if ($testing == 0) { print STDERR "Opening output file $output...\n"; }
      open(OUTPUT,">$output") || die "ERROR: convert_gff_to_genbank: cannot open $output\n";
      print OUTPUT "LOCUS       $scaffold   $scaffold_len bp  DNA\n";
      print OUTPUT "FEATURES             Location/Qualifiers\n";
      print OUTPUT "     source          1..$scaffold_len\n";

       # FIND THE GENES IN THE SCAFFOLD:
       if (!($GENES->{$scaffold}))
       {
          $errormsg         = "ERROR: convert_gff_to_genbank: do not know genes for scaffold $scaffold\n";
          $errorcode        = 8; # ERRORCODE=8 (SHOULDN'T GET HERE, SO CAN'T TEST FOR THIS) 
          return($errorcode,$errormsg);
       }
       $genes               = $GENES->{$scaffold};
       @genes               = split(/\,/,$genes);
       for ($i = 0; $i <= $#genes; $i++)
       {
          $gene             = $genes[$i];
          # GET THE EXONS IN THIS GENE:
          if (!($EXONS->{$gene}))
          {
             $errormsg      = "ERROR: convert_gff_to_genbank: do not know exons in gene $gene\n";
             $errorcode     = 3; # ERRORCODE=3 (TESTED FOR) 
             return($errorcode,$errormsg);
          }
          $exons            = $EXONS->{$gene};
          @exons            = split(/\,/,$exons);
          $new_exons        = "";
          for ($j = 0; $j <= $#exons; $j++)
          {
             $exon          = $exons[$j];
             @temp          = split(/=/,$exon);
             $exon_start    = $temp[0];
             $exon_end      = $temp[1];
             $exon_strand   = $temp[2];
             $new_exons     = $new_exons.",".$exon_start."..".$exon_end;
          }
          $new_exons        = substr($new_exons,1,length($new_exons)-1);
          if ($#exons > 0) { $new_exons = "join(".$new_exons.")"; }  
          if ($exon_strand eq '-') { $new_exons = "complement(".$new_exons.")";}

          print OUTPUT "     CDS             $new_exons\n";
      }
      # PRINT OUT THE FASTA SEQUENCE:
      if (!($SEQ->{$scaffold}))
      {
         $errormsg         = "ERROR: convert_gff_to_genbank: do not know sequence for $scaffold\n";
         $errorcode        = 9; # ERRORCODE=9 (TESTED FOR) 
         return($errorcode,$errormsg);
      }
      $seq                 = $SEQ->{$scaffold};
      $seq                 =~ tr/[A-Z]/[a-z]/;
      ($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg) = &count_bases($seq);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
      print OUTPUT "BASE COUNT     $num_as a   $num_cs c  $num_gs g  $num_ts t\n";
      print OUTPUT "ORIGIN\n";
      # PRINT OUT THE SEQUENCE, eg.:
      #    1 gagctcacat taactattta cagggtaact gcttaggacc agtattatga ggagaattta
      #   61 cctttcccgc ctctctttcc aagaaacaag gagggggtga aggtacggag aacagtattt
      #  121 cttctgttga aagcaactta gctacaaaga taaattacag ctatgtacac tgaaggtagc
      # ...
      # 9421 aaaaaaaaaa aaaaatcgat gtcgactcga gtc	
      $offset              = 0;
      my $index            = 1;
      while ($offset < $scaffold_len)
      {
         printf OUTPUT "%9d", $index ; # PRINT $index RIGHT-JUSTIFIED, IN A REGION 9 SPACES WIDE
         $a_line           = substr($seq,$offset,60); 
         my $length_a_line = length($a_line); 
         my $offset2       = 0; 
         while ($offset2 < $length_a_line)
         {
            $a_chunk       = substr($a_line,$offset2,10);
            print OUTPUT " $a_chunk";
            $offset2       = $offset2 + 10; 
         }
         print OUTPUT "\n";
         $offset           = $offset + 60;
         $index            = $index + $length_a_line;
      }
      print OUTPUT "//\n"; 
      close(OUTPUT);
   }   
 
   return($errorcode,$errormsg);   
}

#------------------------------------------------------------------#

# TEST &read_gene_positions

sub test_read_gene_positions
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $GENES;                       # HASH TABLE TO STORE THE GENES ON SCAFFOLDS
   my $EXONS;                       # HASH TABLE TO STORE THE EXONS IN GENES 
   my $gff;                         # GFF FILE 
   my $sorted_gff;                  # SORTED GFF FILE 
 
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(GFF);
   # SORT THE GFF FILE SO THAT WE HAVE gene,mRNA,CDS FOR EACH GENE:
   ($sorted_gff,$errorcode,$errormsg) = &sort_input_gff($gff);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }    
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($sorted_gff);
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   if ($GENES->{"Ratt_Curated.Sratti_scf00001_Chr00"} ne "ratti_train3,ratti_train4") { print STDERR "ERROR: test_read_gene_positions: failed test1b\n"; exit;}
   if ($GENES->{"Ratt_Curated.Sratti_scf00001_Chr01"} ne "ratti_train414") { print STDERR "ERROR: test_read_gene_positions: failed test1c\n"; exit;}
   if ($EXONS->{"ratti_train3"} ne "676104=676526=+,676766=677515=+,677568=677702=+") { print STDERR "ERROR: test_read_gene_positions: failed test1d\n"; exit;}
   if ($EXONS->{"ratti_train4"} ne "699685=699951=-,700002=700943=-") { print STDERR "ERROR: test_read_gene_positions: failed test1e\n"; exit;}
   if ($EXONS->{"ratti_train414"} ne "1650617=1650755=+,1650801=1651867=+") { print STDERR "ERROR: test_read_gene_positions: failed test1f\n"; exit;}
   system "rm -f $gff";
   system "rm -f $sorted_gff";

   # TEST READING IN A GFF FILE LIKE THE TRAINING GENES:
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00008.1453643	source	CDS	115993	116307	.	+	.	T1.PTRK.scaffold.00008.1453643_1_0;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print GFF "PTRK.scaffold.00008.1453643	source	CDS	116360	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_1;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print GFF "PTRK.scaffold.00008.1453643	source	gene	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1\n";
   print GFF "PTRK.scaffold.00008.1453643	source	mRNA	115993	116413	.	+	.	T1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1\n";
   print GFF "PTRK.scaffold.00008.1453643	source	CDS	535391	535554	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_0;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print GFF "PTRK.scaffold.00008.1453643	source	CDS	535599	537591	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_1;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print GFF "PTRK.scaffold.00008.1453643	source	CDS	537642	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_2;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print GFF "PTRK.scaffold.00008.1453643	source	gene	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1\n";
   print GFF "PTRK.scaffold.00008.1453643	source	mRNA	535391	537971	.	+	.	T2.PTRK.scaffold.00008.1453643.embl_1_mRNA;Parent=T2.PTRK.scaffold.00008.1453643.embl_1\n";
   close(GFF);
   # SORT THE GFF FILE SO THAT WE HAVE gene,mRNA,CDS FOR EACH GENE:
   ($sorted_gff,$errorcode,$errormsg) = &sort_input_gff($gff);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }    
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($sorted_gff);
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test3 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   if ($GENES->{"PTRK.scaffold.00008.1453643"} ne "T1.PTRK.scaffold.00008.1453643_1,T2.PTRK.scaffold.00008.1453643.embl_1") { print STDERR "ERROR: test_read_gene_positions: failed test3b\n"; exit;}
   if ($EXONS->{"T1.PTRK.scaffold.00008.1453643_1"} ne "115993=116307=+,116360=116413=+") { print STDERR "ERROR: test_read_gene_positions: failed test3c\n"; exit;}
   if ($EXONS->{"T2.PTRK.scaffold.00008.1453643.embl_1"} ne "535391=535554=+,535599=537591=+,537642=537971=+") { print STDERR "test_read_gene_positions: failed test3d\n"; exit;}
   system "rm -f $gff";
   system "rm -f $sorted_gff";

}

#------------------------------------------------------------------#

# READ IN THE POSITIONS OF GENES ON SCAFFOLDS:

sub read_gene_positions
{
   my $gff                 = $_[0]; # GFF FILE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my @temp;                        # 
   my $scaffold;                    # SCAFFOLD NAME
   my $feature;                     # FEATURE TYPE
   my $start;                       # START 
   my $end;                         # END 
   my $strand;                      # STRAND
   my $name                = "none";# FEATURE NAME
   my $score;                       # FEATURE SCORE
   my %EXONS               = ();    # HASH TABLE TO STORE EXONS IN A GENE
   my %GENES               = ();    # HASH TABLE TO STORE GENES ON A SCAFFOLD 
   my $exon;                        # COORDINATES OF AN EXON 
   my %SEEN                = ();    # HASH TABLE TO RECORD WHETHER WE SAW A GENE BEFORE 

   # READ IN THE GFF FILE:
   # NOTE: THE GFF FILE HAS A 'gene' FEATURE FIRST FOR THE GENE, AND THIS
   #       IS FOLLOWED BY CDS & mRNA FEATURES FOR THE GENE (THE sort_gff FUNCTION ENSURES THIS).
   open(GFF,"$gff") || die "ERROR: read_gene_positions: cannot open $gff\n";
   while(<GFF>)
   {
      $line                = $_;   
      chomp $line;
      if (substr($line,0,1) ne '#') # IF IT'S NOT A COMMENT LINE
      {
         @temp             = split(/\t+/,$line);
         $scaffold         = $temp[0];
         $feature          = $temp[2];
         $start            = $temp[3];
         $end              = $temp[4];
         $score            = $temp[5];
         $strand           = $temp[6];
         if    ($feature eq 'gene')
         {
            # FIND THE GENE NAME:
            $name          = $temp[8];
            if ($name =~ /ID=/)
            { 
               @temp       = split(/ID=/, $name); # eg. ID=ratti_train414
               $name       = $temp[1];            # eg. ratti_train414 
            }
            if (!($GENES{$scaffold}))
            {
               $GENES{$scaffold} = $name;
            }
            else
            {
               $GENES{$scaffold} = $GENES{$scaffold}.",".$name;
            }
         }
         elsif ($feature eq 'CDS' || $feature eq 'cds' || $feature eq 'First' || $feature eq 'Internal' || $feature eq 'Terminal' || $feature eq 'Single')
         # exonerate USES 'cds', cegma 'First'/'Terminal'/'Internal'/'Single'
         {
            if ($name eq 'none')
            {
               $errormsg   = "ERROR: read_gene_positions: name $name line $line in gff $gff name $name\n";
               $errorcode  = 14; # ERRORCODE=14 (NOTE: THIS SHOULDN'T HAPPEN, SO CAN'T TEST FOR) 
               return(\%GENES,\%EXONS,$errorcode,$errormsg);
            }
            $exon       = $start."=".$end."=".$strand;
            if (!($EXONS{$name}))
            {
               $EXONS{$name} = $exon;
            }
            else
            {
               $EXONS{$name} = $EXONS{$name}.",".$exon;
            }
         }
      }
   }
   close(GFF);

   return(\%GENES,\%EXONS,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &read_assembly

sub test_read_assembly
{
   my $outputdir           = $_[0]; # DIRECTORY WHERE WE CAN WRITE OUTPUT FILES
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES
   my $assembly;                    # TEMPORARY ASSEMBLY FILE NAME 
   my $SEQ;                         # HASH TABLE WITH SEQUENCES OF SCAFFOLDS
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR 

   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_assembly: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY "\n";
   print ASSEMBLY ">seq3\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY "TTTTT\n";
   print ASSEMBLY ">seq4\n";
   print ASSEMBLY ">seq5\n";
   print ASSEMBLY " AAA AA \n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($SEQ->{'seq1'} ne 'AAAAA' || $SEQ->{'seq2'} ne 'AAAAATTTTT' || $SEQ->{'seq3'} ne 'AAAAATTTTT' || defined($SEQ->{'seq4'}) || 
       $SEQ->{'seq5'} ne 'AAAAA' || $errorcode != 0) 
   { 
      print STDERR "ERROR: test_read_assembly: failed test1\n"; 
      exit;
   }
   system "rm -f $assembly";
 
   # TEST ERRORCODE=4:
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 4) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";

   # TEST ERRORCODE=5:
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 5) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";

}

#------------------------------------------------------------------#

# READ IN THE ASSEMBLY FILE:
# SUBROUTINE SYNOPSIS: read_assembly(): read sequences in a fasta file into a hash

sub read_assembly       
{
   my $input_assembly      = $_[0]; # THE INPUT ASSEMBLY FILE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my $scaffold            = "none";# NAME OF SCAFFOLD
   my @temp;                        # 
   my $seq;                         # SEQUENCE OF SCAFFOLD 
   my %SEQ                 = ();    # HASH TABLE FOR STORING SCAFFOLD SEQUENCE

   $seq                    = "";
   open(ASSEMBLY,"$input_assembly") || die "ERROR: read_assembly: cannot open $input_assembly\n";
   while(<ASSEMBLY>)
   {
      $line                = $_;
      chomp $line;
      if (substr($line,0,1) eq '>')
      {
         if ($seq ne "")
         {
            if ($scaffold eq 'none')
            {
               $errormsg   = "ERROR: read_assembly: do not know name of scaffold\n";
               $errorcode  = 5; # ERRORCODE=5 (TESTED FOR)
               return(\%SEQ,$errorcode,$errormsg);
            }
            if ($SEQ{$scaffold}) 
            {
               $errormsg   = "ERROR: read_assembly: already know sequence for scaffold $scaffold\n";
               $errorcode  = 4; # ERRORCODE=4 (TESTED FOR)
               return(\%SEQ,$errorcode,$errormsg);
            }
            $SEQ{$scaffold}= $seq;
            $seq           = "";
         }
         @temp             = split(/\s+/,$line);
         $scaffold         = $temp[0];
         $scaffold         = substr($scaffold,1,length($scaffold)-1);
      }
      else
      {
         $line             =~ s/\s+//g; # REMOVE SPACES
         $seq              = $seq.$line; 
      }
   }
   close(ASSEMBLY); 
   if ($seq ne "")
   {
      if ($scaffold eq 'none')
      {
         $errormsg         = "ERROR: read_assembly: do not know name of scaffold\n";
         $errorcode        = 5; # ERRORCODE=5 (TESTED FOR)
         return(\%SEQ,$errorcode,$errormsg);
      }

      if ($SEQ{$scaffold}) 
      {
         $errormsg         = "ERROR: read_assembly: already know sequence for scaffold $scaffold\n";
         $errorcode        = 4; # ERRORCODE=4 (TESTED FOR)
         return(\%SEQ,$errorcode,$errormsg);
      }
      $SEQ{$scaffold}      = $seq;
      $seq                 = "";
   }

   return(\%SEQ,$errorcode,$errormsg);

}

#------------------------------------------------------------------#

# READ IN THE LENGTHS OF CONTIGS IN THE ASSEMBLY:
# SUBROUTINE SYNOPSIS: read_scaffold_length(): store lengths of sequences from a fasta file in a hash

sub read_scaffold_lengths
{
   my $assembly            = $_[0]; # FILE CONTAINING THE ASSEMBLY
   my %SCAFFOLDLEN         = ();    # HASH TABLE FOR STORING LENGTHS OF SCAFFOLDS
   my $line;                        #  
   my @temp;                        # 
   my $length;                      # LENGTH OF A SCAFFOLD
   my $name                = "none";# NAME OF A SCAFFOLD
   my $seq;                         # SEQUENCE OF A SCAFFOLD
   my $errorcode           = 0;     # THIS IS RETURNED AS 0 IF NO ERROR OCCURRED. 
   my $errormsg            = "none";# THIS IS RETURNED AS 'none' IF NO ERROR OCCURRED.
   my %SEEN                = ();    # HASH TABLE TO RECORD WHICH SCAFFOLD NAMES WE HAVE SEEN.
 
   open(ASSEMBLY,"$assembly") || die "ERROR: read_scaffold_lengths: cannot open $assembly\n"; 
   while(<ASSEMBLY>)
   {
      $line                = $_;
      chomp $line;
      if (substr($line,0,1) eq '>')
      {
         @temp             = split(/\s+/,$line);
         $name             = $temp[0];
         $name             = substr($name,1,length($name)-1);
         if ($SEEN{$name})
         {
            $errormsg      = "ERROR: read_scaffold_lengths: seen scaffold $name already\n"; 
            $errorcode     = 1; # ERRORCODE=1 (TESTED FOR)
            return(\%SCAFFOLDLEN, $errorcode, $errormsg);
         }
         $SEEN{$name}      = 1;
      } 
      else
      {
         $seq              = $line;
         # REMOVAL SPACES FROM THE LINE, EITHER INTERNAL SPACES OR SPACES AT EITHER END:
         $seq              = $line;
         $seq              =~ s/\s+//g;
         # STORE THE LENGTH OF THE SEQUENCE, OR UPDATE THE STORED LENGTH:
         if ($seq eq '') { $length = 0;           }
         else            { $length = length($seq);}
         if ($name eq 'none') 
         { 
            $errormsg      = "ERROR: read_scaffold_lengths: name is $name\n"; 
            $errorcode     = 2; # ERRORCODE=2 (TESTED FOR)
            return(\%SCAFFOLDLEN, $errorcode, $errormsg);
         }
         if (!($SCAFFOLDLEN{$name})){ $SCAFFOLDLEN{$name} = $length;}
         else {$SCAFFOLDLEN{$name} =  $SCAFFOLDLEN{$name} + $length;}   
      }
   }
   close(ASSEMBLY);

   return(\%SCAFFOLDLEN,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &read_scaffold_lengths

sub test_read_scaffold_lengths
{
   my $outputdir           = $_[0];  # DIRECTORY TO WRITE OUTPUT IN.
   my $assembly;                     # FILE CONTAINING ASSEMBLY.      
   my $LEN;                          # HASH TABLE CONTAINING LENGTHS OF SEQUENCES.
   my $errorcode;                    # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR. 
   my $errormsg;                     # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR.

   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY "\n";
   print ASSEMBLY ">seq3\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY "TTTTT\n";
   print ASSEMBLY ">seq4\n";
   print ASSEMBLY ">seq5\n";
   print ASSEMBLY " AAA AA \n";
   close(ASSEMBLY);
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($LEN->{'seq1'} != 5 || $LEN->{'seq2'} != 10 || $LEN->{'seq3'} != 10 || defined($LEN->{'seq4'}) || $LEN->{'seq5'} != 5 || $errorcode != 0) { print STDERR "ERROR: test_read_scaffold_lengths: failed test1\n"; exit;}
   system "rm -f $assembly";

   # TEST ERRORCODE=1:
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 1) { print STDERR "ERROR: test_read_scaffold_lengths: failed test2\n"; exit;}
   system "rm -f $assembly";

   # TEST ERRORCODE=2:
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 2) { print STDERR "ERROR: test_read_scaffold_lengths: failed test3\n"; exit;}
   system "rm -f $assembly";
   
}

#------------------------------------------------------------------#

# SUBROUTINE TO MAKE A FILE NAME FOR A TEMPORARY FILE:

sub make_filename
{
   my $outputdir             = $_[0]; # OUTPUT DIRECTORY TO WRITE TEMPORARY FILE NAME TO
   my $found_name            = 0;     # SAYS WHETHER WE HAVE FOUND A FILE NAME YET
   my $filename              = "none";# NEW FILE NAME TO USE 
   my $errorcode             = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg              = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $poss_filename;                 # POSSIBLE FILE NAME TO USE
   my $random_number;                 # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
 
   while($found_name == 0)
   {
      $random_number         = rand();
      $poss_filename         = $outputdir."/tmp".$random_number;
      if (!(-e $poss_filename))
      {
         $filename           = $poss_filename;
         $found_name         = 1;
      } 
   } 
   if ($found_name == 0 || $filename eq 'none')
   {
      $errormsg              = "ERROR: make_filename: found_name $found_name filename $filename\n";
      $errorcode             = 6; # ERRORCODE=6 
      return($filename,$errorcode,$errormsg);
   }

   return($filename,$errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &print_error

sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR

   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 

   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}

#------------------------------------------------------------------#

# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT

   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }

   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }

   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 

   return($errormsg,$errorcode);
}

#------------------------------------------------------------------#



