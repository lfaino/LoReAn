package HelminthGenomeAnalysis::AvrilFileUtils;

use strict;
use warnings;
use Math::Round; # HAS THE nearest() FUNCTION
use Carp::Assert; # HAS THE assert() FUNCTION 
use Scalar::Util qw(looks_like_number);

use base 'Exporter';
our @EXPORT_OK = qw( make_filename check_if_files_are_identical write_array_to_file );

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: make_filename(): SUBROUTINE TO MAKE A FILE NAME FOR A TEMPORARY FILE:

sub make_filename
{
   my $outputdir             = $_[0]; # OUTPUT DIRECTORY TO WRITE TEMPORARY FILE NAME TO
   my $found_name            = 0;     # SAYS WHETHER WE HAVE FOUND A FILE NAME YET
   my $filename              = "none";# NEW FILE NAME TO USE 
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

   # DIE IF WE HAVE NOT FOUND A FILE NAME: 
   assert($found_name != 0 && $filename ne 'none'); # THIS SHOULD NEVER HAPPEN, PROGRAM WILL DIE IF IT DOES.

   return($filename); 
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: check_if_files_are_identical(): SUBROUTINE TO CHECK WHETHER TWO FILES ARE IDENTICAL:

sub check_if_files_are_identical
{
   my $file                  = $_[0]; # FILE WE WANT TO TEST
   my $expected_file         = $_[1]; # FILE THAT $file SHOULD BE IDENTICAL TO
   my $differences;                   # DIFFERENCES BETWEEN $file AND $expected_file 
   my $line;                          # 
   my $length_differences;            # LENGTH OF $differences
   my $files_are_identical   = 0;     # SAYS WHETHER THE TWO FILES ARE IDENTICAL OR NOT

   # THROW AN EXCEPTION IF $file DOES NOT EXIST:
   throw Error::Simple("ERRORCODE=1: check_if_files_are_identical: file $file does not exist") if !(-e "$file"); 
   # THROW AN EXCEPTION IF $expected_file DOES NOT EXIST:
   throw Error::Simple("ERRORCODE=2: check_if_files_are_identical: file $expected_file does not exist") if !(-e "$file"); 
 
   $differences              = "";
   open(TEMP,"diff $file $expected_file |");
   while(<TEMP>)
   {
      $line                  = $_;
      $differences           = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);

   if ($length_differences == 0) { $files_are_identical = 1;} # FILES ARE IDENTICAL
   else                          { $files_are_identical = 0;} # FILES ARE NOT IDENTICAL

   return($files_are_identical);

}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: write_array_to_file(): SUBROUTINE TO WRITE AN ARRAY TO A FILE, AS A SINGLE COLUMN.

sub write_array_to_file
{
   my $file                  = $_[0]; # NAME OF FILE TO WRITE TO
   my $array                 = $_[1]; # POINTER TO ARRAY
   my $i;                             # 

   open(FILE,">$file") || die "ERROR: write_array_to_file: cannot open $file\n";
   for ($i = 0; $i < @$array; $i++)
   {
      print FILE "$array->[$i]\n";
   }
   close(FILE);
   
   return(1);
}

#------------------------------------------------------------------#

1;
