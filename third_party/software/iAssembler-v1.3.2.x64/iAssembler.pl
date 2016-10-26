#!/usr/bin/perl -w
=head1 NAME

 iAssembler - 

=head1 DESCRIPTION
 
 Author : Yi zheng
 E-mail : bioinfo@cornell.edu
 Update : 07/20/2012

=head1 UPDATE INFO

* iAssembler v1.32 - 07/20/12. Changes from previous version:
	1. Fixed a bug in parsing megablast result (query length info)

* iAssembler v1.31 - 03/28/12. Changes from previous version:
        1. support MIRA V3.4.0. iAssembler will not support the old version of MIRA from this update.

* iAssembler v1.30 - 05/04/11. Changes from previous version:
	1. Add a function to correct unigene base errors
	2. Add headers to the output SAM file:

* iAssembler v1.22 - 03/28/10. Changes from previous version:
        1. Fixed -h parameters that previously was not working on MIRA and cap3.
	2. Remove 10bp buffer settings for -e parameter
	3. Fixed failed-run problem for running the same project at the second time 
	4. Fixed failed-run problem of mira caused by high megahub rate, set mnr=yes.

* iAssembler v1.21 - 12/02/10. Changes from previous version:
	1. Changed minimum number of -e parameters to 6, to be compatible with cap3.

* iAssembler v1.2 - 08/02/10. Changes from previous version:
	1. Compatible with MIRA version 3.x

* iAssembler v1.1 - 06/23/10. Changes from previous version:
	1. Fixed the error that caused EST clustering to fail for datasets containing highly redundant sequences
	2. Fixed several other small bugs

* iAssembler v1.0 - 05/21/10. Changes from previous version:
	1. Added an output file in SAM format. The file contains the alignment information of each sequence read to its corresponding unigene and can be views by several visualization programs such as Tablet and IGV.
	2. Combined percent identity cutoff for clustering (-x) and assembly (-p) into a single parameter (-p). Parameter -x is disabled
	3. Disabled clustering using blastn. Currently only megablast is used for clustering. Parameter -b now has different meaning (see below)
	4. Added -b parameter which specifies the number of threads used for MIRA assembly program
	5. Added -d parameter to control whether to generate program log files

* iAssembler v1.0 (beta) - 04/13/10


#######################################################################
 May 06 2010: 3 bugs fixed
 1. endless loop at error correct step.
 2. using same percent identity parameter for clustering and assembling.
 3. pass long sequence to next step.
 4. using megablast in perl assembler
 5. using default parameter run megablast
 6. delete temp file when assemble using mira cap3 and perl
#######################################################################

=cut
use strict;
use Cwd;
use File::Basename;
use IO::File;
use Getopt::Std;

my $usage = q'
VERSION: v1.32
USAGE: 
	Perl  iAssembler.pl  [parameters]

	Section 1: Input parameters
	-i  [String]   Name of the input sequence file in FASTA format (required) 
	-q  [String]   Name of the quality file in FASTA format (default: none)

	Section 2: Assembly parameters
	-a  [Integer]  number of CPUs used for megablast clustering (default = 1)
	-b  [Integer]  number of CPUs used for MIRA assembly program (default = 1)
	-e  [Integer]  maximum length of end clips (6~100; default = 30)
	-h  [Integer]  minimum overlap length (>=30; default = 40)
	-p  [Integer]  minimum percent identify for sequence clustering and assembly (95~100; default = 97) 

	Section 3: Output parameters
	-u  [String]   prefix used for IDs of the assembled unigenes (default = UN)
		       iAssembler names the resulted unigenes with a prefix and trailing numbers, e.g., UN00001
	-l  [Integer]  length of the trailing numbers in unigene IDs (>= default; defalut = number characters of the maximum number assigned to unigenes)
		       For example, if the maximum trailing number assigned to the resulted unigenes is 5000, then the default of -l is 4 (\'5000\' has 5 characters). In this case users can set a number greater than or equal to 4.
	-s  [Integer]  start number of unigene ID trailing number (>= 1; default = 1)
	-o  [String]   Name of the output directory (default = "input file name" + "_output") 
	-d  	       Produce log files. With this parameter will produce log files in the output folder. 
';

#################################################################
# Input Parameters						#
#################################################################
my %options;

getopts('i:q:a:b:e:h:p:u:l:s:o:d', \%options) || die "$usage\nError getting options!";

# Section 1: Input parameters
my $fasta_file  = $options{i} || die $usage."\nYou must provide input file (-i)!\n\n";
my $quality_file= $options{q} || "";

# Section 2: Assembly parameters
my $cpus        = $options{a} || 1;
my $cpus_mira   = $options{b} || 1;
my $max_end_clip= $options{e} || 30;
my $min_overlap = $options{h} || 40;
my $min_identity= $options{p} || 97;

# Section 3: Output parameters
my $id_prefix   = $options{u} || "UN";
my $id_length   = $options{l} || 0;
my $id_start    = $options{s} || 1;
my $output_dir  = $options{o} || $fasta_file."_output";
my $debug_mode  = $options{d} || 0;

#################################################################
# start report                                                  #
#################################################################
my $start_time = localtime();
print qq'
iAssembler starts at: $start_time
';

#################################################################
# replace the parameters using manual config file               #
#################################################################
my ($blast_program, $blast_param);
$blast_program = "megablast";
$blast_param = "-F F -a $cpus -e 1e-5 -W 20";

#################################################################
# Check input parameters & Store them to %ENV                   #
#################################################################
if ($cpus <= 0 || $cpus > 99 )
{ die $usage."\n Error at number of CPUs used for blast program (default = 1, From 1 to 99)"; }
if ($cpus_mira <= 0 || $cpus_mira > 99 )
{ die $usage."\n Error at number of CPUs used for MIRA program (default = 1, From 1 to 99)"; }

if (    $max_end_clip <= 100 &&
        $max_end_clip >= 6 &&
        $min_overlap  >= 30 &&
        $min_identity >= 95 &&
        $min_identity <= 100  )
{
        $ENV{'max_end_clip'} = $max_end_clip;
        $ENV{'min_overlap' } = $min_overlap;
        $ENV{'min_identity'} = $min_identity;
}
else { die "Error! Input parameters have some errors!\n"; }

if ($id_length > 15) { die "The length of ID is too long.\n"; }
if ($id_start  < 1){ die "The start of ID must be positive int.\n"; }
if ($debug_mode != 0 && $debug_mode != 1) { die "Error at debug_mode.\n"; }
if ($output_dir eq "" && !$output_dir) { die "Error at output folder.\n"; }
if ($blast_program ne "megablast" ) {die "Error at megablast.\n"; }

$ENV{'debug_mode'} = $debug_mode;

#################################################################
# store megablast parameters to $ENV                            #
#################################################################
$ENV{'blast_program'} = $blast_program;  #blastn or megablast;
$ENV{'blast_param'}   = $blast_param;

#################################################################
# Set directory for iAssembler                                  #
# 1. current_dir: the path to executive command                 #
# 2. working_dir: the path for temp files                       #
# 3. output_dir : the path for output result                    #
# 4. program_dir: the path of iAssembler                        #
# 5. program_bin_dir: the path of other assemblers and script   #
# 6. log_dir: the path of log files                             #
#################################################################
my $current_dir = getcwd;
my $working_dir = $current_dir."/".$fasta_file."_Assembly";
$output_dir  = $current_dir."/".$output_dir;
my $log_dir  = $output_dir."/log";
my $log_file = $log_dir."/main.log";

my ($program_dir, $program_bin_dir);
if ($0 =~ m{^/}) { $program_dir = dirname($0);}
else { $program_dir = dirname("$current_dir/$0"); }
$program_bin_dir = $program_dir."/bin";
$ENV{'PATH'} = $program_dir."/bin".":".$ENV{'PATH'};

if ($working_dir =~ m/\/\// || $output_dir =~ m/\/\// || $log_dir =~ m/\/\//) { die "Error at some folders\n"; }

if (-e $working_dir) { system ("rm -rf $working_dir") && die "Can not delete exist output directroy: $working_dir . $!\n"; }
mkdir($working_dir) || die "Can not creat output directory! $!\n";

if (-e $output_dir)  { system ("rm -rf $output_dir") && die "Can not delete exist output directroy: $output_dir . $!\n"; }
mkdir($output_dir) || die "Can not creat output directory! $!\n";

if ($debug_mode == 1) { mkdir($log_dir) || die "Can not creat log directory! $!\n"; }

$ENV{"current_dir"} = $current_dir;
$ENV{"working_dir"} = $working_dir;
#$ENV{"program_dir"} = $program_dir;
$ENV{"program_bin_dir"} = $program_bin_dir;
$ENV{"log_dir"} = $log_dir;
$ENV{"log_file"} = $log_file;

#################################################################
# the initial iAssembler verion                                 #
#################################################################
$ENV{"mira_version"} = "2.9.43";

#foreach my $a (sort keys %ENV){ print "$a \t $ENV{$a} \n"; } die;
#########################################################################
# Produce Auto Pipeline File						#
#########################################################################
my ($pipeline_file, $pipeline_manual);

if ($pipeline_manual)
{
	$pipeline_file = $pipeline_manual;
}
else
{
	$pipeline_file = $working_dir."/pipeline.auto";
	if (-s $pipeline_file) { unlink($pipeline_file);}
	my $paf = IO::File->new(">".$pipeline_file) || die "Can no open auto pipeline file $!\n";

	#########################################################
	# decide which version of mira                          #
	#########################################################
	my $version_file = $working_dir."/mira_version";
	system("mira | head > $version_file");

	my $version_number;

	my $mfh = IO::File->new($version_file) || die "Can not open mira version file\n";
	
	while(<$mfh>)
	{
		if (/This is MIRA/)
		{
			my @vm = split(/\s+/, $_);
			$version_number = $vm[3];
		}
	}
	$mfh->close;

	write_log("\nUsing MIRA $version_number\n\n") if $debug_mode == 1;

	my $mira_size = -s "mira";
	my $mira_min_overlap = $min_overlap-1;
	# mira 2.9.43 x32 and x64
	if ( $version_number =~ /V2\.9/ )
	{
		$ENV{"mira_version"} = $version_number;
		print $paf qq'
#################################################################
# iAssembly Pipeline Config File V0.1				#
# update: 08-26-2009						#
# Author: Yi Zheng  yz357 at cornell.edu			#
#################################################################

#################################################################
# Format:							#
#								#
# Part    1 --- the order of subpipeline run			#
# Name    mira --- subpipeline name				#
# Program mira_pipeline -subpipeline program name		#
# Circle  0 --- 						#
# 	0 means run did not run this subpipeline		#
# 	1 means run once					#
# Parameter1	xxx						#
# Parameter2	xxx						#
# 	if have two parameters, it run 2 cycles, 		#
@	circle is disabled;					#
#################################################################


Part	1
Name	mira
Program	mira_pipeline
Circle	2
Parameter1	-job=denovo,est,normal,sanger -notraceinfo -GENERAL:kcim=yes,not=$cpus_mira -CO:fnicpst=yes -CL:cpat=no,pec=no,pvlc=no,qc=no,bsqc=no,mbc=no,emlc=yes,mlcr=0,smlc=0,emrc=yes,mrcr=0,smrc=0 -AL:mrs=$min_identity,mo=$mira_min_overlap -AS:mrl=30,bdq=30 -SK:mnr=yes,pr=$min_identity
Parameter2	-job=denovo,est,normal,sanger -notraceinfo -GENERAL:kcim=yes,not=$cpus_mira -CO:fnicpst=yes -CL:cpat=no,pec=no,pvlc=no,qc=no,bsqc=no,mbc=no,emlc=yes,mlcr=0,smlc=0,emrc=yes,mrcr=0,smrc=0 -AL:mrs=$min_identity,mo=$mira_min_overlap -AS:mrl=30,bdq=30 -SK:mnr=yes,pr=$min_identity
Parameter3	-job=denovo,est,normal,sanger -notraceinfo -GENERAL:kcim=yes,not=$cpus_mira -CO:fnicpst=yes -CL:cpat=no,pec=no,pvlc=no,qc=no,bsqc=no,mbc=no,emlc=yes,mlcr=0,smlc=0,emrc=yes,mrcr=0,smrc=0 -AL:mrs=$min_identity,mo=$mira_min_overlap -AS:mrl=30,bdq=30 -SK:mnr=yes,pr=$min_identity
Parameter4	-job=denovo,est,normal,sanger -notraceinfo -GENERAL:kcim=yes,not=$cpus_mira -CO:fnicpst=yes -CL:cpat=no,pec=no,pvlc=no,qc=no,bsqc=no,mbc=no,emlc=yes,mlcr=0,smlc=0,emrc=yes,mrcr=0,smrc=0 -AL:mrs=$min_identity,mo=$mira_min_overlap -AS:mrl=30,bdq=30 -SK:mnr=yes,pr=$min_identity
Input1	est1.fa
Output2	mira.cmf
Output1	mira.uniseq

Part	2
Name	cap3
Program	cap3_pipeline
Circle	1
Parameter1	-y $max_end_clip -p $min_identity -o $min_overlap -s 251 -f 6
Input1	mira.uniseq
Input2	mira.cmf
Output1	cap3.uniseq
Output2	cap3.cmf

';
	}
	# mira 3.0.5 x32 and x64
	elsif( $version_number =~ m/V3/ )
	{
		$ENV{"mira_version"} = $version_number;
		print $paf qq'
#################################################################
# iAssembler Pipeline Config File V0.1              		#	
# update: 08-26-2009                                        	#
# Author: Yi Zheng  yz357 at cornell.edu                    	#
#################################################################


#################################################################
# Format:                                                    	#
#                                                             	#
# Part    1 --- the order of subpipeline run              	#
# Name    mira --- subpipeline name                         	#
# Program mira_pipeline -subpipeline program name             	#
# Circle  0 ---                                             	#
#       0 means run did not run this subpipeline           	#
#       1 means run once                                   	#
# Parameter1    xxx                                  		#
# Parameter2    xxx                                    		#
#       if have two parameter, it run 2 cycles, 		#
#	circle is disabled;              			#
#################################################################

Part	1
Name	mira
Program	mira_pipeline
Circle	2
Parameter1	-job=denovo,est,accurate,454 -notraceinfo -noclipping -GE:not=$cpus_mira 454_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no:mrl=30 COMMON_SETTINGS -AS:nop=4 -SK:not=$cpus_mira:pr=$min_identity -CL:pec=no 454_SETTINGS -CL:cpat=no:pvlc=no:qc=no:bsqc=no:mbc=no:emlc=no:emrc=no:msvs=no -AL:mo=$mira_min_overlap:mrs=$min_identity
Parameter2	-job=denovo,est,accurate,454 -notraceinfo -noclipping -GE:not=$cpus_mira 454_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no:mrl=30 COMMON_SETTINGS -AS:nop=4 -SK:not=$cpus_mira:pr=$min_identity -CL:pec=no 454_SETTINGS -CL:cpat=no:pvlc=no:qc=no:bsqc=no:mbc=no:emlc=no:emrc=no:msvs=no -AL:mo=$mira_min_overlap:mrs=$min_identity
Parameter3	-job=denovo,est,accurate,454 -notraceinfo -noclipping -GE:not=$cpus_mira 454_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no:mrl=30 COMMON_SETTINGS -AS:nop=4 -SK:not=$cpus_mira:pr=$min_identity -CL:pec=no 454_SETTINGS -CL:cpat=no:pvlc=no:qc=no:bsqc=no:mbc=no:emlc=no:emrc=no:msvs=no -AL:mo=$mira_min_overlap:mrs=$min_identity
Parameter4	-job=denovo,est,accurate,454 -notraceinfo -noclipping -GE:not=$cpus_mira 454_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no:mrl=30 COMMON_SETTINGS -AS:nop=4 -SK:not=$cpus_mira:pr=$min_identity -CL:pec=no 454_SETTINGS -CL:cpat=no:pvlc=no:qc=no:bsqc=no:mbc=no:emlc=no:emrc=no:msvs=no -AL:mo=$mira_min_overlap:mrs=$min_identity
Input1	est1.fa
Output2	mira.cmf
Output1	mira.uniseq

Part	2
Name	cap3
Program	cap3_pipeline
Circle	1
Parameter1	-y $max_end_clip -p $min_identity -o $min_overlap -s 251 -f 6
Input1	mira.uniseq
Input2	mira.cmf
Output1	cap3.uniseq
Output2	cap3.cmf
';
	}
        else
        {
                die "iAssembler can not identify the version of MIRA, please check it.\n";
        }

	$paf->close;
}

#die;
#########################################################################
# Debug mode 								#
# if you want debug, set $debug_mode = 1				#
# It will not run pipline, just print the commond			#
#########################################################################

my $debug_line = "======================================================";

#convert input config file to a pipeline text 
#and then run the pipeline text.
my $pipeline_text = parse_pipeline_cfg($pipeline_file);

my $last_out = parse_pipeline_text($pipeline_text);

my @last_out = split(/#/, $last_out);
my $last_out_uniseq; my $last_out_cmf;
if ($last_out[1] && $last_out[2]) { $last_out_uniseq = $last_out[1]; $last_out_cmf = $last_out[2];}else{ die "Error! No Last Out!\n";}

#########################################################################
# process the results of perl assembly then correct errors		#
#########################################################################
my $b_unigene   = $working_dir."/b_unigene_seq.fasta";
my $b_cmf       = $working_dir."/b_contig_member";
my $b_mp        = $working_dir."/b_unigene_mp";
my $b_id_map    = $working_dir."/b_member_position_stat";
my $b_sam       = $working_dir."/b_unigene.sam";

my $last_command = "last_pipeline -e $fasta_file -u $last_out_uniseq -c $last_out_cmf -p $id_prefix -s $id_start -l $id_length -o $b_unigene -t $b_cmf -b $b_mp -d $b_id_map -a $b_sam";

system($last_command);

my $o_unigene	= $output_dir."/unigene_seq.fasta";
my $o_cmf	= $output_dir."/contig_member";
my $o_mp	= $output_dir."/unigene_mp";
my $id_map	= $output_dir."/member_position_stat";
my $o_sam	= $output_dir."/unigene.sam";

my $base_correct_command = "base_correct -u $b_unigene -i $b_sam -c $b_cmf -p $b_mp -o $o_unigene -s $o_sam -t $o_cmf -b $o_mp -d $id_map";
#print $last_command."\n";
system($base_correct_command);

my $end_time = localtime();
print qq'
iAssembler finishes at: $end_time
';

#delete Assembly Folder
system("rm -rf $working_dir");
chdir($current_dir);
unlink('formatdb.log');
###############################################################################
# kentnf: parse_pipeline_cfg, all line to $all_pipeline			      #
###############################################################################

=head1 SUBROUTINE PART

=head2 parse_pipeline_ctg

 Function: parse pipeline contig file, convert pipeline contig file to iAssembler recongized command.

 Input: pipeline config file (file name)

 Return:ouput is command that iAssembler can perform.

=cut
sub parse_pipeline_cfg
{
	my $usage = "USAGE: &get_pipeline(pipeline_config_file)\n\n";

	my $pipeline_file = shift || die $usage;
	# Step1 check the file exist;
	unless (-s $pipeline_file) { die "\nError! $pipeline_file does not exists! $usage";}

	# Step2 convert the contig file to pipeline text; 

	my $pf_fh = IO::File->new($pipeline_file) || die "Can not open $pipeline_file $!\n";

	my $all_pipeline;

	while(my $pf_line = <$pf_fh>)
	{
		if($pf_line !~ m/^#/ && $pf_line =~ /\w+/)
		{
			$all_pipeline.=$pf_line;
		}
	}
	$pf_fh->close;

	write_log("\nPipeline Content Begin\n$debug_line\n$all_pipeline\n$debug_line\nPipeline Content End\n\n") if $debug_mode == 1; 
	return $all_pipeline;
}

###############################################################################
# sub: parse_pipeline_cfg, all line to $all_pipeline                          #
###############################################################################
=head2 parse_pipeline_text


=cut
sub parse_pipeline_text
{
	my $usage = "USAGE: &parse_pipeline_text(pipeline_text)\n\n";

	my $all_pipeline = shift || die $usage;

	my @part = split(/Part/, $all_pipeline);

	my %pipeline = ();

	my $ipart;

	my $last_output;

	# pipeline text to hash, one part one member, order by part number
	for($ipart=1; $ipart<@part; $ipart++)
	{
		my @line = split(/\n/, $part[$ipart]);

		$line[0] =~ s/\t//;

		if (defined $pipeline{$line[0]})
		{
			die "Error! Repeat order at $line[0] $line[1]\n";
		}
		else
		{
			$pipeline{$line[0]} = $part[$ipart];
		}
	}

	my %all_inout = ();

	# every part run pipeline;
	foreach my $order (sort keys %pipeline)
	{
		# prepare
		my @pline = split(/\n/, $pipeline{$order});

		my $name; my $program; my $circle; my %para_hash = (); my $para_num = 0; 

		my %input = (); my %output = ();
	
		my $ip;

		for($ip=1; $ip<@pline; $ip++)
		{
			my @mm = split(/\t/, $pline[$ip]);
			
			#process part name
			if ($mm[0] eq "Name")
			{
				$name = $mm[1];
			}

			#using which program(pipeline) to process data
			if($mm[0] eq "Program")
			{
				$program = $mm[1];
			}

			#how many circles
			if($mm[0] eq "Circle")
                        {
                                $circle = $mm[1];
                        }

			#process the Parameters
			if ($mm[0] =~ m/Parameter/)
			{
				$para_num++;

				if (defined $para_hash{$mm[0]})
				{
					die "Error! Repeat parameter at para $name --> $mm[0] --> $para_hash{$mm[0]}\n";
				}
				else
				{
					$para_hash{$mm[0]} = $mm[1];
				}
			}

			#the input the input can repeat in the pipeline
			if ($mm[0] =~ m/Input/)
			{
				# input to hash
				if (defined $input{$mm[0]})
				{
					die "Error! Repeat parameter at input $name --> $mm[0] --> $input{$mm[0]}\n";
				}
				else
				{
					$input{$mm[0]} = $working_dir."/".$mm[1];
				}

			}

			#out put
			if ($mm[0] =~ m/Output/)
			{
				#check unique
				if (defined $all_inout{$mm[1]})
				{
					die "Error! Repeat all output $name --> $mm[0] --> $mm[1]\n"
				}
				else
				{
					$all_inout{$mm[1]} = 1
				}

				# output to hash
				if (defined $output{$mm[0]})
				{
					die "Error! Repeat parameter at output $name --> $mm[0] --> $mm[1]\n";
				}
				else
				{
					$output{$mm[0]} = $working_dir."/".$mm[1];
				}
			}

		}
		
		# == change input ==
		my $input_pk = &pack_hash(\%input);
		
		if ($order == 1)
		{
			if ($quality_file =~ m/\w+/)
			{
				$input_pk = "#$fasta_file#$quality_file";
			}
			else
			{
				$input_pk = "#$fasta_file";
			}
		}
		if ($program eq "member_position_pipeline")
		{
			$input_pk = $input_pk."#".$fasta_file;
		}

		# == change circle ==
		if ($para_num >=2)
		{
			$circle = 1;
		}

		# == change output == 
		my $output_pk = &pack_hash(\%output);

		# == change parameter ==
		my $para_pk = &pack_hash(\%para_hash);
		$para_pk =~ s/_/ZZZ/ig;
		$para_pk =~ s/\s+/_/ig;

		my $command = "$program 1$input_pk 2$output_pk $circle 4$para_pk";
		#print $command."\n";

		$last_output = $output_pk;

		write_log("Command: $command\n") if $debug_mode == 1;
		system($command) && die "Error at command: $command\n"; 
	}
	#print $last_output."\n";
	return $last_output;
}

#################################################################
# pack the hash to command line, for pass command to script	#
#################################################################
sub pack_hash
{
	my $usage = "USAGE: package_value = &pack_hash{ref_hash}\n";
	my $hash = shift || die $usage;	

	my %hash = %$hash;

	my $package = "";

	foreach my $key (sort keys %hash)
	{
		if ($hash{$key})
		{
			$package .= "#".$hash{$key};
		}
	}

	return $package;
}

#################################################################
# sub: replace underline with space for blast parameters	#
#################################################################
sub parse_underline
{
	my $char = shift;
	$char =~ s/_/ /ig;
	return $char;
}

#################################################################
# sub: delete Folder 						#
#################################################################
sub deleteFolder
{
	my $path = shift;
	chdir $path;
	#get all the files in that directory.
	@_=<*>;
  	for(@_){
    		if(-d $_){
      		#if the destination file is a directory, go recursion.
      		deleteFolder($_);
		}else
		{
      			unlink;
		}
  	}
	#Go up and del the destination directory.
	chdir "../";
	rmdir $path;
}

=head2 write_log

=cut
sub write_log
{
	my $content = shift;
	my $fh = IO::File->new(">>".$log_file) || die "Can not write info to log file $! \n";
	print $fh $content;
	$fh->close;
}
