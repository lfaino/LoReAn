#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Config::Std;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "

Synopsis:

gtf2gff3 --cfg gtf2gff3_MY_CONFIG.cfg gtf_file > gff3_file

gtf2gff3 --help # for a more detailed help message.

Description:

This script will convert GTF formatted files to valid GFF3 formatted
files.  It will map the column 3 (\"type\" column) to valid SO, but
because any non standard term may appear in that column in GTF files,
you may edit the config file to provide your own GTF feature to SO
mapping.  The script will also build gene models from exons, CDSs and
other features given in the GTF file.  It is currently tested on Ensemble
and Twinscan GTF, and it should work on any other files that follow the
same specification.  It does not work on GTF from the UCSC table browser
because those files use the same ID for gene and transcript, so it is
impossible to group multiple transcripts to a gene.  See the README that
came with the script for more info.

Options:

  --cfg

    Provide the filename for a config file.  See the configuration
    file provided with this script for format details.  Use this
    configuration file to modify the behavior of the script. If no
    config file is given it looks for ./gtf2gff3.cfg, ~/gtf2gff3.cfg
    or /etc/gtf2gff3.cfg in that order.

  --print_cfg

    Print a copy of the config variables as currently set (including
    those set by --cfg) to STDOUT.  This output can be used as the
    starting point for creating a custom config file.

  --help  Provide a more detailed help message.

";

my ($help, $cfg_file, $print_cfg);
my $opt_success;
eval{$opt_success = GetOptions('help'      => \$help,
			       'cfg=s'     => \$cfg_file,
                               'print_cfg' => \$print_cfg)};

handle_message('FATAL', 'invalid_options', $@) unless $opt_success;

if ($help) {
    print `perldoc $0`;
    exit;
}

my $file = shift;
if (! $print_cfg) {
    if (! defined $file) {
	handle_message('FATAL', 'missing_gtf_file');
    }
    elsif (! -e $file) {
	handle_message('FATAL', 'file_does_not_exist', $file);
    }
    elsif (! -r $file) {
	handle_message('FATAL', 'file_not_readable', $file);
    }
}

my $home = $ENV{HOME} || '.'; #Apache has no home, so avoid the error.

my @cfg_files = ($cfg_file,
		 './gtf2gff3.cfg',
		 "$home/gtf2gff3.cfg",
		 '/etc/gtf2gff3.cfg',
		 );

@cfg_files = grep {-e $_ if $_} @cfg_files;
$cfg_file = shift @cfg_files;
$cfg_file ||= '';

################################################################################
# LOAD CONFIGURATION OR USE DEFAULTS
################################################################################

my %config;
if ($cfg_file) {
    handle_message('FATAL', 'unreadable_file', $cfg_file) unless -r $cfg_file;
    eval {
	require Config::Std;
    };
    if ($@ =~ /Can\'t locate/) {
	handle_message('WARN', 'missing_dependency', 'Config::Std not installed. Using internal defaults.');
    }
    else {
	read_config($cfg_file => %config);
    }
}

# Use INPUT_FEATURE_MAP to map GTF feature types (column 3 in GTF) to valid SO types.
# Don't edit the SO terms below. Mapping must be many to one.  That means that
# exon_this and exon_that could both map to the SO term exon, but exon_this
# could not map to multiple SO terms.
#                               GTF Term     (maps to)       #SO Term
$config{INPUT_FEATURE_MAP} ||= {"3'-UTR"                   =>  'three_prime_UTR',
				'3UTR'                     =>  'three_prime_UTR',
				"5'-UTR"                   =>  'five_prime_UTR',
				'5UTR'                     =>  'five_prime_UTR',
				ARS                        =>  'ARS',
				binding_site               =>  'binding_site',
				BLASTN_HIT                 =>  'nucleotide_match',
				CDS                        =>  'CDS',
				CDS_motif                  =>  'nucleotide_motif',
				CDS_parts                  =>  'mRNA_region',
				centromere                 =>  'centromere',
				chromosome                 =>  'chromosome',
				conflict                   =>  'conflict',
				Contig                     =>  'contig',
				exon                       =>  'exon',
				five_prime_utr             =>  'five_prime_UTR',
				five_prime_UTR             =>  'five_prime_UTR',
				gene                       =>  'gene',
				insertion                  =>  'insertion',
				intron                     =>  'intron',
				LTR                        =>  'long_terminal_repeat',
				misc_feature               =>  'sequence_feature',
				misc_RNA                   =>  'transcript',
				mRNA                       =>  'mRNA',
				nc_primary_transcript      =>  'nc_primary_transcript',
				ncRNA                      =>  'ncRNA',
				nucleotide_match           =>  'nucleotide_match',
				polyA_signal               =>  'polyA_signal_sequence',
				polyA_site                 =>  'polyA_site',
				promoter                   =>  'promoter',
				pseudogene                 =>  'pseudogene',
				real_mRNA                  =>  'mRNA',
				region                     =>  'region',
				repeat_family              =>  'repeat_family',
				repeat_region              =>  'repeat_region',
				repeat_unit                =>  'repeat_region',
				rep_origin                 =>  'origin_of_replication',
				rRNA                       =>  'rRNA',
				Selenocysteine             =>  'selenocysteine',
				snoRNA                     =>  'snoRNA',
				snRNA                      =>  'snRNA',
				source                     =>  'sequence_feature',
				start_codon                =>  'start_codon',
				start_codon                =>  'start_codon',
				stop_codon                 =>  'stop_codon',
				stop_codon                 =>  'stop_codon',
				telomere                   =>  'telomere',
				three_prime_utr            =>  'three_prime_UTR',
				three_prime_UTR            =>  'three_prime_UTR',
				transcript_region          =>  'transcript_region',
				transcript                 =>  'transcript',
				transposable_element_gene  =>  'transposable_element',
				transposable_element       =>  'transposable_element',
				tRNA                       =>  'tRNA',
				tss                        =>  'TSS',
				TSS                        =>  'TSS',
				tts                        =>  'transcription_end_site',
				TTS                        =>  'transcription_end_site',
				UTR                        =>  'UTR',
                               };

our $INPUT_FEATURE_MAP = $config{INPUT_FEATURE_MAP};

# Maps tags used internally to tags found in your GTF file
# Don't edit the code tags.
# Note that the gene_id and transcript_id tags tell the script
# who the parents of a feature are.
#                           Code Tag  (comes from)  #GTF Tag
$config{GTF_ATTRB_MAP} ||= {gene_id		  => 'gene_id',
			    gene_name		  => 'gene_name',
			    trnsc_id		  => 'transcript_id',
			    trnsc_name		  => 'transcript_name'};

our $GTF_ATTRB_MAP  = $config{GTF_ATTRB_MAP};
# Maps tags used internally to output GFF3 attribute tags.
# Also, when LIMIT_ATTRB is set to 1 only these tags will be
# Output to the GFF3 attributes column.
#                            Code Tag      #GFF3 Tag
$config{GFF3_ATTRB_MAP} ||= {gene_id    => 'gene_id',
			     gene_name  => 'gene_name',
			     trnsc_id   => 'transcript_id',
			     trnsc_name => 'transcript_name',
			     id         => 'ID',
			     parent     => 'Parent',
			     name       => 'Name'};
our $GFF3_ATTRB_MAP = $config{GFF3_ATTRB_MAP};

# Limit the attribute tags printed to only those in the GFF3_ATTRB_MAP
$config{MISC}{LIMIT_ATTRB} ||= 0;
our $LIMIT_ATTRB     = $config{MISC}{LIMIT_ATTRB};
# A perl regexp that splits the attributes column into seperate attributes.
$config{MISC}{ATTRB_DELIMITER} ||= qr{\s*;\s*};
our $ATTRB_DELIMITER = $config{MISC}{ATTRB_DELIMITER};
# A perl regexp that captures the tag value pairs.
$config{MISC}{ATTRB_REGEX} ||= qr{^\s*(\S+)\s+(\"?[^\"]*\"?)\s*$};
our $ATTRB_REGEX = $config{MISC}{ATTRB_REGEX};
# If CDSs are annotated in the GTF file, are the start codons already included (1=yes 0=no)
$config{MISC}{START_IN_CDS} = defined $config{MISC}{START_IN_CDS} ? $config{MISC}{START_IN_CDS} : 1;
our $START_IN_CDS   = $config{MISC}{START_IN_CDS};
# If CDSs are annotated in the GTF file, are the stop codons already included (1=yes 0=no)
$config{MISC}{STOP_IN_CDS} = defined $config{MISC}{STOP_IN_CDS} ? $config{MISC}{STOP_IN_CDS} : 0;
our $STOP_IN_CDS    = $config{MISC}{STOP_IN_CDS};
# Use the following value (+ or -) for a features strand as the default if an invalid value is passed.
$config{MISC}{DEFAULT_STRAND} = defined $config{MISC}{DEFAULT_STRAND} ? $config{MISC}{DEFAULT_STRAND} : '+';
our $DEFAULT_STRAND = $config{MISC}{DEFAULT_STRAND};
################################################################################

print_config(\%config) if $print_cfg;

my @META_DATA;
my ($genes, $features) = parse_gtf($file);

$genes = build_genes($genes);

print_gff3($genes, $features);

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub parse_gtf {

    my $file = shift;
    open (my $IN, '<', $file) or
	handle_message('FATAL', 'cant_open_file_for_reading', $file);

    # %genes will have the following structure:
    # $genes->gene_id->trnsc_id->feature_type->feature
    #                |
    #                ->other->feature_type->feature
    #
    my (%genes, %features);

  LINE:
    while (my $line = <$IN>) {

	chomp $line;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;

	next LINE unless $line;

	#Handle full line comments, meta-comments and blank lines;
	if ($line =~ /^\#/) {
	    if ($line =~ /^\#\#gff-version\s+3/) {
		handle_message('FATAL', 'file_claims_to_be_gff3', $line);
	    }
	    $line =~ s/^\#+/\#/;
	    push @META_DATA, $line;
	    next LINE;
	}

	#Grab any comments
	my ($comment) = $line =~ s/(\#.*)$//;

	my @columns = split /\t/, $line;

	handle_message('FATAL', 'too_many_columns_for_gtf_data', $line)
	    unless scalar @columns == 9;

	my ($seqname, $source, $feature_type, $start, $end, $score,
	    $strand, $frame, $attrb_text) = @columns;

	handle_message('FATAL', 'invalid_start', $line)
	    unless $start =~ /^\d+$/;
	handle_message('FATAL', 'invalid_end', $line)
	    unless $end =~ /^\d+$/;
	handle_message('WARN', 'missing_attribute_text', $line)
	    unless $attrb_text;

	my $attributes = parse_attributes($attrb_text);
	if (! scalar keys %{$attributes}) {
	    handle_message('WARN', 'failed_to_parse_attributes', $line);
	}

	if (exists $INPUT_FEATURE_MAP->{$feature_type}) {
	    $feature_type = $INPUT_FEATURE_MAP->{$feature_type};
	}
	else {
	    handle_message('WARN', 'cant_map_feature_type', $feature_type);
	    $feature_type = 'sequence_feature';
	}

	#Note here that we're mapping between GTF/GFF and GFF3
	#nomeclature for the hash keys
	my $feature = {seq_id     => $seqname,
		       source     => $source,
		       type       => $feature_type,
		       start      => $start,
		       end        => $end,
		       score      => $score,
		       strand     => $strand,
		       phase      => $frame,
		       attributes => $attributes,
		       comment    => $comment};

	#Feature has a gene_id
	if (exists $attributes->{$GTF_ATTRB_MAP->{gene_id}}) {

	    my $gene_id = $attributes->{$GTF_ATTRB_MAP->{gene_id}}[0];

	    #Feature has a gene ID and transcript ID
	    if (exists $attributes->{$GTF_ATTRB_MAP->{trnsc_id}}) {
		my $trnsc_id = $attributes->{$GTF_ATTRB_MAP->{trnsc_id}}[0];
		# fix broken UCSC GTF!
		#if ($gene_id eq $trnsc_id) {
		#    ($gene_id = $trnsc_id) =~ s/\.\d+$//;
		#    $attributes->{gene_id}[0] = $gene_id;
		#}
		push @{$genes{$seqname}{$gene_id}{$trnsc_id}{$feature_type}},
		$feature;
	    }
	    # Feature is a gene
	    elsif ($feature->{type} eq 'gene') {
		# Do nothing
	    }
	    else {
		push @{$genes{$seqname}{$gene_id}{other}{$feature}},
		$feature;
	    }

	}
	else {
	    push @{$features{$feature_type}}, $feature;
	}
    }
    return (\%genes, \%features);
}

#-----------------------------------------------------------------------------

sub parse_attributes {

    my $attrb_text = shift;

    #Skip it if it's undefined or empty
    if (! defined $attrb_text || $attrb_text =~ /^\s+$/) {
	    handle_message('WARN', 'empty_attribute_text');
	return {};
    }

    #Split the attributes
    my @attrb_array = split $ATTRB_DELIMITER, $attrb_text;

    #Grab key, value pairs
    my %attributes;
  ATT:
    for my $attrb (@attrb_array) {
	my ($key, $value);

	if ($attrb =~ $ATTRB_REGEX) {
	    ($key, $value) = ($1, $2);
	    $value =~ s/\"//g if $value !~ /\s/;
	}
	elsif (! $attrb) {
	    handle_message('WARN', 'empty_attribute',
			   'Check $ATTRB_DELIMITER in the config file ' .
			   'to match: ' . $attrb_text);
	    next ATT;
	}
	else {
	    handle_message('WARN', 'attribute_regex_doesnt_match',
			   'Check $ATTRB_REGEX in the config file ' .
			   'to match: ' . $attrb);
	    next ATT;
	}
	push @{$attributes{$key}}, $value;
    }
    return \%attributes;
}

#-----------------------------------------------------------------------------

sub build_genes {

    my $genes = shift;

    my @genes;
    for my $seq_id (keys %{$genes}) {
	for my $gene_id (keys %{$genes->{$seq_id}}) {

	    my $gene = build_gene($genes->{$seq_id}{$gene_id});
	    push @genes, $gene;
	}
    }
    return \@genes;
}

#-----------------------------------------------------------------------------

sub build_gene {

    my $trnscs = shift;

    my @trnscs;
  TRN:
    for my $trnsc_id (keys %{$trnscs}) {

	if ($trnsc_id eq 'other') {
	    #Handle non-transcript gene features
	    handle_message('WARN', 'non-transcript_gene_feature',
			   'Non-transcript gene features'  .
			   'are not currently supported.'  .
			   'Please contact the author for' .
			   'support.');
	    next TRN;
	}
	my $features = $trnscs->{$trnsc_id};
	my $trnsc = build_trnsc($features);
	if ($trnsc) {
	    push @trnscs, $trnsc;
	}
	else {
	    handle_message('WARN', 'transcript_not_built', $trnsc_id);
	}
    }

    my $gene = validate_and_build_gene(\@trnscs);
    return $gene;
}

#-----------------------------------------------------------------------------

sub build_trnsc {

    my $features = shift;

    my $exons      = $features->{exon}	        || [];
    my $CDSs       = $features->{CDS}	        || [];
    my $start      = $features->{start_codon}	|| [];
    my $stop       = $features->{stop_codon}	|| [];
    my $five_UTRs  = $features->{five_prime_UTR}    || [];
    my $three_UTRs = $features->{three_prime_UTR}   || [];

    #We require at least CDSs or exons
    unless (scalar @{$CDSs} || scalar @{$exons}) {
	handle_message('WARN', 'must_have_exons_or_CDS_to_build_transcript',
		       join ", ", keys %{$features});
	return;
    }

    #Make start codons if they don't exist
    $start = process_start($exons, $CDSs, $start, $five_UTRs)
	if ! scalar @{$start};
    #Make stop codons if they don't exist
    $stop  = process_stop($exons, $CDSs, $stop, $three_UTRs)
	if ! scalar @{$stop};

    #Make UTRs, exons and/or CDSs if they don't exist
    ($five_UTRs,
     $exons,
     $CDSs,
     $three_UTRs) = process_exon_CDS_UTRs($exons, $CDSs,
					  $start, $stop,
					  $five_UTRs,
					  $three_UTRs);

    @{$features}{qw(exon five_prime_UTR start_codon CDS stop_codon
		    three_prime_UTR)} =
		    ($exons, $five_UTRs, $start, $CDSs,
		     $stop, $three_UTRs);

    #Build and validate the transcripts
    my $trnsc = validate_and_finish_trnsc($features);

    return $trnsc;
}

#-----------------------------------------------------------------------------

sub process_start {
    my ($exons, $CDSs, $start, $five_UTRs) = @_;

    my ($start_codon_start, $start_codon_end);

    my $strand;
    #Get the strand for exons or CDSs
    if (scalar @{$exons}) {
	$strand = strand($exons->[0]->{strand});
    }
    elsif (scalar @{$CDSs}) {
	$strand = strand($CDSs->[0]->{strand});
    }
    else {
	handle_message('WARN', 'need_exons_or_cdss_to_build_transcripts');
    }

    #Already have start
    if (scalar @{$start}) {
	return $start;
    }

    ################################################################################
    # Don't infer start codons from CDS unless backed by UTR or exon infered UTR
    ################################################################################

    #Build start from CDS
    elsif (scalar @{$CDSs} && scalar @{$exons}) {

	$CDSs = sort_features($CDSs, $strand);
	my $first_CDS = $CDSs->[0];

	$exons = sort_features($exons, $strand);
	my $first_exon = $exons->[0];

	if ($strand == 1) {
	    if ($first_exon->{start} < $first_CDS->{start}) {
		$start_codon_start = $first_CDS->{start};
		#Move start codon out of CDS if config says to
		if ($START_IN_CDS == 0) {
		    $start_codon_start -= 3;
		}
	    }
	}
	elsif ($strand == -1) {
	    if ($first_exon->{end} > $first_CDS->{end}) {
		$start_codon_start = $first_CDS->{end} - 2;
		#Move start codon out of CDS if config says to
		if ($START_IN_CDS == 0) {
		    $start_codon_start += 3;
		}
	    }
	}
	$start_codon_end   = $start_codon_start + 2
	    if defined $start_codon_start;
    }
    #Here we assume that if you have 5' UTR and CDS that the start codon
    #must be at the begining of the first CDS.  To be more rigorous we should
    #check the coordinates and be sure that the 5' UTR and CDS are contiguous.
    elsif (scalar @{$CDSs} && scalar @{$five_UTRs}) {

	$CDSs = sort_features($CDSs, $strand);
	my $first_CDS = $CDSs->[0];

	if ($strand == 1) {
	    $start_codon_start = $first_CDS->{start};
	    #Move start codon out of CDS if config says to
	    if ($START_IN_CDS == 0) {
		$start_codon_start -= 3;
	    }
	}
	elsif ($strand == -1) {
	    $start_codon_start = $first_CDS->{end} - 2;
	    #Move start codon out of CDS if config says to
	    if ($START_IN_CDS == 0) {
		$start_codon_start += 3;
	    }
	}
	$start_codon_end   = $start_codon_start + 2;
    }
    #Build start from UTRs - I haven't seen an example of this yet. Maybe I should
    #only create a start codon if it falls within an exon or CDS.
    elsif (scalar @{$five_UTRs}) {
	handle_message('FATAL', 'untested_code', 'The subroutine '       .
		       'process_start has run into data it doesnt know ' .
		       'how to handle. Please contact the aurthor for '  .
		       'support if you see this message');
	$five_UTRs = sort_features($five_UTRs, $strand);
	my $last_five_UTR = $five_UTRs->[-1];
	$start_codon_start = $strand == 1 ? $last_five_UTR->{start} + 1 :
	    $last_five_UTR->{end} - 3;
	$start_codon_end   = $start_codon_start + 2;
    }
    #No CDSs or UTRs - let's pretend it's non-coding
    elsif (scalar @{$exons}) {
	return [];
    }
    else {
	return [];
	# We don't really need to die here some annotations have
	# stop_codon, but not start_codon
    }

    return [] if (! defined $start_codon_start ||
		  ! defined $start_codon_end);

    $start = [{start => $start_codon_start,
	       end   => $start_codon_end,
	       type  => 'start_codon',
	       score => '.',
	       phase => '0'}];

    return $start;
}

#-----------------------------------------------------------------------------

sub process_stop {
    my ($exons, $CDSs, $stop, $three_UTRs) = @_;

    my ($stop_codon_start, $stop_codon_end);

    #Get the strand
    my $strand;
    if (scalar @{$exons}) {
	$strand = strand($exons->[0]->{strand});
    }
    elsif (scalar @{$CDSs}) {
	$strand = strand($CDSs->[0]->{strand});
    }
    else {
	handle_message('WARN', 'need_exons_or_cdss_to_build_transcripts');
    }

    #If we already have a stop then return it.
    if (scalar @{$stop}) {
	return $stop;
    }

    ################################################################################
    # Don't infer stop codons from CDS unless backed by UTR or exon infered UTR
    ################################################################################

    #Build stop from CDSs
    elsif (scalar @{$CDSs} && scalar @{$exons}) {

	$CDSs = sort_features($CDSs, $strand);
	my $last_CDS = $CDSs->[-1];

	$exons = sort_features($exons, $strand);
	my $last_exon = $exons->[-1];

	if ($strand == 1) {
	    if ($last_exon->{end} > $last_CDS->{end}) {
		$stop_codon_start = $last_CDS->{end} - 2;
		#Move stop codon out of CDS is config says to
		if ($STOP_IN_CDS == 0) {
		    $stop_codon_start += 3;
		}
	    }
	}
	elsif ($strand == -1) {
	    if ($last_exon->{start} < $last_CDS->{start}) {
		$stop_codon_start = $last_CDS->{start};
		#Move stop codon out of CDS is config says to
		if ($STOP_IN_CDS == 0) {
		    $stop_codon_start -= 3;
		}
	    }
	}

	$stop_codon_end   = $stop_codon_start + 2
	    if defined $stop_codon_start;
    }
    elsif (scalar @{$CDSs} && scalar @{$three_UTRs}) {

	$CDSs = sort_features($CDSs, $strand);
	my $last_CDS = $CDSs->[-1];

	if ($strand == 1) {
	    $stop_codon_start = $last_CDS->{end} - 2;
	    #Move stop codon out of CDS is config says to
	    if ($STOP_IN_CDS == 0) {
		$stop_codon_start += 3;
	    }
	}
	elsif ($strand == -1) {
	    $stop_codon_start = $last_CDS->{start};
	    #Move stop codon out of CDS is config says to
	    if ($STOP_IN_CDS == 0) {
		$stop_codon_start -= 3;
	    }
	}
	$stop_codon_end   = $stop_codon_start + 2;
    }
    #Build stop from UTRs
    elsif (scalar @{$three_UTRs}) {
	handle_message('FATAL', 'untested_code', 'The subroutine '       .
		       'process_stop has run into data it doesn\'t know ' .
		       'how to handle. Please contact the aurthor for '  .
		       'support if you see this message');
	my $strand = strand($three_UTRs->[0]);
	$three_UTRs = sort_features($three_UTRs, $strand);
	my $first_three_UTR = $three_UTRs->[0];
	$stop_codon_start = $strand == 1 ? $first_three_UTR->{start} - 3 :
	    $first_three_UTR->{end} + 1;

	$stop_codon_end   = $stop_codon_start + 2;
    }
    elsif (scalar @{$exons}) {
	#Treating this as a non-coding transcript
	return [];
    }
    else {
	# We don't really need to die here.  Some features may have
	# start_codon, but not stop_codon
	handle_message('WARN', 'cant_infer_stop');
	return [];
	#die "FATAL : Invalid feature set : process_stop\n";
    }

    return [] if (! defined $stop_codon_start ||
		  ! defined $stop_codon_end);

    $stop = [{start => $stop_codon_start,
	      end   => $stop_codon_end,
	      type  => 'stop_codon',
	      score => '.',
	      phase => '0'}];

    return $stop;
}

#-----------------------------------------------------------------------------

sub process_exon_CDS_UTRs {
    my ($exons, $CDSs, $start, $stop,
	$five_UTRs, $three_UTRs) = @_;

    #Check what features we already have so we don't rebuild them
    my $have_five_UTRs++  if scalar @{$five_UTRs};
    my $have_exons++      if scalar @{$exons};
    my $have_start++      if scalar @{$start};
    my $have_CDSs++       if scalar @{$CDSs};
    my $have_stop++       if scalar @{$stop};
    my $have_three_UTRs++ if scalar @{$three_UTRs};

    #If CDSs already exist make sure that they include start and stop codons
    if (scalar @{$CDSs}) {
	include_terminal_codons($CDSs, $start, $stop);
    }

    #If we already have everything, then return
    if ($have_exons && $have_CDSs && ($have_five_UTRs || $have_three_UTRs)) {
	return ($five_UTRs, $exons, $CDSs, $three_UTRs);
    }
    #Build CDSs && UTRs
    elsif ($have_exons && ($have_start || $have_stop || $have_CDSs ||
			   $have_three_UTRs || $have_five_UTRs)) {
	# Make sure that evaluate exons can handle only a start OR a stop!!!
	for my $exon (@{$exons}) {
	    my ($five_UTR, $CDS, $three_UTR) =
		evaluate_exon($exon, $start->[0], $stop->[0]);
	    push @{$five_UTRs}, $five_UTR
		if scalar keys %{$five_UTR} && ! $have_five_UTRs;
	    push @{$CDSs}, $CDS if scalar keys %{$CDS} && ! $have_CDSs;
	    push @{$three_UTRs}, $three_UTR
		if scalar keys %{$three_UTR} && ! $have_three_UTRs;
	}
    }
    #Build exons
    elsif (!$have_exons && $have_CDSs){
	$exons = build_exons($five_UTRs, $start, $CDSs, $stop, $three_UTRs);
    }
    #Treat as non-coding even if we have a start or stop but not both
    elsif ($have_exons && ! $have_CDSs && ! $have_five_UTRs &&
	   ! $have_three_UTRs && ! $have_start && ! $have_stop) {
	#Treating this as a non_coding transcript
	return ($five_UTRs, $exons, $CDSs, $three_UTRs);
    }
    else {
	handle_message('FATAL', 'invalid_feature_set',
		       'The subroutine process_exon_CDS_UTRs was given ' .
		       'a set of feature that it can\'t work with.');
    }

    return ($five_UTRs, $exons, $CDSs, $three_UTRs);
}

#-----------------------------------------------------------------------------

sub include_terminal_codons {
    my ($CDSs, $start, $stop) = @_;

    my $strand = strand($CDSs->[0]{strand});

    $CDSs = sort_features($CDSs, $strand);
    my $first_CDS = $CDSs->[0];
    my $last_CDS  = $CDSs->[-1];

    if ($strand == 1) {
	if (scalar @{$start}) {
	    $first_CDS->{start} = $start->[0]{start}
	    if $start->[0]{end} + 1 == $first_CDS->{start};
	}
	#START_IN_CDS was considered when start was built, so it could introduce errors
	#to consider it again here
	#		elsif (defined $START_IN_CDS) {
	#			$first_CDS->{start} -= 3 if ! $START_IN_CDS;
	#		}

	if (scalar @{$stop}) {
	    $last_CDS->{end} = $stop->[0]{end}
	    if $stop->[0]{start} - 1 == $last_CDS->{end};
	}

	#STOP_IN_CDS was considered when stop was built, so it could introduce errors
	#to consider it again here
	#		elsif (defined $STOP_IN_CDS) {
	#			$last_CDS->{end} += 3 if ! $STOP_IN_CDS;
	#		}
    }
    elsif ($strand == -1) {
	if (scalar @{$start}) {
	    $first_CDS->{end} = $start->[0]{end}
	    if $start->[0]{start} - 1 == $first_CDS->{end};
	}
	#START_IN_CDS was considered when start was built, so it could introduce errors
	#to consider it again here
	#		elsif (defined $START_IN_CDS) {
	#			$first_CDS->{end} += 3 if ! $START_IN_CDS;
	#		}

	if (scalar @{$stop}) {
	    $last_CDS->{start} = $stop->[0]{start}
	    if $stop->[0]{end} + 1 == $last_CDS->{start};
	}
	#STOP_IN_CDS was considered when stop was built, so it could introduce errors
	#to consider it again here
	#		elsif (defined $STOP_IN_CDS) {
	#			$last_CDS->{start} -= 3 if ! $STOP_IN_CDS;
	#               }
    }
}

#-----------------------------------------------------------------------------

sub evaluate_exon {
    my ($exon, $start, $stop) = @_;

    my $strand = strand($exon);

    my %five_UTR;
    my %CDS;
    my %three_UTR;

    #######################################
    # Allow success if missing either start or stop
    #######################################

    if ($strand == 1) {
	#Exon is fully 5' UTR
	if (defined $start->{start} &&
	    $exon->{end} <= $start->{start}) {
	    $five_UTR{start} = $exon->{start};
	    $five_UTR{end}   = $exon->{end};
	}
	#Exon stradles start codon
	elsif (defined $start->{start}           &&
	       $exon->{start} <  $start->{start} &&
	       $exon->{end}   >  $start->{start}     ) {
	    $five_UTR{start} = $exon->{start};
	    $five_UTR{end}   = $start->{start} - 1;
	    if (defined $stop->{start}) {
		$CDS{start}      = $start->{start};
		$CDS{end}        = $exon->{end};
		$CDS{phase}      = 0;
	    }
	}
	#Exon is fully CDS
	elsif (defined $start->{start}           &&
	       defined $stop->{end}              &&
	       $exon->{start} >= $start->{start} &&
	       $exon->{end}   <= $stop->{end}        ) {
	    $CDS{start}      = $exon->{start};
	    $CDS{end}        = $exon->{end};
	    #Set phase if this is first CDS
	    if ($exon->{start} == $start->{end} + 1  ||
		$exon->{start} == $start->{start}    ) {
		$CDS{phase} = 0;
	    }
	}
	#Exon stradles stop codon
	if (defined $stop->{end}              &&
	    $exon->{start} <  $stop->{end}    &&
	    $exon->{end}   >  $stop->{end}        ) {
	  $three_UTR{start} = $stop->{end} + 1;
	  $three_UTR{end}   = $exon->{end};
	  if (defined $start->{start}) {
	    $CDS{start}       = $exon->{start};
	    $CDS{end}         = $stop->{end};
	  }
	}
	#Exon is fully 3' UTR
	elsif (defined $stop->{end} &&
	       $exon->{start} >= $stop->{end}) {
	    $three_UTR{start} = $exon->{start};
	    $three_UTR{end}   = $exon->{end};
	}
	else {
	    # die "Fatal error in evaluate_exon\n";
	}
    }
    else {
	#Exon is fully 5' UTR
	if (defined $start->{end} &&
	    $exon->{start} >= $start->{end}) {
	    $five_UTR{start} = $exon->{start};
	    $five_UTR{end}   = $exon->{end};
	}
	#Exon stradles start codon
	elsif (defined $start->{end}            &&
	       $exon->{start} < $start->{end}   &&
	       $exon->{end}   > $start->{end}      ) {
	    $five_UTR{start} = $start->{end} + 1;
	    $five_UTR{end}   = $exon->{end};
	    if (defined $stop->{start}) {
		$CDS{start} = $exon->{start};
		$CDS{end}   = $start->{end};
	    }
	    $CDS{phase} = 0;
	}
	#Exon is fully CDS
	elsif (defined $start->{end}              &&
	       defined $stop->{start}             &&
	       $exon->{end}  <= $start->{end}  &&
	       $exon->{start} >= $stop->{start}    ) {
	    $CDS{start} = $exon->{start};
	    $CDS{end} = $exon->{end};
	    #Set phase 0 if this is first CDS
	    if ($exon->{end} == $start->{start} - 1  ||
		$exon->{end} == $start->{end}    ) {
		$CDS{phase} = 0;
	    }
	}
	#Exon stradles stop codon
	if (defined $stop->{start}           &&
	       $exon->{end}   > $stop->{start}  &&
	       $exon->{start} < $stop->{start}     ) {
	    $three_UTR{start} = $exon->{start};
	    $three_UTR{end}   = $stop->{start} - 1;
	    if (defined $start->{end}) {
		$CDS{start} = $stop->{start};
		$CDS{end}   = $exon->{end};
	    }
	}
	#Exon is fully 3' UTR
	elsif (defined $stop->{start} &&
	       $exon->{end} <= $stop->{start}) {
	    $three_UTR{start} = $exon->{start};
	    $three_UTR{end}   = $exon->{end};
	}
	else {
	    # die "Fatal error in evaluate_exon\n";
	}
    }
    return (\%five_UTR, \%CDS, \%three_UTR);
}

#-----------------------------------------------------------------------------

sub build_exons {
    my ($five_UTRs, $start, $CDSs, $stop, $three_UTRs) = @_;

    my $strand = $CDSs->[0]{strand};

    my @exons;

    #Make an exon for every 5' UTR
    for my $five_UTR (@{$five_UTRs}) {
	push @exons, {start => $five_UTR->{start},
		      end   => $five_UTR->{end},
		      type  => 'exon',
		      score      => '.',
		      strand     => $strand,
		      phase      => '.'};
    }

    #Make and exon for every CDS
    $CDSs = sort_features($CDSs, $strand);
    for my $CDS (@{$CDSs}) {
	my $exon = {start => $CDS->{start},
		    end   => $CDS->{end},
		    type  => 'exon',
		    score      => '.',
		    strand     => $strand,
		    phase      => '.'};

    push @exons, $exon;
    }

    #Make an exon for every 3' UTR
    for my $three_UTR (@{$three_UTRs}) {
	push @exons, {start 	 => $three_UTR->{start},
		      end   	 => $three_UTR->{end},
		      type  	 => 'exon',
		      score      => '.',
		      strand     => $strand,
		      phase      => '.'};
    }

    # Merge any contiguous exons (i.e. from UTR & CDS neighbors)
    my $exons = sort_features(\@exons, 1);
    my @merged_exons = shift @{$exons};
    while (my $exon = shift @{$exons}) {
	if ($exon->{start} <= $merged_exons[-1]{end} + 1) {
	    $merged_exons[-1]{end} = $exon->{end};
	}
	else {
	    push @merged_exons, $exon;
	}
    }
    return \@merged_exons;
}

#-----------------------------------------------------------------------------

sub validate_and_finish_trnsc {
    my $features = shift;

    my ($seq_id, $source, $strand, $gene_id, $gene_name,
	$trnsc_id, $trnsc_name);

    #Find default parameters from either exons or CDSs
  TYPE:
    for my $type ( qw(exon CDS) ) {
	for my $feature (@{$features->{$type}}) {
	    $seq_id     ||= $feature->{seq_id};
	    $source     ||= $feature->{source};
	    $strand     ||= $feature->{strand};
	    $gene_id    ||=
		$feature->{attributes}{$GTF_ATTRB_MAP->{gene_id}};
	    $gene_name  ||=
		$feature->{attributes}{$GTF_ATTRB_MAP->{gene_name}};
	    $trnsc_id   ||=
		$feature->{attributes}{$GTF_ATTRB_MAP->{trnsc_id}};
	    $trnsc_name ||=
		$feature->{attributes}{$GTF_ATTRB_MAP->{trnsc_name}};
	    last TYPE if ! grep {! defined $_} ($seq_id, $source,
						$strand, $gene_name,
						$trnsc_id, $trnsc_name);
	}
    }


    #Flag if we have any coding features (CDS, start_codon, stop_codon)
    my $coding_flag;
    #Min and max for transcript boundaries
    my ($min, $max);
    for my $feature_type (keys %{$features}) {
	my $count;

	#For keeping track of the phase for the next CDS.
	my $next_phase = '.';

	#Sort the features
	$features->{$feature_type} =
	    sort_features($features->{$feature_type},
			  $strand);
	for my $feature (@{$features->{$feature_type}}) {

	    #min and max to calculate transcript boundaries
	    $min = ! defined $min        ?
		$feature->{start}        :
		$min > $feature->{start} ?
		$feature->{start}        :
		$min;

	    $max = ! defined $max        ?
		$feature->{end}          :
		$max < $feature->{end}   ?
		$feature->{end}          :
		$max;

	    #Set the flag if we see indications of coding features
	    if (grep {$feature_type eq $_} qw(CDS start_codon
					      stop_codon)) {
		$coding_flag++ if  $feature->{start};
	    }

	    #Calculate CDS phases.
	    if ($feature_type eq 'CDS') {
		($feature, $next_phase) = CDS_phase($feature, $next_phase);
	    }

	    #Set parameters to defaults for all features that
	    #don't already have them set
	    if (! defined $feature->{seq_id}) {
		$feature->{seq_id} = $seq_id;
	    }
	    elsif ($feature->{seq_id} ne $seq_id) {
		handle_message('WARN', 'transcript_seq_id_conflict',
			      "$trnsc_id->[0] ($feature->{seq_id}) ne $feature->{attributes}{id}[0] ($seq_id)");
	    }
	    if (! defined $feature->{source}) {
		$feature->{source} = $source;
	    }
	    elsif ($feature->{source} ne $source) {
		handle_message('WARN', 'transcript_source_conflict',
			      "$trnsc_id->[0] ($feature->{source}) ne $feature->{attributes}{id}[0] ($source)");
	    }
	    if (! defined $feature->{type}) {
		$feature->{type} = $feature_type;
	    }
	    elsif ($feature->{type} ne $feature_type) {
		handle_message('WARN', 'transcript_type_conflict',
			      "$trnsc_id->[0] ($feature->{type}) ne $feature->{attributes}{id}[0] ($feature_type)");
	    }
	    if (! defined $feature->{strand}) {
		$feature->{strand} = $strand;
	    }
	    elsif ($feature->{strand} ne $strand) {
		handle_message('WARN', 'transcript_strand_conflict',
			      "$trnsc_id->[0] ($feature->{strand}) ne $feature->{attributes}{id}[0] ($strand)");
	    }
	    if (! defined $feature->{score}) {
		$feature->{score} = '.';
	    }
	    if (! defined $feature->{phase}) {
		$feature->{phase} = '.';
	    }
	    #Set attributes
	    $feature->{attributes}{parent} = $trnsc_id;
	    $feature->{attributes}{id}     = ["$feature_type:" .
					      $trnsc_id->[0]   .
					      ":" . ++$count];
	}
    }

    my $trnsc_type = $coding_flag ? 'mRNA' : 'transcript';

    my $attributes = {parent      => $gene_id,
		      parent_name => $gene_name,
		      id          => $trnsc_id,
		      name        => $trnsc_name};

    my $trnsc = {seq_id     => $seq_id,
		 source     => $source,
		 type       => $trnsc_type,
		 start      => $min,
		 end        => $max,
		 score      => '.',
		 strand     => $strand,
		 phase      => '.',
		 attributes => $attributes,
		 features   => $features};

    return $trnsc;
}

#-----------------------------------------------------------------------------

sub CDS_phase {
    my ($feature, $next_phase) = @_;

    #If phase isn't already valid assign it
    if (! defined $feature->{phase} ||
	$feature->{phase} !~ /^0|1|2$/) {
	$feature->{phase} = $next_phase;
    }
    #If phase is valid, calculate the next phase
    if ($feature->{phase} =~ /^0|1|2$/) {
	my $length = ($feature->{end} -
		      $feature->{start}) + 1;
	# my $hang_3 = $length % 3; # 3' overhang
	# my $hang_5 = 3 - $hang_3; # 5' overhang
	#
	# #The next phase is equal to this phase
	# #plus the modulus 3 of the length wrapped
	# #at 2.
	# $next_phase = $feature->{phase} + $hang_5;
	# $next_phase -= 3 if $next_phase > 2;

	# This was update 5/24/10 in response to an
	# e-mail from Leighton Prichard regarding
	# errors in the GFF3 spec.  The code above
	# calculates the phase correctly, but the
	# formula suggested by Leighton is cleaner.

	$next_phase = ($feature->{phase} - $length) % 3;

    }

    return ($feature, $next_phase);
}

#-----------------------------------------------------------------------------

sub validate_and_build_gene {
    my $trnscs = shift;

    #Get parameter defaults
    my $seq_id     = $trnscs->[0]{seq_id};
    my $source     = $trnscs->[0]{source};
    my $strand     = $trnscs->[0]{strand};
    my $gene_id    = $trnscs->[0]{attributes}{parent};
    my $gene_name  = $trnscs->[0]{attributes}{parent_name};

    my %sources;
    #Get gene boundaries and check all transcripts for agreement with
    #parameter defaults
    my ($min, $max);
    for my $trnsc (@{$trnscs}) {
	$min = ! defined $min      ?
	    $trnsc->{start}        :
	    $min > $trnsc->{start} ?
	    $trnsc->{start}        :
	    $min;

	$max = ! defined $max      ?
	    $trnsc->{end}          :
	    $max < $trnsc->{end}   ?
	    $trnsc->{end}          :
	    $max;

	# Make sure feature seqid matches gene seqid
	handle_message('WARN', 'gene_seq_id_conflict',
		       "$gene_id->[0] ($seq_id) ne $trnsc->{attributes}{id}[0] ($trnsc->{seq_id})")
	    if $seq_id ne $trnsc->{seq_id};

	# Make sure feature source matches gene source
	handle_message('WARN', 'gene_source_conflict',
		       "$gene_id->[0] ($source) ne $trnsc->{attributes}{id}[0] ($trnsc->{source})")
	    if $source ne $trnsc->{source};
	$sources{$trnsc->{source}}++;

	# Make sure feature strand matches gene strand
	handle_message('WARN', 'gene_strand_conflict',
		       "$gene_id->[0] ($strand) ne $trnsc->{attributes}{id}[0] ($trnsc->{strand})")
	    if $strand ne $trnsc->{strand};

	# Make sure feature gene_id matches gene gene_id
	handle_message('WARN', 'gene_id_conflict',
		       "$gene_id->[0] (if $gene_id->[0] ne $trnsc->{attributes}{id}[0] (Parent: $trnsc->{attributes}{parent}[0])")
	    if $gene_id->[0] ne $trnsc->{attributes}{parent}[0];

    }
    my $attributes = {id   => $gene_id,
		      name => $gene_name || $gene_id};

    # If multiple sources exist for this gene default to '.'
    $source = scalar keys %sources > 1 ? '.' : $source;

    my $gene = {seq_id     => $seq_id,
		source     => $source,
		type       => 'gene',
		start      => $min,
		end        => $max,
		score      => '.',
		strand     => $strand,
		phase      => '.',
		attributes => $attributes,
		trnscs     => $trnscs};

    return $gene;
}

#-----------------------------------------------------------------------------

sub print_gff3 {
    my ($genes, $features) = @_;

    #GFF3 Header here
    print "##gff-version 3\n";

    if (scalar @META_DATA) {
      print join "\n", @META_DATA;
      print "\n";
    }

    for my $gene (@{$genes}) {
	print_gene($gene);
    }

    print_features($features) if scalar keys %{$features};
}

#-----------------------------------------------------------------------------

sub print_gene {
    my $gene = shift;

    my $attrb_text = make_attribute_text($gene->{attributes});
    #$attrb_text   .= ";Name=$gene->{attributes}{id}[0]";
    print join "\t", ($gene->{seq_id},
		      $gene->{source},
		      $gene->{type},
		      $gene->{start},
		      $gene->{end},
		      $gene->{score},
		      $gene->{strand},
		      $gene->{phase},
		      $attrb_text,
		      );
    print " " . $gene->{comment} if $gene->{comment};
    print "\n";

    for my $trnsc (@{$gene->{trnscs}}) {
	print_trnsc($trnsc);
    }
}

#-----------------------------------------------------------------------------

sub print_trnsc {
    my $trnsc = shift;

    my $attrb_text = make_attribute_text($trnsc->{attributes});
    print join "\t", ($trnsc->{seq_id},
		      $trnsc->{source},
		      $trnsc->{type},
		      $trnsc->{start},
		      $trnsc->{end},
		      $trnsc->{score},
		      $trnsc->{strand},
		      $trnsc->{phase},
		      $attrb_text,
		      );
    print " " . $trnsc->{comment} if $trnsc->{comment};
    print "\n";

    print_features($trnsc->{features});
}

#-----------------------------------------------------------------------------

sub print_features {
    my $features = shift;

    my @sorted_feature_types = sort_feature_types($features);
    for my $feature_type (@sorted_feature_types) {
	$features->{$feature_type} =
	    sort_features($features->{$feature_type}, 1);
	for my $feature (@{$features->{$feature_type}}) {
	    my $attrb_text =
		make_attribute_text($feature->{attributes});
	    print join "\t", ($feature->{seq_id},
			      $feature->{source},
			      $feature->{type},
			      $feature->{start},
			      $feature->{end},
			      $feature->{score},
			      $feature->{strand},
			      $feature->{phase},
			      $attrb_text,
			      );
	    #print " " . $feature->{comment} if
	    #$feature->{comment};
	    print "\n";
	    print '';
	}
    }
}

#-----------------------------------------------------------------------------

sub make_attribute_text {
    my $attributes = shift;

    #Only print the attributes listed in $LIMIT_ATTRB
    my @tags = $LIMIT_ATTRB ?
	grep {exists $GFF3_ATTRB_MAP->{$_}} keys %{$attributes} :
	keys %{$attributes};

    my %order = (id     => 'A',
		 name   => 'B',
		 parent => 'C',
		 );

    my @pairs;
    for my $tag (sort {($order{$a} || $a) cmp ($order{$b} || $b)} @tags) {
	next if $tag eq 'gene_id';
	next if $tag eq 'transcript_id';
	next unless $attributes->{$tag};
	my $tag_text = $GFF3_ATTRB_MAP->{$tag} || $tag;
	my $value_text = join ',', @{$attributes->{$tag}};
	my $pair_text = join '=', ($tag_text, $value_text);
	push @pairs, $pair_text;
    }
    my $attrb_text = join ';', @pairs;
    return $attrb_text;
}

#-----------------------------------------------------------------------------

sub sort_features {
    my ($features, $strand) = @_;

    $strand = strand($strand) || strand($features->[0]);

    #Make sure we get the array ref that we wanted
    if (ref $features eq 'ARRAY') {
	#If we get an array ref with more than one element - sort it
	if (scalar @{$features} > 1) {
	    #Sort + strand features
	    if ($strand == 1) {
		my @sorted_features = sort {$a->{start} <=> $b->{start}}
		@{$features};
		$features = \@sorted_features;
		return $features;
	    }
	    #Sort - strand features
	    elsif ($strand == -1) {
		my @sorted_features = sort {$b->{end} <=> $a->{end}}
		@{$features};
		$features = \@sorted_features;
		return $features;
	    }
	}
	#Empty or 1 length array just return it
	else {
	    return $features;
	}
    }
    #Don't allow someone to misuse sort_features - require array ref
    else {
	handle_message('FATAL', 'array_reference_required',
		       'The subroutine sort_features requires ' .
		       'an array reference, but was given a '   .
		       ref $features);
    }
}

#-----------------------------------------------------------------------------

sub sort_feature_types {
    my ($feature_types) = @_;

    #Make sure we get the hash ref that we wanted
    if (ref $feature_types eq 'HASH') {
	my @types = keys %{$feature_types};
	@types = grep{scalar @{$feature_types->{$_}} > 0} @types;
	my $strand;
	for my $type (@types) {
	    if (defined $feature_types->{$type}[0]{strand}) {
		$strand = $feature_types->{$type}[0]{strand};
		last;
	    }
	}
	$strand = strand($strand);
	handle_message('FATAL', 'cant_determine_strand',
		       'The sort_feature subroutine can not ' .
		       'determine the strand of the feature')
	    if ! defined $strand;

	#If we get an hash ref with more than one element - sort it
	if (scalar keys %{$feature_types} > 1) {
	    #Sort + strand features
	    if ($strand == 1) {
		my @sorted_types =
		    sort {
			min_feature($feature_types->{$a}, 'start') <=>
			    min_feature($feature_types->{$b}, 'start')
			    ||
			    min_feature($feature_types->{$a}, 'end') <=>
			    min_feature($feature_types->{$b}, 'end')
			} @types;
		return @sorted_types;
	    }
	    #Sort - strand features
	    elsif ($strand == -1) {
		my @sorted_types =
		    sort {
			max_feature($feature_types->{$b}, 'start') <=>
			    max_feature($feature_types->{$a}, 'start')
			    ||
			    max_feature($feature_types->{$b}, 'end') <=>
			    max_feature($feature_types->{$a}, 'end')
			} @types;
		return @sorted_types;
	    }
	}
	#Empty or 1 length hash just return it
	else {
	    return keys %{$feature_types};
	}
    }
    #Don't allow someone to misuse sort_features - require array ref
    else {
	handle_message('FATAL', 'hash_reference_required',
		       'The subroutine sort_feature_types '       .
		       'requires a hash reference but was given ' .
		       ref $features);

    }
}

#-----------------------------------------------------------------------------

sub strand {
    my $feature = shift;

    my $strand = ref $feature ? $feature->{strand} : $feature;

    #If strand is undefined or invalid...
    if (! defined $strand || $strand !~ /\+|-|1|-1/) {
	# ...and if we have a default value then use it...
	if ($DEFAULT_STRAND) {
	    $strand ||= 'undef';
	    handle_message('WARN', 'setting_invalid_strand_to_default_value',
			   "($strand -> $DEFAULT_STRAND)");
	    $strand = $DEFAULT_STRAND;
	}
	else {
	    # ...otherwise allow undefined to pass through..
	    return $strand if ! defined $strand;
	}
    }

    if ($strand eq '+') {
	return 1;
    }
    elsif ($strand eq '-') {
	return -1;
    }
    elsif ($strand == 1 || $strand == -1) {
	return $strand;
    }
    # ...finally if all else fails die
    handle_message('FATAL', 'invalid_value_for_strand',
		   "Set \$DEFAULT_STRAND in config file to handle ($strand)");
}

#-----------------------------------------------------------------------------

sub min_feature {
    my ($features, $terminus) = @_;

    my $min;

    for my $feature (@{$features}) {
	$min ||= $feature->{$terminus};
	$min = $feature->{$terminus} if $min > $feature->{$terminus};
    }
    return $min;
}

#-----------------------------------------------------------------------------

sub max_feature {
    my ($features, $terminus) = @_;

    my $max;

    for my $feature (@{$features}) {
	$max ||= $feature->{$terminus};
	$max = $feature->{$terminus} if $max < $feature->{$terminus};
    }
    return $max;
}

#-----------------------------------------------------------------------------

sub handle_message {

    my ($level, $code, @messages) = @_;

    $level ||= 'FATAL';
    $code  ||= 'unknown_warning';
    my $message = join ' ', @messages;
    $message ||= 'No additional message';

    $message = join ' : ', ($level, $code, $message);
    chomp $message;
    $message .= "\n";

    if ($level eq 'FATAL') {
	print STDERR $message;
	die;
    }
    else {
	print STDERR $message;
    }
}

#-----------------------------------------------------------------------------

sub print_config {

    my $config = shift;

    my $config_text = "#This config file allows the user to customize the gtf2gff3
#converter.

[INPUT_FEATURE_MAP]
#Use INPUT_FEATURE_MAP to map your GTF feature types (column 3 in GTF) to valid SO types.
#Don't edit the SO terms below.
#Mapping must be many to one.  That means that exon_this and exon_that could both
#map to the SO term exon, but exon_this could not map to multiple SO terms.

#GTF Term                  #SO Term
";

    for my $key (sort keys %{$config{INPUT_FEATURE_MAP}}) {
        my $value = $config{INPUT_FEATURE_MAP}{$key};
        $config_text .= "$key = $value\n";
    }   

    $config_text .= "\n[GTF_ATTRB_MAP]
#Maps tags used internally to tags found in your GTF file
#Don't edit the code tags.
#Note that the gene_id and transcript_id tags tell the script
#who the parents of a feature are.

#Code Tag    #GTF Tag
";

    for my $key (sort keys %{$config{GTF_ATTRB_MAP}}) {
        my $value = $config{GTF_ATTRB_MAP}{$key};
        $config_text .= "$key = $value\n";
    }   

    $config_text .= "\n[GFF3_ATTRB_MAP]
#Maps tags used internally to output GFF3 attribute tags.
#Also, when LIMIT_ATTRB is set to 1 only these tags will be
#Output to the GFF3 attributes column.

#Code Tag  #GFF3 Tag
";

    for my $key (sort keys %{$config{GFF3_ATTRB_MAP}}) {
        my $value = $config{GFF3_ATTRB_MAP}{$key};
        $config_text .= "$key = $value\n";
    }   

    $config_text .= "\n[MISC]
# Limit the attribute tags printed to only those in the GFF3_ATTRB_MAP
LIMIT_ATTRB     = $config{MISC}{LIMIT_ATTRB}
#A perl regexp that splits the attributes column into seperate attributes.
ATTRB_DELIMITER = $config{MISC}{ATTRB_DELIMITER}
#A perl regexp that captures the tag value pairs.
ATTRB_REGEX     = $config{MISC}{ATTRB_REGEX}
#If CDSs are annotated in the GTF file, are the start codons already included (1=yes 0=no)
START_IN_CDS    = $config{MISC}{START_IN_CDS}
#If CDSs are annotated in the GTF file, are the stop codons already included (1=yes 0=no)
STOP_IN_CDS     = $config{MISC}{STOP_IN_CDS}
#Use the following value (+ or -) for a features strand as the default if an invalid value is passed.
#DEFAULT_STRAND  = $config{MISC}{DEFAULT_STRAND}

";

    print $config_text;
    exit(0);

}

#-----------------------------------------------------------------------------

=head1 NAME

gtf2gff3

=head1 VERSION

This document describes version 0.1

=head1 SYNOPSIS

gtf2gff3 --cfg gtf2gff3_MY_CONFIG.cfg gtf_file > gff3_file

=head1 DESCRIPTION

This script will convert GTF formatted files to valid GFF3 formatted
files.  It will map the value in column 3 (\"type\" column) to valid
SO, but because many non standard term may appear in that column in GTF
files, you may edit the config file to provide your own GTF feature to
SO mapping.  The script will also build gene models from exons, CDSs
and other features given in the GTF file.  It is currently tested on
Ensemble and Twinscan GTF, and it should work on any other files that
follow the same specification.  It does not work on GTF from the UCSC
table browser because those files use the same ID for gene and
transcript, so it is impossible to group multiple transcripts to a
gene.  See the README that came with the script for more info.

=head1 OPTIONS:

=over

=item --cfg

Provide the filename for a config file.  See the configuration file
provided with this script for format details.  Use this configuration
file to modify the behavior of the script. If no config file is given
it looks for ./gtf2gff3.cfg, ~/gtf2gff3.cfg or /etc/gtf2gff3.cfg in
that order.

=item --print_cfg

Print a copy of the config variables as currently set (including those
set by --cfg) to STDOUT.  This output can be used as the starting
point for creating a custom config file.

=item --help

Provide a detailed man page style help message and then exit.

=back

=head1 DIAGNOSTICS

=over

=item C<< ERROR : Missing or non-standard attributes : parse_attributes >>

A line in the GTF file did not have any attributes, or it's attributes column
was unparsable.

=item C<< ERROR : Non-transcript gene feature not supported.  Please contact the author for support : build_gene >>

This warning indicates that a line was skipped because it contained a
non-transcript gene feature, and the code is not currently equipped to
handle this type of feature.  This probably isn't too hard to add, so
contact me if you get this error and would like to have these features
supported.

=item C<< ERROR : Must have at least exons or CDSs to build a transcript : build_trnsc >>

Some feature had a transcript_id and yet there were no exons or CDSs
associated with that transcript_id so the script failed to build a
transcript.

=item C<< ERROR : seq_id conflict : validate_and_finish_trnsc >>

Found two features within the same transcript that didn't share the
same seq_id.

=item C<< ERROR : source conflict : validate_and_finish_trnsc >>

Found two features within the same transcript that didn't share the
same source.

=item C<< ERROR : type conflict : validate_and_finish_trnsc >>

Found two features within the same transcript that were expected to
share the same type and yet they didn't.

=item C<< ERROR : strand conflict : validate_and_finish_trnsc >>

Found two features within the same transcript that didn't share the
same strand.

=item C<< ERROR : seq_id conflict : validate_and_build_gene >>

Found two features within the same gene that didn't share the same
seq_id.

=item C<< ERROR : source conflict : validate_and_build_gene >>

Found two features within the same gene that didn't share the same
source.

=item C<< ERROR : strand conflict : validate_and_build_gene >>

Found two features within the same gene that didn't share the same
strand.

=item C<< ERROR : gene_id conflict : validate_and_build_gene >>

Found two features within the same gene that didn't share the same
gene_id.

=item C<< FATAL : Can't open GTF file : file_name for reading. >>

Unable to open the GTF file for reading.

=item C<< FATAL : Need exons or CDSs to build transcripts : process_start >>

A start_codon feature was annotated and yet there were no exons or
CDSs associated with that transcript_id so the script failed.

=item C<< FATAL : Untested code in process_start.  Contact the aurthor for support. >>

The script is written to infer a start codon based on the presence of
a 5' UTR, but we had no example GTF of this type when we wrote the
code, so we killed process rather than run untested code.  Contact the
author for support.

=item C<< FATAL : Invalid feature set : process_start >>

We tried to consider all possible ways of infering a start codon or
infering a a non-coding gene, and yet we've failed.  Your combination
of gene features doesn't make sense to us.  You should never get
this error, and if you do, we'd really like to see the GTF file that
generated it.  Please contact the author for support.

=item C<< FATAL : Need exons or CDSs to build transcripts : process_stop >>

A stop_codon feature was annotated and yet there were no exons or
CDSs associated with that transcript_id so the script failed.

=item C<< FATAL : Untested code in process_stop.  Contact the aurthor for support. >>

The script is written to infer a stop codon based on the presence of
a 3' UTR, but we had no example GTF of this type when we wrote the
code, so we killed process rather than run untested code.  Contact the
author for support.

=item C<< FATAL : Invalid feature set : process_stop >>

We tried to consider all possible ways of infering a stop codon or
infering a a non-coding gene, and yet we've failed.  Your combination
of gene features doesn't make sense to us.  You should never get
this error, and if you do, we'd really like to see the GTF file that
generated it.  Please contact the author for support.

=item C<< FATAL : Invalid feature set : process_exon_CDS_UTR >>

We tried to consider all possible ways of infering exons, CDSs and
UTRs and yet we've failed.  Your combination of gene features doesn't
make sense to us.  You really should ever get this error, and if you
do, we'd really like to see the GTF file that generated it.  Please
contact the author for support.

=item C<< FATAL : Array reference required : sort_features. >>

A user shouldn't be able to trigger this error.  It almost certainly
indicates a software bug.  Please contact the author.

=item C<< FATAL : Can't determine strand in : sort_feature_types. >>

This may indicate that your GTF file does not indicate the strand for
features that require it.  It may also indicate a software bug.
Please contact the author.

=item C<< FATAL : Hash reference required : sort_feature_types. >>

A user shouldn't be able to trigger this error.  It almost certainly
indicates a software bug.  Please contact the author.

=item C<< FATAL : Invalid value passed to strand : strand.  >>

This may indicate that your GTF file does not indicate the strand for
features that require it.  Consider using the DEFAULT_STRAND paramater
in the config file.  It may also indicate a software bug.  Please
contact the author.

=back

=head1 CONFIGURATION AND ENVIRONMENT

A configuration file is provided with this script.  The script will
look for that configuration file in ./gtf2gff3.cfg, ~/gtf2gff3.cfg
or /etc/gtf2gff3.cfg in that order.  If the configuration file is not
found in one of those locations and one is not provided via the --cfg
flag it will try to choose some sane defaults, but you really should
provide the configuration file.  See the supplied configuration file
itself as well as the README that came with this package for format
and details about the configuration file.

=head1 DEPENDENCIES

This script requires the following perl packages that are available
from CPAN (www.cpan.org).

Getopt::Long;
Config::Std;

d1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
<barry.moore@genetics.utah.edu>

=head1 AUTHOR

Barry Moore
<barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2007, University of Utah

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

