#!/usr/bin/env perl

###################################################################
## Installation Instructions:
#
#     Set env variable AATPATH to directory containing the AAT tools
#     including the BS blosum matrix file.
#      
#     * if AATPATH env var is not set, utilities must be available 
#       via the global PATH setting.    
# 
#     * if the BS file is not within AATPATH, the code line below:
#       my $BS = "";  must be set to the complete path to the BS file.
#
###################################################################

# script written by Brian Haas (bhaas@tigr.org) 9/25/2003

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;

my @args = @ARGV;
my $commandLine = join (" ", @args);

## globals:
my ($dds, $dps, $ext, $filter, $nap, $gap2, $show); # AAT binaries


# set paths to utilities and required BS file.
my $AATPATH = $ENV{AATPATH} || $FindBin::Bin;
my $MATRIX_DIR = $ENV{MATRIX_DIR};

my $BS = ($MATRIX_DIR) ? "$MATRIX_DIR/BS" : "BS"; #blosum matrix file.  Set 
$dds = "$AATPATH/dds";
$dps = "$AATPATH/dps";
$ext = "$AATPATH/ext";
$filter = "$AATPATH/filter";
$nap = "$AATPATH/nap";
$gap2 = "$AATPATH/gap2";
$show = "$AATPATH/show";
my $collapse = "$AATPATH/extCollapse.pl";

if ($AATPATH && $BS eq "") {
    $BS = "$AATPATH/BS";
}
else {
    $BS = "$FindBin::Bin/../matrices/BS";
}


my $verbose; #global, print lots to stdout to follow run.
my $btab_only = 0;
my $DONT_DELETE = 0;

main: {
    ## Process commandline arguments:
    #options to AAT utilities
    my ($dds_opts, $dps_opts, $ext_opts, $filter_opts, $nap_opts, $gap2_opts) = ("","","","","","");
    # pipeline specification:
    my ($transcriptPipeline, # dds/gap2 
	$proteinPipeline); #dps/nap
    

    # input files
    my ($queryDB, $searchDB, $unmaskedQueryDB);

    # other options:
    my ($help, $AAThelp);
    my @args = @ARGV;

    &GetOptions('q=s'=>\$queryDB,
		's=s'=>\$searchDB,
		'dds=s'=>\$dds_opts,
		'dps=s'=>\$dps_opts,
		'ext=s'=>\$ext_opts,
		'filter=s'=>\$filter_opts,
		'nap=s'=>\$nap_opts,
		'gap2=s'=>\$gap2_opts,
		'N'=>\$transcriptPipeline,
		'P'=>\$proteinPipeline,
		'h'=>\$help,
		'v'=>\$verbose,
		'AAThelp'=>\$AAThelp,
		'b'=>\$btab_only,
		'X'=>\$DONT_DELETE,
		'unmasked=s'=>\$unmaskedQueryDB);
    

    
    if ( (!@args) || $help) { &usage_and_exit(); };


    unless ($transcriptPipeline xor $proteinPipeline) {
	print STDERR "Please specify \n"
	    . "-P for protein alignment pipeline, or\n"
		. "-N for nucleotide (transcript) alignment pipeline\n\n"
		    . "use -h or no options for usage info.\n\n";
	exit(1);
    }
    
    if ($AAThelp) {
	&advanced_help_and_exit($transcriptPipeline);
    }
    
    unless ($queryDB && $searchDB) {
	print STDERR "Please specify query and search databases using options -q and -s\n\n";
	exit(1);
    }
    

    ## Get filename tokens for output filenames:
    my $searchDBtoken = $searchDB;
    $searchDBtoken =~ s/^.*\/(\S+)$/$1/;
    my $queryDBtoken = $queryDB;
    $queryDBtoken =~ s/^.*\/(\S+)$/$1/;

    ## write errors to errlog file
    my $errlog = "$searchDBtoken.$queryDBtoken.$$.AAT.errlog";
    open (ERRLOG, ">$errlog") or die $!;
    
    ## Run pipeline:
    if ($transcriptPipeline) {
	&run_transcript_pipeline($searchDB, $queryDB, $unmaskedQueryDB, $searchDBtoken, $dds_opts, $ext_opts, $filter_opts, $gap2_opts);
    } else {
	&run_protein_pipeline($searchDB, $queryDB, $unmaskedQueryDB, $searchDBtoken, $dps_opts, $ext_opts, $filter_opts, $nap_opts);
    }
    
    close ERRLOG;
    
    if (-s $errlog) {
	## Errors were encountered, inform user and exit with non-zero status:
	print STDERR "Errors were encountered while running the AAT pipeline.\n"
	    . "Check the file ($errlog) for details.\n\n";
	exit(2);
    } else {
	## no errors were reported, remove the errlog file.
	unlink $errlog;
	exit(0);
    }
}


####
sub run_transcript_pipeline {
    my ($searchDB, $queryDB, $unmaskedQueryDB, $searchDBtoken, $dds_opts, $ext_opts, $filter_opts, $gap2_opts) = @_;
    

    my $fastaReader = Fasta_reader->new($queryDB);

    my $unmaskedReader;
    if ($unmaskedQueryDB) {
	$unmaskedReader = Fasta_reader->new($unmaskedQueryDB);
    }


    while (my $sequence = $fastaReader->next()) {
	
	my $acc = $sequence->{accession};

	my ($unmaskedSeqFile);
	if ($unmaskedReader) {
	    my $unmaskedSeqObj = $unmaskedReader->next();
	    my $unmasked_acc = $unmaskedSeqObj->{accession};
	    if ($unmasked_acc ne $acc) {
		die "Error, unmasked acc: ($unmasked_acc) != masked acc ($acc) ";
	    }
	    $unmaskedSeqFile = $unmaskedSeqObj->write_fasta_file($unmasked_acc . "_unmasked");
	}

	print STDOUT "\n## Processing (AAT-transcript) " . $sequence->{header} . "\n"; 

	## write temp sequence file:
	my $querySeq = $sequence->write_fasta_file();
    
	eval {
	    ## run dds:
	    my $ddsout = "$acc.$searchDBtoken.dds";
	    my $cmd = "$dds $querySeq $searchDB $dds_opts 2>&1 > $ddsout";
	    &process_cmd($cmd);
	    
	    ## run ext
	    my $extout = "$acc.$searchDBtoken.$$.ext";
	    $cmd = "$ext $ddsout $ext_opts 2>&1 > $extout";
	    &process_cmd($cmd);
	    
	    ## Collapse ext overlaps.
	    my $collapsedEXTout = "$acc.$searchDBtoken.$$.extCol";
	    $cmd = "$collapse $extout > $collapsedEXTout";
	    &process_cmd($cmd);
	    
	    ## run filter
	    my $filterout = "$acc.$searchDBtoken.$$.filter";
	    $cmd = "$filter $collapsedEXTout $filter_opts 2>&1 > $filterout";
	    &process_cmd($cmd);
	    
	    ## run gap2
	    my $gap2out = "$acc.$searchDBtoken.gap2";
	    if ($unmaskedSeqFile) {
		$querySeq = $unmaskedSeqFile;
	    }
	    
	    $cmd = "$gap2 $querySeq $searchDB $filterout $gap2_opts 2>&1 > $gap2out";
	    &process_cmd($cmd);
	    
	    if (-e "$querySeq.gap2.btab") {
		rename("$querySeq.gap2.btab", "$acc.$searchDBtoken.gap2.btab");
	    } else {
		## create empty file for btab results
		`touch $acc.$searchDBtoken.gap2.btab`; 
	    } 
	    
	    # remove temp files.
	    unless ($DONT_DELETE) {
		unlink ($ddsout, $extout, $filterout, $querySeq, $collapsedEXTout);
		if ($unmaskedSeqFile) {
		    unlink ($unmaskedSeqFile);
		}
		if ($btab_only) {
		    unlink ($gap2out);
		}
	    }
	    
	};
	
	
	if ($@) {
	    my $err = "----\nError encountered while running pipeline:\n"
		. "CMD: $commandLine\n"
		. "transcript (dds/gap2) pipeline\n"
		. "genomic acc: $acc\n"
		. "searchDB: $searchDB\n"
		. "queryDB: $queryDB\n"
		. "dds_opts: $dds_opts\n"
		. "ext_opts: $ext_opts\n"
		. "filter_opts: $filter_opts\n"
		. "gap2_opts: $gap2_opts\n"
		. $@;
	    print STDERR $err;
	    print ERRLOG $err;
	    
	}
    }
    
    $fastaReader->finish();

}


####
sub run_protein_pipeline {
    my ($searchDB, $queryDB, $unmaskedQueryDB, $searchDBtoken, $dps_opts, $ext_opts, $filter_opts, $nap_opts) = @_;
     
    print "params: $searchDB, $queryDB, $unmaskedQueryDB, $searchDBtoken, $dps_opts, $ext_opts, $filter_opts, $nap_opts\n";
    

    my $fastaReader = Fasta_reader->new($queryDB);


    my $unmaskedReader;
    if ($unmaskedQueryDB) {
	$unmaskedReader = Fasta_reader->new($unmaskedQueryDB);
    }
    


    while (my $sequence = $fastaReader->next()) {
	
	my $acc = $sequence->{accession};

	my ($unmaskedSeqFile);
	if ($unmaskedReader) {
	    my $unmaskedSeqObj = $unmaskedReader->next();
	    my $unmasked_acc = $unmaskedSeqObj->{accession};
	    if ($unmasked_acc ne $acc) {
		die "Error, unmasked acc: ($unmasked_acc) != masked acc ($acc) ";
	    }
	    $unmaskedSeqFile = $unmaskedSeqObj->write_fasta_file($unmasked_acc . "_unmasked");
	    print "unmaskedFile: $unmaskedSeqFile\n";
	}
	
	print STDERR "\n## Processing (AAT-protein) " . $sequence->{header} . "\n"; 
	
	## write temp sequence file:
	my $querySeq = $sequence->write_fasta_file();
	print "querySeq: $querySeq\n";
	
	eval {
	    ## run dps
	    my $dpsout = "$acc.$searchDBtoken.dps";
	    my $cmd = "$dps $querySeq $searchDB $BS $dps_opts 2>&1 > $dpsout";
	    &process_cmd($cmd);
	    
	    ## run ext
	    my $extout = "$acc.$searchDBtoken.$$.ext";
	    $cmd = "$ext $dpsout $ext_opts 2>&1 > $extout";
	    &process_cmd($cmd);

	    ## Collapse ext overlaps.
	    my $collapsedEXTout = "$acc.$searchDBtoken.$$.extCol";
	    $cmd = "$collapse $extout > $collapsedEXTout";
	    &process_cmd($cmd);
	    
	    ## run filter:
	    my $filterout = "$acc.$searchDBtoken.$$.filter";
	    $cmd = "$filter $collapsedEXTout $filter_opts 2>&1 > $filterout";
	    &process_cmd($cmd);
	    
	    ## run nap:
	    my $napout = "$acc.$searchDBtoken.nap";
	    if ($unmaskedSeqFile) {
		$querySeq = $unmaskedSeqFile;
	    }
	    $cmd = "$nap $querySeq $searchDB $filterout $BS $nap_opts 2>&1 > $napout";
	    &process_cmd($cmd);
	    
	    if (-e "$querySeq.nap.btab") {
		rename ("$querySeq.nap.btab", "$acc.$searchDBtoken.nap.btab");
	    } else {
		## create empty btab file for empty results
		`touch $acc.$searchDBtoken.nap.btab`;
	    } 
	    
	    # remove temp files.
	    unless ($DONT_DELETE) {
		unlink ($dpsout, $extout, $filterout, $querySeq, $collapsedEXTout);
		if ($unmaskedSeqFile) {
		    unlink ($unmaskedSeqFile);
		}
		if ($btab_only) {
		    unlink ($napout);
		}
	    }
	    


	};
	
	if ($@) {
	    my $err = "----\nError encountered while running pipeline:\n"
		. "CMD: $commandLine\n"
		. "protein (dps/nap) pipeline\n"
		. "genomic acc: $acc\n"
		. "searchDB: $searchDB\n"
		. "queryDB: $queryDB\n"
		. "dps_opts: $dps_opts\n"
		. "ext_opts: $ext_opts\n"
		. "filter_opts: $filter_opts\n"
		. "nap_opts: $nap_opts\n"
		. $@;

	    print STDERR $err;
	    print ERRLOG $err;
	}
    }
    
    $fastaReader->finish();
}



####
sub process_cmd {
    my $cmd = shift;

    print STDOUT "CMD: $cmd\n";
    my $response = `$cmd`;
    my $ret = $?;
    print STDOUT "ret($?)\n";
    if ($?) {
	print STDERR "Error, command failed: $cmd, returns($ret), response: $response\n";
	die "ERROR running cmd: $cmd, returns($ret), reponse: $response\n";
    }
        
    return (0);
}




####
sub usage_and_exit {
    print <<_EOHELP;

########### AAT Usage ##############################################################
#
# Required:

#  Flags:
#
#  -N  :nucleotide (transcript) spliced alignment pipeline (dds/ext/filter/nap)
#
#  -P  :protein spliced alignment pipeline (dps/ext/filter/nap)
# 
#  -b  :btab files only (deletes the nap alignment outputs)
#  -X  :don't delete the intermediate files of the pipeline
#
#  Parameters:
#
#  -q 'param' :query database (fasta or multi-fasta of genomic sequence(s) )
#
#  -s 'param' :search database (nucleotide or protein fasta or multi-fasta sequence database)
#
#  --unmasked 'param'  :if -q provides a masked sequence, here you can provide the unmasked sequence.  The 
#                       masked sequence will be used in the hit-generation stage, and the unmasked sequence
#                       will be used for generating the global alignments.
#
# Optional:
# 
#  --dds 'params' :parameters passed to dds
#
#  --gap2 'params'  :parameters passed to gap2
#
#  --dps 'params'  :parameters passed to dps
#
#  --nap 'params' :parameters passed to nap
#
#  --ext 'params' :parameters passed to ext
#
#  --filter 'params' :parameters passed to filter
#
#  --AAThelp :provides help menu for all programs in pipeline (must specify -P or -N)
# 
#  Examples:
#  
#   protein pipeline:
#        AAT.pl -P -q genomic_db -s protein_db --dps '-f 100 -i 30 -a 200' --filter '-c 10' --nap '-x 10'
#
#   nucleotide (transcript) pipeline:
#        AAT.pl -N -q genomic_db -s cDNA_transcript_db --dds '-f 100 -i 20 -o 75 -p 70 -a 2000' --filter '-c 10' --gap2 '-x 1'
#
##########################################################################################



_EOHELP

    ;

    exit(1);
}


#### 
sub advanced_help_and_exit {
    my $transcriptPipeline = shift;
    my @cmds;
    if ($transcriptPipeline) {
	@cmds = qw (dds ext filter gap2);
    } else {
	@cmds = qw (dps ext filter nap);
    }
    
    foreach my $cmd (@cmds) {
	# run command with no params so just help options print.
	print "\n\n************ $cmd *****************\n";
	system $cmd;
    }
    exit(1);

}




##########################################################

# lightweight fasta reader capabilities:
package Fasta_reader;

use strict;

sub new {
    my ($packagename, $fastaFile) = @_;
    my $self = { fastaFile => $fastaFile,
		 fileHandle => undef };
    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
    $self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;
    my $seqobj = undef;
    if ($next_text_input) {
	$next_text_input =~ s/^>|>$//g; #remove trailing > char.
	my ($header, @seqlines) = split (/\n/, $next_text_input);
	my $sequence = join ("", @seqlines);
	$seqobj = Sequence->new($header, $sequence);
    }
    
    $/ = $orig_record_sep; #reset the record separator to original setting.
    
    return ($seqobj); #returns null if not instantiated.
    
}

#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}




##############################################
package Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;
    
    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
    $acc =~ s/\W/_/g; #convert non-word chars in acc.
    
    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
		 filename => undef };
    bless ($self, $packagename);
    return ($self);
}

####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;
    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});
    $sequence =~ s/(\w{60})/$1\n/g;
    my $time = time();
    
    my $tempfile = $filename;
    unless ($tempfile) {
	$tempfile = $accession;
    }
    $tempfile =~ s/\W/_/g;
            
    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP ">$header\n$sequence";
    close TMP;
    return ($tempfile);
}

1; #EOM


