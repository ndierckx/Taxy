#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2023-2024}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#           Taxy -Taxonomy assignment of ASVs
#           nicolasdierckxsens@hotmail.com
use Getopt::Long;
use strict;
use MCE::Child;
use MCE::Channel;

my $reads12 = "";
my $reads1 = "";
my $reads2 = "";
my $nt_path = "";
my $tax_db = "";
my $output_path = "";

my %blast_files;
my %blast_files_id;
my $go_to_blast_files = "";

my $paired = "SE";
my $batch_file = "";
my $use_hashes_from_file = "";
my $use_quality = "";
my $taxonomy_only = "";
my $direct_ASV_output = "";
my %hash;
my %filehandle;
my %reads;
my $config = "";
my $project = "";
my $blast_batch_size = '300';
my $maxProcs = '1';

#Read the config file----------------------------------------------------------------------------------------------

GetOptions (
            "c=s" => \$config,
            ) or die "Incorrect usage!\n";

open(CONFIG, $config) or die "Error:Can't open the configuration file, please check the manual!\n\nUsage: perl eDNA.pl -c config.txt\n";

while (my $line = <CONFIG>)
{
    $line =~ tr/\r//d;
    $line =~ s/\R/\012/;
    if ($line =~ m/.*Project name\s+\=\s+(.*?)(Combined reads.*)*$/)
    {
        $project = $1;
        chomp $project;
        my $project_tmp = $project;
        my $ggg;
        if ($project =~ m/batch\:(.*)/)
        {
            my $batch_file_tmp = $1;
            if ($batch_file eq "")
            {
                $batch_file = $batch_file_tmp;
                print "Batch file detected...\n\n";
                open(BATCH, $batch_file) or die "Error: $!\nCan't open the batch file, please check the manual!\n";
                $ggg = "yes";
            }
            while (my $line = <BATCH>)
            {
                $project = $line;
                chomp $project;
                last;
            }
            if ($project eq $project_tmp || $project eq "")
            {
                goto EXIT;
            }
            elsif ($ggg ne "yes")
            {
                print "\n\n------------------------------\n------------------------------\n";
                print "        NEXT SAMPLE\n";
                print "------------------------------\n------------------------------\n\n\n";
            }
        }
    }
    if ($line =~ m/.*Combined reads or ASVs\s+\=\s+(.*?)(Forward reads.*)*$/)
    {
        $reads12 = $1;
        chomp $reads12;
        if ($reads12 eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $reads12 = $line;
                chomp $reads12;
                last;
            }
        }    
    }
    if ($line =~ m/.*Forward reads\s+\=\s+(.*?)(Reverse reads.*)*$/)
    {
        $reads1 = $1;
        chomp $reads1;
        if ($reads1 eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $reads1 = $line;
                chomp $reads1;
                last;
            }
        }    
    }     
    if ($line =~ m/.*Reverse reads\s+\=\s+(.*?)(Keep read ids.*)*$/)
    {
        $reads2 = $1;
        chomp $reads2;
        if ($reads2 eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $reads2 = $line;
                chomp $reads2;
                last;
            }
        }
    }
    if ($line =~ m/.*Keep read ids\s+\=\s+(.*?)(Nucleotide database.*)*$/)
    {
        $direct_ASV_output = $1;
        chomp $direct_ASV_output;
        if ($direct_ASV_output eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $direct_ASV_output = $line;
                chomp $direct_ASV_output;
                last;
            }
        }
    }
    if ($line =~ m/.*Nucleotide database\s+\=\s+(.*?)(Taxonomy database.*)*$/)
    {
        $nt_path = $1;
        chomp $nt_path;
        if ($nt_path eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $nt_path = $line;
                chomp $nt_path;
                last;
            }
        }
    }
    if ($line =~ m/.*Taxonomy database\s+\=\s+(.*?)(Taxonomy only.*)*$/)
    {
        $tax_db = $1;
        chomp $tax_db;
        if ($tax_db eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $tax_db = $line;
                chomp $tax_db;
                last;
            }
        }
    }
    if ($line =~ m/.*Taxonomy only\s+\=\s+(.*?)(Threads.*)*$/)
    {
        $taxonomy_only = $1;
        chomp $taxonomy_only;
        if ($taxonomy_only eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $taxonomy_only = $line;
                chomp $taxonomy_only;
                last;
            }
        }
    }
    if ($line =~ m/.*Threads\s+\=\s+(.*?)(Output path.*)*$/)
    {
        $maxProcs = $1;
        chomp $maxProcs;
        if ($maxProcs eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $maxProcs = $line;
                chomp $maxProcs;
                last;
            }
        }
    }
    if ($line =~ m/.*Output path\s+\=\s+(.*?)$/)
    {
        $output_path = $1;
        chomp $output_path;
        last;
    }
}

close CONFIG;

my $output_path_test = substr $output_path, -1 , 1;

if ($output_path_test ne "\\" && $output_path_test ne "/"  &&  $output_path ne "")
{
    die "\nCan not recognize the output path, it should end with a directory separator: $output_path, $!\n";                         
}
elsif ($output_path ne "")
{
    unless (-d $output_path)
    {
        mkdir $output_path;
    }
}


my $directory = $output_path."tmp_sequences";
mkdir $directory;

if ($paired ne "PE" && $paired ne "SE")
{
    die "\nPaired has to be 'SE' or 'PE', please check the configuration file!\n";
}
if ($reads12 ne "" && $reads1 ne "")
{
    die "\nYou can't give a path for a combined dataset and a forward and reverse set!\n If you have both, only use the forward and reverse path in the config file\n";
}

my $chnl;  
$chnl = MCE::Channel->new( impl => 'Simple' );

#spin up worker early before creating big hash---------
mce_child
{
    local $SIG{__WARN__} = sub {};
    while ( my ($cmd, @args) = $chnl->recv ) {
        local ($?, $!);
        system($cmd, @args);
        $chnl->send2($?, $!);
    }
};

sub syscmd {
    my $cmd = shift;
    return unless $cmd;

    $chnl->send($cmd, @_);
    my ($status, $errmsg) = $chnl->recv2;
    
    if ($status == -1) {
        print "SYSTEM: failed to execute ($cmd): $errmsg\n";
    }
    elsif ($status & 127) {
        printf "SYSTEM: $cmd died with signal %s, %s coredump\n",
            ($status & 127), ($status & 128) ? 'with' : 'without';
    }
    else {
        #printf "SYSTEM: $cmd exited with status %d\n", $status >> 8;
    }
}  

my %tax_db;
undef %tax_db;

if ($tax_db eq "dgdg")
{
    open(TAXDB, $tax_db) or die "No input file found, make sure the path is correct $!\n";
    my $parent_id = "";
    my $tax_id = "";
    my $tax_name = "";
    my $rank = "";
    
    while (my $line = <TAXDB>)
    {
        chomp($line);
        my $id_check = substr $line, 0, 2;
        my $parent_check = substr $line, 0, 9;
        my $name_check = substr $line, 0, 15;
        my $rank_check = substr $line, 0, 4;
        
        if ($id_check eq "ID")
        {
            $tax_id = substr $line, 28;     
        }
        elsif ($parent_check eq "PARENT ID")
        {
            $parent_id = substr $line, 28;        
        }
        elsif ($rank_check eq "RANK")
        {
            $rank = substr $line, 28;        
        }
        elsif ($name_check eq "SCIENTIFIC NAME")
        {
            $tax_name = substr $line, 28;
            $tax_db{$tax_id} = $parent_id.";".$tax_name.";".$rank;
        }
    }
    close TAXDB;
}

if ($tax_db ne "")
{
    open(TAXDB, $tax_db) or die "No input file found, make sure the path is correct $!\n";
    
    while (my $line = <TAXDB>)
    {
        chomp($line);
        my $taxon_check = substr $line, 0, 7;
        my $parent_id = "";
        my $tax_id = "";
        my $tax_name = "";
        my $rank = "";
        
        if ($taxon_check eq "<taxon ")
        {
            my @line = split /" /, $line;
            
            foreach my $part_tmp (@line)
            {
                my $id_check = substr $part_tmp, 0, 5;
                my $id_check2 = substr $part_tmp, 0, 12;
                if ($id_check eq "taxId")
                {
                    $tax_id = substr $part_tmp, 7;
                }
                elsif ($id_check eq "paren")
                {
                    $parent_id = substr $part_tmp, 13;                
                }
                elsif ($id_check eq "rank=")
                {
                    $rank = substr $part_tmp, 6;      
                }
                elsif ($id_check2 eq "<taxon scien")
                {
                    $tax_name = substr $part_tmp, 23;
                }
            }
            if ($tax_id eq "1")
            {
                $parent_id = "0";
            }
            $tax_db{$tax_id} = $parent_id.";".$tax_name.";".$rank;
        }
    }
    close TAXDB;
}

my %taxonomy_ranks;
undef %taxonomy_ranks;

$taxonomy_ranks{"subspecies"} = '1';
$taxonomy_ranks{"species"} = '2';
$taxonomy_ranks{"subgenus"} = '3';
$taxonomy_ranks{"genus"} = '4';
$taxonomy_ranks{"family"} = '5';
$taxonomy_ranks{"order"} = '6';
$taxonomy_ranks{"class"} = '7';
$taxonomy_ranks{"phylum"} = '8';
$taxonomy_ranks{"kingdom"} = '9';
$taxonomy_ranks{"superkingdom"} = '10';

my $output_file4  = $output_path."log.txt";
open(OUTPUT_LOG, ">" .$output_file4) or die "\nCan't open file $output_file4, $!\n";

if ($go_to_blast_files eq "yes")
{
    goto BLAST_FILES;
}
my @reads_tmp = undef;

if ($reads12 eq "")
{
    @reads_tmp = ($reads1, $reads2);
    if ($reads1 eq $reads2)
    {     
        if ($batch_file ne "")
        {
            print "\nThe two input files are identical, please check the configuration file!\n";
            #print OUTPUT4 "\nThe two input files are identical, please check the configuration file!\n";
            goto BATCH;
        }
        else
        {
            die "\nThe two input files are identical, please check the configuration file!\n";
        }
    }
}
else
{
    @reads_tmp = ($reads12);
}

my $check_zip = substr $reads_tmp[0], -2;
my $check_zip2 = substr $reads_tmp[0], -3;
my $firstLine = "";
my $secondLine = "";
my $thirdLine = "";
my $fourthLine = "";
my $fifthLine = "";

if ($check_zip eq "gz" && $use_hashes_from_file eq "")
{
    open (my $FILE, '-|', 'gzip', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
    $firstLine = <$FILE>;
    chomp $firstLine;
    $secondLine = <$FILE>;
    chomp $secondLine;
    $thirdLine = <$FILE>;
    chomp $thirdLine;
    $fourthLine = <$FILE>;
    chomp $fourthLine;
    $fifthLine = <$FILE>;
    chomp $fifthLine;
    close $FILE;
}
elsif ($check_zip2 eq "bz2" && $use_hashes_from_file eq "")
{
    open (my $FILE, '-|', 'bzip2', '-dc', $reads_tmp[0]) or die "Can't open file $reads_tmp[0], $!\n";
    $firstLine = <$FILE>;
    chomp $firstLine;
    $secondLine = <$FILE>;
    chomp $secondLine;
    $thirdLine = <$FILE>;
    chomp $thirdLine;
    $fourthLine = <$FILE>;
    chomp $fourthLine;
    $fifthLine = <$FILE>;
    chomp $fifthLine;
    close $FILE;
}
elsif ($use_hashes_from_file eq "")
{
    open(INPUT, $reads_tmp[0]) or die "No input file found, make sure it are fastq files $!\n";
    $firstLine = <INPUT>;
    chomp $firstLine;
    $secondLine = <INPUT>;
    chomp $secondLine;
    $thirdLine = <INPUT>;
    chomp $thirdLine;
    $fourthLine = <INPUT>;
    chomp $fourthLine;
    $fifthLine = <INPUT>;
    chomp $fifthLine;
    close INPUT;
}


select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nReading Input...";

my $firstLine_reverse = "";
if ($reads12 eq "" && $use_hashes_from_file eq "")
{
    if ($check_zip eq "gz")
    {
        open (my $FILE, '-|', 'gzip', '-dc', $reads_tmp[1]) or die "Can't open file $reads_tmp[1], $!\n";
        $firstLine_reverse = <$FILE>;
        chomp $firstLine_reverse;
        close $FILE;
    }
    elsif ($check_zip2 eq "bz2")
    {
        open (my $FILE, '-|', 'bzip2', '-dc', $reads_tmp[1]) or die "Can't open file $reads_tmp[1], $!\n";
        $firstLine_reverse = <$FILE>;
        chomp $firstLine_reverse;
        close $FILE;
    }
    else
    {
        open(INPUT2, $reads_tmp[1]) or die "\n\nNo input file found, make sure it are fastq files $!\n";
        $firstLine_reverse = <INPUT2>;
        chomp $firstLine_reverse;
        close INPUT2;
    }
}
my $no_quality_score = substr $thirdLine, 0, 1;
my $type_of_file = "";
my $type_of_file2 = "";
my $code_before_end = substr $firstLine, -2,1;
my $code_before_end0 = substr $firstLine, -1,1;
my $code_before_end_reverse = substr $firstLine_reverse, -2,1;
my $code_before_end0_reverse = substr $firstLine_reverse, -1,1;

my $SRA = "";
if ($paired eq "SE")
{
    $type_of_file = '0';
}
elsif (($code_before_end eq "/" || $code_before_end eq "R" || $code_before_end eq "#") && $firstLine_reverse ne $firstLine && $code_before_end0 eq "1")
{
    $type_of_file = '-1';
}
elsif ($code_before_end eq ":" && $code_before_end0 eq "1" && $firstLine_reverse ne $firstLine && $code_before_end_reverse eq ":" && $code_before_end0_reverse eq "2")
{
    $type_of_file = '-1';
}
elsif($firstLine =~ m/.*(_|\s)(1)(:\w.*\d+:*(\s.*|\w+)*\s*\t*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = "yes";
}
elsif ($code_before_end eq ":" && $code_before_end0 eq "1" && $firstLine_reverse ne $firstLine && $firstLine_reverse eq "")
{
    $type_of_file = '-1';
}
elsif($firstLine =~ m/.*\s(1)(\S*)$/ && $firstLine_reverse ne $firstLine)
{
    my $firstLine_tmp = $firstLine;
    my $test_space = $firstLine_tmp =~ tr/ //;
    $type_of_file = -length($2)-1;
    if ($test_space eq '1')
    {
        $type_of_file = "split";
    }
}
elsif($firstLine =~ m/.*_(1)(:N.*)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($fifthLine =~ m/.*\.(1)(\s(.+)\s.+)$/ && $firstLine_reverse ne $firstLine)
{
    $type_of_file = -length($2)-1;
}
elsif($firstLine =~ m/.*(length=\d+\s*)$/)
{
    $type_of_file = "identical";
    $type_of_file2 = -length($1);
}
elsif($firstLine_reverse eq $firstLine)
{
    $type_of_file = "identical";
}
elsif($firstLine =~ m/\S*\.(1)(\s(\d+)\s.*)$/ && $firstLine_reverse ne $firstLine)
{
    my $test1 = $3;
    if($fifthLine =~ m/\S*\.(1)(\s(\d+)(\s.*))$/ && $firstLine_reverse ne $firstLine)
    {
         my $test2 = $3;
         if ($test2 eq $test1)
         {
            $type_of_file = -length($2)-1;
         }
         else
         {
            $type_of_file = -length($4);
            $SRA = "yes";
         }
    }
}
elsif($reads12 ne "")
{
    print "\n\nCOMBINED FILE NOT SUPPORTED, PLEASE TRY SEPERATE FILES FOR THE FORWARD AND REVERSE READS!\n\n";
    print OUTPUT4 "\n\nCOMBINED FILE NOT SUPPORTED, PLEASE TRY SEPERATE FILES FOR THE FORWARD AND REVERSE READS!\n\n";
    if ($batch_file ne "")
    {
        goto BATCH;
    }
    else
    {
        exit;
    }
}
else
{
    print "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    print OUTPUT4 "\n\nTHE INPUT READS HAVE AN INCORRECT FILE FORMAT!\nPLEASE SEND ME THE ID STRUCTURE!\n\n";
    if ($batch_file ne "")
    {
        goto BATCH;
    }
    else
    {
        exit;
    }
}

my $check_line_end = $secondLine;
chomp($check_line_end);
$check_line_end =~ tr/\r//d;
$check_line_end =~ s/\R/\012/;
my $last_character = substr $check_line_end, -1;
my $space_at_end = "";
if ($last_character =~ m/\s|\t/g)
{
    #print OUTPUT5 $last_character." LAST2\n";
    $space_at_end = "yes";
}
print "...OK\n";

my $out_of_memory = "";
my $code_new = '1';
my $file_count= '0';
my $file1_count = '0';
my $file2_count = '0';
my $count_paired = '0';
my $skipped_reads = '0';
my $keys_hash = '0';
my $count_batch = '0';
my $chop = "";
my $no_pair1 = '0';
my $no_pair2 = '0';
my $last_print_key = "";
my $no_merge = '0';
my $no_id_match = '0';
my $identical_reads = '0';

my $FILE;
my $FILE2;
my $total_ASVs = '0';

foreach my $reads_tmp (@reads_tmp)
{  
    $file_count++;
    if ($file_count eq "1")
    {
        if ($check_zip eq "gz")
        {
            open ($FILE, '-|', 'gzip', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
            my $count_ASVs = qx(zgrep -c ">" $reads_tmp);
            chomp($count_ASVs);
            $total_ASVs += $count_ASVs;
        }
        elsif ($check_zip2 eq "bz2")
        {
            open ($FILE, '-|', 'bzip2', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
            my $count_ASVs = qx(zgrep -c ">" $reads_tmp);
            chomp($count_ASVs);
            $total_ASVs += $count_ASVs;
        }
        else
        {
            open($FILE, $reads_tmp) or die "\n\nCan't open file $reads_tmp, $!\n";
            my $count_ASVs = qx(grep -c ">" $reads_tmp);
            chomp($count_ASVs);
            $total_ASVs += $count_ASVs;
        }
    }
    elsif ($file_count eq "2")
    {
        if ($check_zip eq "gz")
        {
            open ($FILE2, '-|', 'gzip', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
        }
        elsif ($check_zip2 eq "bz2")
        {
            open ($FILE2, '-|', 'bzip2', '-dc', $reads_tmp) or die "Can't open file $reads_tmp, $!\n";
        }
        else
        {
            open($FILE2, $reads_tmp) or die "\n\nCan't open file $reads_tmp, $!\n";
        }
    }
}
$blast_batch_size = int($total_ASVs/$maxProcs);
print "\nSequences per batch : ".$blast_batch_size."\n";

my $N = '0';
my $f = "";
my $code1 = "";
my $code2 = "";
my $id_tmp0 = "";
$code_new = '1';
my $value1 = "";
my $value2 = "";
my $quality1 = "";
my $quality2 = "";
$out_of_memory = "";
my $output_file1 = "";
my %identical_reads;
my $first_id = "";
my $read_merged = "";
my %original_ASVs;
undef %original_ASVs;

FILE_LINE:while (my $line1 = <$FILE> or my $line2 = <$FILE2>)
{   
    $file1_count++;
    $file2_count++;

    if ($out_of_memory eq "yes")
    {    
        last FILE_LINE;
    }
    
    chomp $line1;
    chomp $line2;
    $line1 =~ tr/\r//d;
    $line1 =~ s/\R/\012/;
    $line2 =~ tr/\r//d;
    $line2 =~ s/\R/\012/;
   
    my $first_char_check = "";
    if ($taxonomy_only ne "")
    {
        my $first_char = substr $line1, 0, 1;
        $first_char =~ tr/actgn/ACTGN/;
        if ($first_char eq "A" || $first_char eq "C" || $first_char eq "T" || $first_char eq "G" || $first_char eq "N")
        {
            $first_char_check = "yes";
        }
    }

    if ($f eq "use_quality")
    {
        if ($use_quality ne "")
        {
            $value1 = $line1;
            $value2 = $line2;
        }
        $f = "yes";
        next FILE_LINE;
    }
    if ($f eq "yes" && $no_quality_score ne "@" && $no_quality_score ne ">" && $taxonomy_only eq "")
    {
        $f = "yes2";
        next FILE_LINE;
    }
    if ($f eq "yes" && ($no_quality_score eq "@" || $no_quality_score eq ">")  && $taxonomy_only eq "")
    {
        $code1 = substr $line1, 1;
        $code2 = substr $line2, 1;
        $f = "no";
        next FILE_LINE;
    }
    if ($f eq "yes2")
    {
        if ($use_quality ne "")
        {
            $f = "no";
            $quality1 = $line1;
            $quality2 = $line2;
        }
        else
        {
            $f = "";
            next FILE_LINE;
        }
    }
    if ($f eq "no")
    {
        if ($use_quality eq "")
        {
            $value1 = $line1;
            $value2 = $line2;
        }

        if ($space_at_end eq "yes")
        {
            chop($value1);
            chop($value2);
        }

        my $code1a = $code1;
        my $code2a = $code2;
        my $code_end1 = substr $code1a, $type_of_file, 1;
        my $code_end2 = substr $code2a, $type_of_file, 1;
        my $code0_SRA1 = "";
        my $code0_SRA2 = "";
        if ($paired eq 'SE')
        {
            $code_end1 = '2';
        }

        if ($SRA eq "yes")
        {
            my $code_SRA1 = substr $code1a, 0, $type_of_file;
            my @code_SRA1 = split / /, $code_SRA1;
            $code0_SRA1 = $code_SRA1[0];
            $code_end1 = substr $code0_SRA1, -1, 1, "";
            
            my $code_SRA2 = substr $code2a, 0, $type_of_file;
            my @code_SRA2 = split / /, $code_SRA2;
            $code0_SRA2 = $code_SRA1[0];
            $code_end2 = substr $code0_SRA2, -1, 1, "";
        }                  
        if ($type_of_file eq "identical")
        {
            $code_end1 = "1";
            $code_end2 = "2";   
        }
        if ($type_of_file eq "yes")
        {
            if($code1a =~ m/.*(_|\s)(1|2)(:\w.*\d+:*(\s.*|\w+)*\s*\t*)$/)
            {
                $code_end1 = $2;
                $type_of_file2 = -length($3)-1;
            }
            if($code2a =~ m/.*(_|\s)(1|2)(:\w.*\d+:*(\s.*|\w+)*\s*\t*)$/)
            {
                $code_end2 = $2;
            }  
        }
        if ($type_of_file eq "split")
        {
            my @split1 = split / /, $code1a;
            $code_end1 = substr $split1[1], 0, 1;
            $type_of_file2 = $split1[0];
            
            my @split2 = split / /, $code2a;
            $code_end2 = substr $split2[1], 0, 1;
        }
        if ($out_of_memory ne "yes")
        {
            if ($use_quality ne "")
            {
                my $quality1a = $quality1;
                my $check_quality1 = $quality1a =~ tr/\!|\"|\#|\$|\%|\&|\'|\(|\)|\*|\+|\,|\-|\.|\/|0|1|2|3|4|5/O/;
                
                my $quality2a = $quality2;
                my $check_quality2 = $quality2a =~ tr/\!|\"|\#|\$|\%|\&|\'|\(|\)|\*|\+|\,|\-|\.|\/|0|1|2|3|4|5/O/;
                
                if ($check_quality1 > 0)
                {
                    my $offset = '0';
                    my $result = index($quality1a, 'O', $offset);
                    
                    while ($result != -1)
                    {                        
                        my $nuc1 = substr $value1, $offset, 1;
                        if ($nuc1 eq "A")
                        {
                            substr $value1, $offset, 1, "1";
                        }
                        elsif ($nuc1 eq "C")
                        {
                            substr $value1, $offset, 1, "2";
                        }
                        elsif ($nuc1 eq "T")
                        {
                            substr $value1, $offset, 1, "3";
                        }
                        elsif ($nuc1 eq "G")
                        {
                            substr $value1, $offset, 1, "4";
                        }
                        $offset = $result + 1;
                        $result = index($quality1a, 'O', $offset);
                    }
                }
                
                if ($check_quality2 > 0)
                {
                    my $offset = '0';
                    my $result = index($quality2a, 'O', $offset);
                    
                    while ($result != -1)
                    {                        
                        my $nuc2 = substr $value2, $offset, 1;
                        if ($nuc2 eq "A")
                        {
                            substr $value2, $offset, 1, "1";
                        }
                        elsif ($nuc2 eq "C")
                        {
                            substr $value2, $offset, 1, "2";
                        }
                        elsif ($nuc2 eq "T")
                        {
                            substr $value2, $offset, 1, "3";
                        }
                        elsif ($nuc2 eq "G")
                        {
                            substr $value2, $offset, 1, "4";
                        }
                        $offset = $result + 1;
                        $result = index($quality2a, 'O', $offset);
                    }
                }
            }
            my $code1_0 = substr $code1a, 0, $type_of_file;
            my $code2_0 = substr $code2a, 0, $type_of_file;
            if ($SRA eq "yes")
            {
                $code1_0 = $code0_SRA1;
                $code2_0 = $code0_SRA2;
            }
         
            if ($type_of_file eq "identical")
            {
               $code1_0 = $code1a;
               $code2_0 = $code2a;
               if ($type_of_file2 ne "")
               {
                    $code1_0 = substr $code1a, 0, $type_of_file2;
                    $code2_0 = substr $code2a, 0, $type_of_file2;
               }
            }
            if ($type_of_file eq "yes")
            {
               $code1_0 = substr $code1a, 0, $type_of_file2;
               $code2_0 = substr $code2a, 0, $type_of_file2;
            }
            if ($type_of_file eq "split")
            {
                $code1_0 = $type_of_file2;
                $code2_0 = $type_of_file2;
            }
#                if ($memory_max_current > $max_memory && $max_memory ne "")
           # {
               # $out_of_memory = "yes";          
           # }
            $keys_hash++;
            
            if ($chop eq "")
            {
                my $last_character = substr $code1_0, -1;       
                if ($last_character =~ m/\s|\t/g)
                {
                    $chop = "yes";
                }
                else
                {
                    $chop = "no";
                }
            }
            
            if ($code1_0 ne $code2_0)
            {
                $no_id_match++;
            }
            
            my $value2_reverse = reverse($value2);
            $value2_reverse =~ tr/ACTG/TGAC/;
            
            my $value_merged = "";
            my $f = 20;
            
            if ($taxonomy_only eq "")
            {
                while ($f < length($value1))
                {
                    my $seq_part = substr $value1, -$f, 20;
                    
                    my $c = 0;
                    while ($c < length($value2_reverse))
                    {
                        my $seq_part2 = substr $value2_reverse, $c, 20;
                        
                        if ($seq_part eq $seq_part2)
                        {
                            my $seq_merge1 = substr $value1, 0, -$f;
                            my $seq_merge2 = substr $value2_reverse, $c;
                            $value_merged = $seq_merge1.$seq_merge2;
                            goto SEQ_MERGED;
                        }
                        $c++;   
                    }
                    
                    $f++;
                }
                $read_merged = $value_merged;
            }
            else
            {
                $read_merged = $value1;
            }
SEQ_MERGED:          
        }
        $f = "yes";
        if ($use_quality ne "")
        {
            $f = "";
        }
        unless (eof)
        {
            next FILE_LINE;
        }
    }
    elsif ($taxonomy_only ne "" && $first_char_check eq "yes")
    {
        $read_merged .= $line1;
        unless (eof)
        {
            next FILE_LINE;
        }
    }
    
    $id_tmp0 = $keys_hash;

    if ($direct_ASV_output ne "")
    {
        $id_tmp0 = $code1;
    }
    
    if ($first_id eq "")
    {
        $first_id = "yes";
        $code1 = substr $line1, 1;
        $code2 = substr $line2, 1;
        if ($direct_ASV_output ne "")
        {
            $id_tmp0 = $code1;
        }
    }
    elsif ($read_merged ne "")
    {
        if ($taxonomy_only eq "")
        {
            substr $read_merged, -33, 33, "";
            substr $read_merged, 0, 25, "";
            #substr $value_merged, -15, 15, "";
            #substr $value_merged, 0, 10, "";
        }
        if ($direct_ASV_output eq "")
        {
            $original_ASVs{$id_tmp0} = $read_merged;
            if (exists($reads{$read_merged}))
            {
                $identical_reads++;
                if (exists($identical_reads{$reads{$read_merged}}))
                {
                    my $count_tmp = $identical_reads{$reads{$read_merged}};
                    $identical_reads{$reads{$read_merged}} = $count_tmp+1;
                }
                else
                {
                    $identical_reads{$reads{$read_merged}} = '2';
                } 
                goto SKIP2;
            }
            else
            {
                $reads{$read_merged} = $keys_hash;    
            }
        }
        else
        {
            $original_ASVs{$id_tmp0} = $read_merged;
            if (exists($reads{$read_merged}))
            {
                $reads{$read_merged} .= exists $reads{$read_merged} ? ",$code1" : $code1;
                goto SKIP2;
            }
            $reads{$read_merged} = $code1;       
        }
       
        my $FILE_READ;
        $blast_files_id{$code1} = undef;
        if ($code1 eq "")
        {
            print OUTPUT_LOG $read_merged." EMPTY_ID\n";
        }

        if ($read_merged eq "" && $taxonomy_only eq "")
        {
            #print $value1."\n";
            #print $value2_reverse."\n";
            $no_merge++;
        }
        elsif ($count_batch eq '0')
        {
            $output_file1  = $directory."/sequence_tmp_NP_".$keys_hash.".fasta";
            $filehandle{$keys_hash} = $FILE_READ;
            $last_print_key = $keys_hash;
            open($filehandle{$keys_hash}, ">" .$output_file1) or die "\nCan't open file $output_file1, $!\n";
            print {$filehandle{$last_print_key}} ">".$id_tmp0."\n";
            print {$filehandle{$last_print_key}} $read_merged."\n";
            $count_batch++;
        }
        elsif ($count_batch < $blast_batch_size)
        {
            print {$filehandle{$last_print_key}} ">".$id_tmp0."\n";
            print {$filehandle{$last_print_key}} $read_merged."\n";
            $count_batch++;  
        }
        else
        {
            print {$filehandle{$last_print_key}} ">".$id_tmp0."\n";
            print {$filehandle{$last_print_key}} $read_merged;
            close $filehandle{$last_print_key};
            my $command = "blastn -query ".$output_file1." -db ".$nt_path." -out ".$output_path."blast_tmp_".$keys_hash.".txt -num_threads 1 -qcov_hsp_perc 65 -max_target_seqs 30 -outfmt \"6 qseqid sseqid staxids pident qcovs qseq length qlen slen mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovhsp stitle\" &";  
            syscmd($command);
            $blast_files{$keys_hash} = $output_path."blast_tmp_".$keys_hash.".txt";        
            $count_batch = '0';   
        }

SKIP2:     
        $code1 = substr $line1, 1;
        $code2 = substr $line2, 1;
        $f = "no";
        if ($use_quality ne "")
        {
            $f = "use_quality";
        }
    } 
}
close $filehandle{$last_print_key};

my $command = "blastn -query ".$output_file1." -db ".$nt_path." -out ".$output_path."blast_tmp_".$keys_hash.".txt -num_threads 1 -qcov_hsp_perc 65 -max_target_seqs 30 -outfmt \"6 qseqid sseqid staxids pident qcovs qseq length qlen slen mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovhsp stitle\"";
syscmd($command);
$blast_files{$keys_hash} = $output_path."blast_tmp_".$keys_hash.".txt";

sleep(10);
BLAST_FILES:
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------Read BLAST files--------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------

my %ASV_table;
undef %ASV_table;

my $ASV_table = $output_path."ASV_Table_".$project.".txt";
open(ASV_TABLE, ">". $ASV_table) or die "\nCan't open file $ASV_table, $!\n";

BLAST0:
my $count_remaining_blast = '0';

#my %running = MCE::Child->list();
#foreach my $bi (keys %running)
#{
    #print OUTPUT_LOG $bi." BI\n";
#}

BLAST1: foreach my $blast_id_tmp (sort {$a <=> $b} keys %blast_files)
{
    my $input_BLAST = $blast_files{$blast_id_tmp};
    if (-s $input_BLAST)
    {
        sleep(3);
        open(INPUT_BLAST, $input_BLAST);
        
        my $prev_read_id = "";
        my $first_match_tax = "";
        my $first_match_pi = "";
        my $first_match_cov = "";
        my $first_match_line = "";
        my $second_match_line = "";
        my $first_match_read_count = '1';
        my %equal_matches;
        undef %equal_matches;
        my %equal_matches_identical;
        undef %equal_matches_identical;
        my $equal_matches = '0';
        my $best_hit_count = '0';
        my $first_match_count = '1';
        my $ASV_matches = '1';

INPUT_BLAST_NP: while (my $line2 = <INPUT_BLAST>)
        {                                                     
            chomp($line2);
            my @line = split /\t/, $line2;
            my $read_id = $line[0];
            my @tax_tmp = split /;/, $line[2];
            my %tax_tmp;
            undef %tax_tmp;
            my $tax_tmp_count = '0';
            foreach my $tax_tmp (@tax_tmp)
            {
                $tax_tmp{$tax_tmp} = undef;
                $tax_tmp_count++;
            }
            my $tax_check = "";
            
            if (exists($tax_tmp{$first_match_tax}) && $tax_tmp_count < 2)
            {
                $tax_check = "yes";
            }          
print OUTPUT_LOG $read_id." READ_ID000 ".$line[3]." PI0000\n";
            if ($prev_read_id eq "")
            {
                $first_match_read_count = '1';
                $ASV_matches++;
                if (exists($identical_reads{$read_id}))
                {
                    $first_match_read_count = $identical_reads{$read_id};
                }
                
                $first_match_pi = $line[3];
                $first_match_cov = $line[4];
                $first_match_line = $line[4]."\t".$line[5]."\t".$line[2];
                my @tax_ids_tmp = split /;/, $line[2];
                foreach my $tax_id_tmp (@tax_ids_tmp)
                {
                    $equal_matches{$tax_id_tmp} = undef;
                    $equal_matches_identical{$tax_id_tmp} = undef;
                    $first_match_tax = $tax_id_tmp;
                    $equal_matches++;
                }
                $prev_read_id = $read_id;
                next INPUT_BLAST_NP;
            }
            elsif ($prev_read_id ne $read_id || eof)
            {
                print OUTPUT_LOG $prev_read_id." READ_ID0\n";
                if ($equal_matches > 1)
                {
                    $equal_matches = '0';
                    my %parent_ids;
                    undef %parent_ids;
                    my $count_matches = keys %equal_matches;
                    print OUTPUT_LOG $prev_read_id." READ_ID ".$count_matches." COUNT_MATCHES\n";

                    my %equal_matches2;
                    undef %equal_matches2;
                    
REDUCE_LEVEL:       foreach my $tax_id_tmp (keys %equal_matches)
                    {      
                        if (exists($tax_db{$tax_id_tmp}))
                        {
                            my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                            my $tax_name = $tax_db[1];

                            if (exists($taxonomy_ranks{$tax_db[2]}))
                            {
                                $equal_matches2{$taxonomy_ranks{$tax_db[2]}}{$tax_id_tmp} += 1;
                            }
                            else
                            {
                                if (exists($tax_db{$tax_db[0]}))
                                {
                                    my @tax_db2 = split /;/, $tax_db{$tax_db[0]};
                                    my $tax_name2 = $tax_db[1];
                                    if (exists($taxonomy_ranks{$tax_db2[2]}))
                                    {
                                        $equal_matches2{$taxonomy_ranks{$tax_db2[2]}}{$tax_db[0]} += 1;
                                    }
                                }
                                else
                                {
                                    print OUTPUT_LOG $tax_db[2]." TAX_ID_missing1\n";
                                    print OUTPUT_LOG $tax_id_tmp." TAX_ID_missing1b\n";
                                    #$equal_matches2{'0'}{$tax_id_tmp} += 1;
                                }
                            }
                        }
                        else
                        {
                            print OUTPUT_LOG $tax_id_tmp." TAX_ID_missing2\n";
                        }
                    }
                    
                    foreach my $tax_rank_tmp (sort {$a <=> $b} keys %equal_matches2)
                    {
                        print OUTPUT_LOG $tax_rank_tmp." TAX_RANK0\n";
                        my $count_lowest_rank_tmp = keys %{$equal_matches2{$tax_rank_tmp}};
                        
                        if ($count_lowest_rank_tmp eq '1')
                        {
                            foreach my $tax_id_tmp (keys %{$equal_matches2{$tax_rank_tmp}})
                            {
                                if (exists($tax_db{$tax_id_tmp}))
                                {
                                    my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                                    $first_match_tax = $tax_id_tmp;
                                }
                                goto COMMON_LEVEL;
                            }
                        }
                    }
                    
                    
                    my $lowest_rank = "";
                    foreach my $tax_rank_tmp (sort {$a <=> $b} keys %equal_matches2)
                    {
                        $lowest_rank = $tax_rank_tmp;
                        print OUTPUT_LOG $tax_rank_tmp." TAX_RANK\n";
                        foreach my $tax_id_tmp (keys %{$equal_matches2{$tax_rank_tmp}})
                        {
                            print OUTPUT_LOG $tax_id_tmp." TAX_ID\n";
                            if (exists($tax_db{$tax_id_tmp}))
                            {
                                my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                                if (exists($parent_ids{$tax_db[0]}))
                                {
                                    my $count_tmp = $parent_ids{$tax_db[0]};
                                    $parent_ids{$tax_db[0]} = $count_tmp+$equal_matches2{$tax_rank_tmp}{$tax_id_tmp};
                                }
                                else
                                {
                                    $parent_ids{$tax_db[0]} = $equal_matches2{$tax_rank_tmp}{$tax_id_tmp};
                                }
                            }
                            elsif ($tax_id_tmp ne '0')
                            {
                                print OUTPUT_LOG $tax_id_tmp." DB_MISSING\n";
                                goto COMMON_LEVEL;
                            }
                            else
                            {
                                print OUTPUT_LOG $tax_id_tmp." DB_MISSING1\n";
                                goto COMMON_LEVEL;
                            }
                            delete $equal_matches2{$tax_rank_tmp}{$tax_id_tmp};
                        }
                        delete $equal_matches2{$tax_rank_tmp};
                        last;
                    }

                    my $count_tmp = keys %equal_matches2; 
                    foreach my $tax_id_tmp (keys %parent_ids)
                    {
                        print OUTPUT_LOG $tax_id_tmp." TAX\n";
                        print OUTPUT_LOG $parent_ids{$tax_id_tmp}." TAX_COUNT\n";
                        if ($parent_ids{$tax_id_tmp} eq $count_matches)
                        {
                            if (exists($tax_db{$tax_id_tmp}))
                            {
                                my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                                $first_match_tax = $tax_id_tmp;
                            }
                            if ($count_tmp eq '0')
                            {
                                goto COMMON_LEVEL;
                            }
                        }
                    }
                    
                    undef %equal_matches;
                    foreach my $tax_rank_tmp (sort {$a <=> $b} keys %equal_matches2)
                    {
                        if ($tax_rank_tmp ne $lowest_rank)
                        {
                            foreach my $tax_id_tmp (keys %{$equal_matches2{$tax_rank_tmp}})
                            {
                                $equal_matches{$tax_id_tmp} += 1;
                                print OUTPUT_LOG $tax_id_tmp." ADD_TAX\n";
                            }
                        }
                    }
                    undef %equal_matches2;

                    foreach my $tax_id_tmp (keys %parent_ids)
                    {
                        $equal_matches{$tax_id_tmp} += 1;
                        print OUTPUT_LOG $tax_id_tmp." ADD_TAX2\n";
                    }
                    $count_matches = keys %equal_matches;
                    undef %parent_ids;
                    if ($count_matches eq '0')
                    {
                        goto COMMON_LEVEL;
                    }
                    goto REDUCE_LEVEL;
                }
COMMON_LEVEL:               
                if (($first_match_read_count > 1 || $taxonomy_only ne "") && $prev_read_id ne "")
                {        
                    my $identical_list = "";
                    if (keys %equal_matches_identical > 1)
                    {
                        foreach my $tax_id_tmp (keys %equal_matches_identical)
                        {
                            if (exists($tax_db{$tax_id_tmp}))
                            {
                                my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                                my $tax_name_tmp = $tax_db[1];
                                if ($identical_list ne "")
                                {
                                    $identical_list .= ";";
                                }
                                $identical_list .= $tax_name_tmp;
                            }
                        }
                        print OUTPUT_LOG $identical_list." ID_LIST0\n";
                    }
                    undef %equal_matches_identical;
                    $ASV_table{$first_match_tax}{$first_match_pi}{$first_match_read_count}{$prev_read_id}{$identical_list} = $first_match_line."*".$second_match_line."*".$first_match_count."\\".$ASV_matches;
                    $second_match_line = "";
                    $first_match_line = "";
                    $first_match_count = '1';
                }
                
                $first_match_read_count = '1';
                if (exists($identical_reads{$read_id}))
                {
                    $first_match_read_count = $identical_reads{$read_id};
                }
                #print ASV_TABLE $line[0]."\t".$line[3]."\t".$line[6]."\t".$line[19]."\t".$read_count."\t".$line[21]."\n";
                
                $first_match_pi = $line[3];
                $first_match_line = $line[4]."\t".$line[5]."\t".$line[2];
                $equal_matches = '0';
                $ASV_matches = '1';
                $first_match_tax = "";
                undef %equal_matches;

                my @tax_ids_tmp = split /;/, $line[2];
                my $second_tax_tmp = "";
                foreach my $tax_id_tmp (@tax_ids_tmp)
                {
                    $equal_matches{$tax_id_tmp} = undef;
                    $equal_matches_identical{$tax_id_tmp} = undef;
                    $equal_matches++;
                    if ($first_match_tax eq "")
                    {
                        $first_match_tax = $tax_id_tmp;
                    }
                    elsif ($second_tax_tmp eq "" && $second_tax_tmp ne $first_match_tax)
                    {
                        $second_tax_tmp = $tax_id_tmp;
                        $second_match_line = $second_tax_tmp."\t".$line[3]."\t".$line[4];
                        print OUTPUT_LOG $read_id." READ_ID1 ".$second_match_line." SECOND0000\n";
                    }
                }
                $prev_read_id = $read_id;
                next INPUT_BLAST_NP;
            }
            elsif ($line[3] > $first_match_pi && $line[3] > 98.5 && $line[4] > 0.8*$first_match_cov)
            {
                print OUTPUT_LOG $line[3]."\t".$line[3]."\t".$line[4]."\t".$line[5]." NEW_FIRST_HIT\n";               

                $first_match_tax = $line[2]; 
                $first_match_pi = $line[3];
                $first_match_cov = $line[4];
                $first_match_line = $line[4]."\t".$line[5]."\t".$line[2];
                $equal_matches = '0';
                $first_match_count = '1';
                $ASV_matches++;
                undef %equal_matches;
                undef %equal_matches_identical;
                my @tax_ids_tmp = split /;/, $line[2];
                foreach my $tax_id_tmp (@tax_ids_tmp)
                {
                    if (exists($equal_matches{$tax_id_tmp}))
                    {}
                    else
                    {
                        $equal_matches{$tax_id_tmp} = undef;
                        $equal_matches++;
                        if ($line[3] >= $first_match_pi)
                        {
                            $equal_matches_identical{$tax_id_tmp} = undef;
                        }
                    } 
                }
            }
            elsif ($tax_check ne "yes" && $line[3] > $first_match_pi-1 && ($first_match_pi < 100 || $line[3] > 99.9) && $line[4] > 0.75*$first_match_cov)
            {                
                my @tax_ids_tmp = split /;/, $line[2];
                foreach my $tax_id_tmp (@tax_ids_tmp)
                {
                    if (exists($equal_matches{$tax_id_tmp}))
                    {}
                    else
                    {
                        $equal_matches{$tax_id_tmp} = undef;
                        $equal_matches++;
                        if ($line[3] >= $first_match_pi)
                        {
                            $equal_matches_identical{$tax_id_tmp} = undef;
                        }
                    } 
                }
                $ASV_matches++;
                print OUTPUT_LOG $read_id." READ_ID1 ".$line[3]." PI\n";
            }
            else
            {
                $ASV_matches++;
            }
            if ($tax_check ne "yes" && $second_match_line eq "")
            {
                my @tax_ids_tmp = split /;/, $line[2];
                my $second_tax_tmp = "";
                foreach my $tax_id_tmp (@tax_ids_tmp)
                {
                    if ($tax_id_tmp ne $first_match_tax && $second_tax_tmp eq "")
                    {
                        $second_tax_tmp = $tax_id_tmp;
                    }
                }
                $second_match_line = $second_tax_tmp."\t".$line[3]."\t".$line[4];
                print OUTPUT_LOG $read_id." READ_ID1 ".$second_match_line." SECOND\n";
            }
            elsif ($tax_check eq "yes")
            {
                $first_match_count++;
            }
            $prev_read_id = $read_id;
        }
#Process last blast hit of file---------        
        if ($equal_matches > 100000000)
        {
            $equal_matches = '0';
            my %parent_ids;
            undef %parent_ids;
            my $count_matches = keys %equal_matches;
            my $missing_ids = '0';
           
REDUCE_LEVEL_FINAL: foreach my $tax_id_tmp (keys %equal_matches)
            {
                if (exists($tax_db{$tax_id_tmp}))
                {
                    my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                    if (exists($parent_ids{$tax_db[0]}))
                    {
                        my $count_tmp = $parent_ids{$tax_db[0]};
                        $parent_ids{$tax_db[0]} = $count_tmp+1;
                    }
                    else
                    {
                        $parent_ids{$tax_db[0]} = 1;
                    }
                }
                elsif ($tax_id_tmp ne '0')
                {
                    print OUTPUT_LOG $tax_id_tmp." DB_MISSING\n";
                    $count_matches--;
                }
                else
                {
                    print OUTPUT_LOG $tax_id_tmp." DB_MISSING1\n";
                    goto COMMON_LEVEL_FINAL;
                }

            }
            foreach my $tax_id_tmp (keys %parent_ids)
            {
                if ($parent_ids{$tax_id_tmp} eq $count_matches)
                {
                    if (exists($tax_db{$tax_id_tmp}))
                    {
                        my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                        $first_match_tax = $tax_id_tmp;
                    }
                    goto COMMON_LEVEL_FINAL;
                }
            }
            undef %equal_matches;
            foreach my $tax_id_tmp (keys %parent_ids)
            {
                $equal_matches{$tax_id_tmp} = undef;
            }
            $count_matches = keys %equal_matches;
            if ($missing_ids eq $count_matches || $count_matches eq '0')
            {
                goto COMMON_LEVEL_FINAL;
            }
            goto REDUCE_LEVEL_FINAL;

COMMON_LEVEL_FINAL:               
            if (($first_match_read_count > 1 || $taxonomy_only ne "") && $prev_read_id ne "")
            {
                my $identical_list = "";
                if (keys %equal_matches_identical > 1)
                {
                    foreach my $tax_id_tmp (keys %equal_matches_identical)
                    {
                        if (exists($tax_db{$tax_id_tmp}))
                        {
                            my @tax_db = split /;/, $tax_db{$tax_id_tmp};
                            my $tax_name_tmp = $tax_db[1];
                            if ($identical_list ne "")
                            {
                                $identical_list .= ";";
                            }
                            $identical_list .= $tax_name_tmp;
                        }
                    }
                    print OUTPUT_LOG $identical_list." ID_LIST\n";
                }
                undef %equal_matches_identical;
                $ASV_table{$first_match_tax}{$first_match_pi}{$first_match_read_count}{$prev_read_id}{$identical_list} = $first_match_line."*".$second_match_line."*".$first_match_count."\\".$ASV_matches;
            }
        }
 #---------------------------------------------------               
        close INPUT_BLAST;
        delete $blast_files{$blast_id_tmp};
        print OUTPUT_LOG $count_remaining_blast." REMAINING_BLAST\n";
    }
    else
    {
        $count_remaining_blast++;
        next BLAST1;
    }     
}

if ($count_remaining_blast > 0)
{
    goto BLAST0;
}

my %ASV_table2;
undef %ASV_table2;

ASV_table1: foreach my $species_tmp (keys %ASV_table)
{
    my $total_read_count = '0';
    my $total_pi = "";
    my $merged_read_number = "";
    my $highest_pi = "";
    my $last_pi = "";
    my $extra_info = "";
    my $last_id = "";
    my $last_identical_list = "";
    
    print OUTPUT_LOG "\n--------------\n".$species_tmp."\n-------------\n";
    
    foreach my $perc_id_tmp (sort {$b <=> $a} keys %{$ASV_table{$species_tmp}})    
    {        
        print OUTPUT_LOG $perc_id_tmp." PI\n";
        if ($highest_pi eq "")
        {
            $highest_pi = $perc_id_tmp;
        }

        foreach my $read_count_tmp (sort {$b <=>$a} keys %{$ASV_table{$species_tmp}{$perc_id_tmp}})
        {                                
            print OUTPUT_LOG $read_count_tmp." COUNT\n";
            foreach my $id_tmp (keys %{$ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}})
            {
                foreach my $identical_list_tmp (keys %{$ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}})
                {
                    print OUTPUT_LOG $id_tmp." ID_TMP\n";
                    if ($direct_ASV_output ne "")
                    {
                        $ASV_table2{$read_count_tmp}{$species_tmp}{$perc_id_tmp}{$id_tmp}{$identical_list_tmp} = $ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}{$identical_list_tmp};
                    }
                    elsif ($perc_id_tmp > 98.5 || $perc_id_tmp > $highest_pi - 0.6)
                    {
                        $total_read_count += $read_count_tmp;
                        $merged_read_number = "yes";
                        if ($perc_id_tmp ne $last_pi)
                        {
                            $total_pi .= $perc_id_tmp.";";
                        }
                        $extra_info .= $ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}{$identical_list_tmp}.";";
                    }
                    elsif ($merged_read_number eq "yes")
                    {
                        chop($total_pi);
                        chop($extra_info);
                        $ASV_table2{$total_read_count}{$species_tmp}{$total_pi}{$id_tmp}{$identical_list_tmp} = $ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}{$identical_list_tmp};
                        $merged_read_number = "";
                        $total_read_count = $read_count_tmp;
                        $total_pi = $perc_id_tmp.";";
                        $extra_info = $ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}{$identical_list_tmp}.";";
                        $highest_pi = $perc_id_tmp;
                    }
                    if ((($read_count_tmp > 1 && $perc_id_tmp > 98.5) || $read_count_tmp > 5 || $taxonomy_only ne "") && $merged_read_number eq "" && $direct_ASV_output eq "")
                    {
                        $merged_read_number = "yes";
                        $total_read_count = $read_count_tmp;
                        $total_pi = $perc_id_tmp.";";
                        $extra_info = $ASV_table{$species_tmp}{$perc_id_tmp}{$read_count_tmp}{$id_tmp}{$identical_list_tmp}.";";
                        $highest_pi = $perc_id_tmp;
                    }
                    $last_pi = $perc_id_tmp;
                    $last_id = $id_tmp;
                    $last_identical_list = $identical_list_tmp;
                    #my $read_seq = $ASV_table{$species_tmp}{$read_count_tmp}{$perc_id_tmp}{$id_tmp};
                }
            }
        }     
    }
    if ($merged_read_number eq "yes")
    {
        chop($total_pi);
        chop($extra_info);
        $ASV_table2{$total_read_count}{$species_tmp}{$total_pi}{$last_id}{$last_identical_list} = $extra_info;
        $merged_read_number = "";
        $total_read_count = '0';
        $total_pi = "";
    }
}

if ($direct_ASV_output eq "")
{
    print ASV_TABLE "Read count\tKingdom\tPhylum\tClass\tFamily\tOrder\tGenus\tSpecies\tPI\tAdditional info\n";
}
else
{
    print ASV_TABLE "ID\tKingdom\tPhylum\tClass\tFamily\tOrder\tGenus\tSpecies\tBest Hit\tBest Hit count\tPI%\tQC%\tRead count\tSecond best\t2nd PI%\t2nd QC%\tBest hit 97%\tFirst Alternative Hit\tSecond Alternative Hit\tList identical matches\tSequence\n";
}
foreach my $read_count_tmp (sort {$b <=> $a} keys %ASV_table2)
{
    foreach my $tax_tmp0 (sort {$a <=> $b} keys %{$ASV_table2{$read_count_tmp}})
    {
        my $total_perc_id = "";
        my $extra_info = "";
        my $tax_tmp = $tax_tmp0;
        if ($direct_ASV_output eq "")
        {
            foreach my $perc_id_tmp (sort {$b <=> $a} keys %{$ASV_table2{$read_count_tmp}{$tax_tmp}})
            {
                $total_perc_id .= $perc_id_tmp.";";
                $extra_info .= $ASV_table2{$read_count_tmp}{$tax_tmp}{$perc_id_tmp}.";";
            }
        }
        
        chop($total_perc_id);
        chop($extra_info);
        my $species_tmp = "";
        my $genus_tmp = "";
        my $order_tmp = "";
        my $family_tmp = "";
        my $class_tmp = "";
        my $phylum_tmp = "";
        my $kingdom_tmp = "";
DROP_RANK:        
        if (exists($tax_db{$tax_tmp}))
        {
            my @tax_db = split /;/, $tax_db{$tax_tmp};
            if ($tax_db[2] eq "species")
            {
                $species_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "genus")
            {
                $genus_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "order")
            {
                $order_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "family")
            {
                $family_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "class")
            {
                $class_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "phylum")
            {
                $phylum_tmp = $tax_db[1];
            }
            elsif ($tax_db[2] eq "kingdom")
            {
                $kingdom_tmp = $tax_db[1];
            }
            if ($tax_db[0] ne "0")
            {
                $tax_tmp = $tax_db[0];
                goto DROP_RANK
            }
        }
        else
        {
            $species_tmp = $tax_tmp;
        }
        
        if ($species_tmp eq "")
        {
            $species_tmp = "unassigned";
        }
        if ($genus_tmp eq "")
        {
            $genus_tmp = "unassigned";
        }
        if ($order_tmp eq "")
        {
            $order_tmp = "unassigned";
        }
        if ($family_tmp eq "")
        {
            $family_tmp = "unassigned";
        }
        if ($class_tmp eq "")
        {
            $class_tmp = "unassigned";
        }
        if ($phylum_tmp eq "")
        {
            $phylum_tmp = "unassigned";
        }
        if ($kingdom_tmp eq "")
        {
            $kingdom_tmp = "unassigned";
        }
        if ($direct_ASV_output eq "")
        {
            print ASV_TABLE $read_count_tmp."\t".$kingdom_tmp."\t".$phylum_tmp."\t".$class_tmp."\t".$family_tmp."\t".$order_tmp."\t".$genus_tmp."\t".$species_tmp."\t".$total_perc_id."\t".$extra_info."\n";
        }
        else
        {
            foreach my $perc_id_tmp (sort {$b <=> $a} keys %{$ASV_table2{$read_count_tmp}{$tax_tmp0}})
            {
                foreach my $id_tmp (sort {$a <=> $b} keys %{$ASV_table2{$read_count_tmp}{$tax_tmp0}{$perc_id_tmp}})
                {
                    foreach my $idential_list_tmp (sort {$a <=> $b} keys %{$ASV_table2{$read_count_tmp}{$tax_tmp0}{$perc_id_tmp}{$id_tmp}})
                    {
                        my $extra_info_tmp = $ASV_table2{$read_count_tmp}{$tax_tmp0}{$perc_id_tmp}{$id_tmp}{$idential_list_tmp};
                        my @extra_info_tmp0 = split /\*/, $extra_info_tmp;
                        my @extra_info_tmp = split /;/, $extra_info_tmp0[0];
                        my @extra_info_tmp2 = split /\t/, $extra_info_tmp[0];
                        my @extra_info_tmp_2nd = split /;/, $extra_info_tmp0[1];
                        my @extra_info_tmp2_2nd = split /\t/, $extra_info_tmp_2nd[0];
                        my $qc = $extra_info_tmp2[0];
                        my $best_hit_tax_id = $extra_info_tmp2[2];
                        my $original_seq = $original_ASVs{$id_tmp};
                        my $secondhit = $extra_info_tmp2_2nd[0];
                        my $secondpi = $extra_info_tmp2_2nd[1];
                        my $secondqc = $extra_info_tmp2_2nd[2];
                        my $first_hit_count = $extra_info_tmp0[2];
                                          
                        my $best_hit_species = "";
                        if (exists($tax_db{$best_hit_tax_id}))
                        {
                            my @tax_db = split /;/, $tax_db{$best_hit_tax_id};
                            $best_hit_species = $tax_db[1];      
                        }
                        else
                        {
                            $best_hit_species = $best_hit_tax_id;
                        }
             
                        my $second_best_hit_species = "";
                        if (exists($tax_db{$secondhit}))
                        {
                            my @tax_db = split /;/, $tax_db{$secondhit};
                            $second_best_hit_species = $tax_db[1];      
                        }
                        else
                        {
                            $second_best_hit_species = $secondhit;
                        }
                        
                        my $best_hit_97 = "";
                        my $first_alt_hit_97 = "";
                        my $second_alt_hit_97 = "";
                        if ($perc_id_tmp >= 97 && $idential_list_tmp eq "")
                        {
                            $best_hit_97 = $best_hit_species;
                        }
                        elsif ($idential_list_tmp eq "")
                        {
                            $first_alt_hit_97 = $best_hit_species;
                            $second_alt_hit_97 = $second_best_hit_species;
                        }
                        
                        print ASV_TABLE $id_tmp."\t".$kingdom_tmp."\t".$phylum_tmp."\t".$class_tmp."\t".$family_tmp."\t".$order_tmp."\t".$genus_tmp."\t".$species_tmp."\t".$best_hit_species."\t".$first_hit_count
                        ."\t".$perc_id_tmp."\t".$qc."\t".$read_count_tmp."\t".$second_best_hit_species."\t".$secondpi."\t".$secondqc."\t".$best_hit_97."\t".$first_alt_hit_97."\t".$second_alt_hit_97."\t".$idential_list_tmp."\t".$original_seq."\n";
                        delete $blast_files_id{$id_tmp};
                    }
                }
            }
        }
    }
}
foreach my $id_tmp2 (keys %blast_files_id)
{
     print ASV_TABLE $id_tmp2."\t\t\t\t\t\t\tNO BLAST HIT\n";
}
END1:
print "\n";
if ($taxonomy_only eq "")
{
    print "Non merged reads    : ".$no_merge."\n";
    print "Unmatched IDs       : ".$no_id_match."\n";
}
    print "Identical reads     : ".$identical_reads."\n";
    print "Total reads         : ".$keys_hash."\n\n";

close $FILE;
if ($file_count eq "2")
{
    close $FILE2;
}
close ASV_TABLE;
close OUTPUT_LOG;

$chnl->end;
MCE::Child->waitall;

print "\nThank you for using Taxy!\n\n";
