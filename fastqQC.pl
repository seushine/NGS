#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;

my $usage = <<"USAGE";

        FastQC Sequence QC

Usage:
        perl $0 [opts]
        
Options:
    --- Input/Output reads/sequences (FASTQ) ---
        -a | -in1 <Forward read/sequence file>
        -b | -in2 <Reverse read/sequence file of paired-end data>
        -o <Trimmed sequence file>
    --- Trimming Options ---
        -l | -leftTrimBases <Integer>
        Number of bases to be trimmed from left end (5' end)
        default: 0
        -r | -rightTrimBases <Integer>
        Number of bases to be trimmed from right end (3' end)
        default: 0
        -min | -minlenCutoff <Integer>
        Reads shorter than given length will be discarded
        default: 0 (i.e. length filtering is OFF)
        -max | -maxlenCutoff <Integer>
        Reads longer than given length will be trimed
        default: 0 (i.e. length filtering is OFF)
        -qt | qualTrimming <Bool>
        Trim bases having PHRED quality score less than [QualCutOff] at 3' end of the read
        default: False
    --- Filter Options ---
        -q | -qualCutOff <Integer>
        Cut-off PHRED quality score for filtering reads
        -p | -percentage <Integer>
        Cut-off PHRED quality percentageuency for filtering reads
        For eg.: -q 20 -f 80, will discard read having 80% bases PHRED quality score less than 20
    --- Other options ---
        -illu | -illumina <Bool>
        PHRED quality score Q64 if True
        -z | -gz <Bool> 
        Output sequence file will be compressed
        -s | -show <Bool>
        Show Running Details
        -h | -help <Bool>
        <Prints this help>

Example:
        perl $0 -in exam.fq -o exam_trimmed.fq -l 5 -r 10 -min 30 -max 200 -q 20 -p 80

Version:
        V0.2

Date:
        20161024 <The Programmer's Day>

USAGE

my ($infile,$infile2,$outfile,$qualTrimming,$illumina,$gzip,$showDetails,$sampleName,$help);
my $rTrimBases = 0;
my $lTrimBases = 0;
my $minlenCutoff = 0;
my $maxlenCutoff = 0;
my $qualCutoff = 0;
my $percentage = 0;
GetOptions(
        'a|in1=s'    => \$infile,
        'b|in2=s'   => \$infile2,
        'o|outfile=s'=> \$outfile,
        'l|leftTrimBases=i'  => \$lTrimBases,
        'r|rightTrimBases=i' => \$rTrimBases,
        'min|minlenCutoff=i' => \$minlenCutoff,
        'max|maxlenCutoff=i' => \$maxlenCutoff,
        'qt|qualTrimming?' => \$qualTrimming,
        'q|qualCutoff=i' => \$qualCutoff,
        'p|percentage=s'   => \$percentage,
        'illu|illumina?'=>\$illumina,
        'z|gz?' => \$gzip,
        's|show?' => \$showDetails,
        'n|name=s' => \$sampleName,
        'h|help?' => \$help,
        );

die $usage if $help;
die "The Input Parameters ERROR!\n\n$usage" if (!$infile||!$outfile);

###### Opts Detection

open LOG, ">$outfile.log" or die "can't write logfile\n";
my $subVal = 33;
if ($illumina){
    $subVal = 64;
    $illumina = "Q64";
}else{
    $illumina = "Q33";
}
if ($qualTrimming){
    $qualTrimming = "Ture";
}else{
    $qualTrimming = "False";
}
if (!defined($sampleName)){
    $sampleName = basename($outfile);
}

my $runDetails = <<"INFO";
Software Running Info:
Sample Name:    $sampleName
Input:    $infile
Output:   $outfile
Trim Parameters:
leftTrimBases:  $lTrimBases
rightTrimBases: $rTrimBases
minlenCutoff:   $minlenCutoff
maxlenCutoff:   $maxlenCutoff
qualTrimming:   $qualTrimming
Filter Parameters:
Quality:  $qualCutoff
Percentage:  $percentage
Other Parameters:
Quality Format: $illumina
INFO
print LOG "$runDetails\n";
if ($showDetails){
    print "$runDetails\n";
}
######
die "Parameters ERROR!\nTrim Bases must between 0bp and 100bp\nUse -h to get more infomations\n\n" if ($lTrimBases<0 or $lTrimBases>100 or $rTrimBases<0 or $rTrimBases>100);
die "Parameters ERROR!\nMinReadsLengthCutOff less than 100bp\nMaxReadsLengthCutOff larger than 22bp\n\n" if ($minlenCutoff>100 or $maxlenCutoff<22);
die "Parameters ERROR!\nQuality Cutoff must between 0 and 40\n\n" if ($qualCutoff<0 or $qualCutoff>40);
die "Parameters ERROR!\nQuality Filter Percentage must between 0 and 100\n\n" if ($percentage<0 or $percentage>100);

######
my $totalLines = 0;
if ($showDetails){
    if (-B $infile){$totalLines = readpipe("zcat $infile|wc -l");}else{$totalLines = readpipe("wc -l $infile");}
    my @temp = split(/\s+/,$totalLines);
    $totalLines = $temp[0]/4;
}
my $temp_progress = 0;

######
if (-B $infile){open IN, "<:gzip","$infile" or die "$!\n";}else{open IN, "$infile" or die "$!\n";}
if ($gzip){open OUT, ">:gzip","$outfile" or die "can't write file: $outfile\n";}
else{open OUT, ">$outfile" or die "can't write file: $outfile\n";}

my $total_reads;
my ($count_base,$count_gc,$count_n,$count_q20,$count_q30) = (0,0,0,0,0);
my $pass_reads;
my ($count_pass_base,$count_pass_gc,$count_pass_n,$count_pass_q20,$count_pass_q30) = (0,0,0,0,0);

while(my $id=<IN>){
    chomp($id);
    chomp(my $seq = <IN>);
    chomp(my $plus = <IN>);
    chomp(my $ascii = <IN>);
    my $pass = "False";
    $total_reads++;
    my $length = length $seq;
    my $length_trim = $length;
    my ($base_gc,$base_n,$q20_base,$q30_base) = &statistic($seq,$ascii);
    $count_base += $length;
    $count_gc += $base_gc;
    $count_n += $base_n;
    $count_q20 += $q20_base;
    $count_q30 += $q30_base;

    if ($lTrimBases!=0 and $lTrimBases<$length_trim){
        $seq = substr($seq,$lTrimBases);
        $ascii = substr($ascii,$lTrimBases);
        $length_trim = length $seq;
    }
    if ($rTrimBases!=0 and $rTrimBases<$length){
        $seq = substr($seq,0,$length_trim-$rTrimBases);
        $ascii = substr($ascii,0,$length_trim-$rTrimBases);
        $length_trim = length $seq;
    }
    next if ($length_trim < $minlenCutoff);
    if ($maxlenCutoff != 0 and $length_trim > $maxlenCutoff){
        $seq = substr($seq,0,$maxlenCutoff);
        $ascii = substr($ascii,0,$maxlenCutoff);
        $length_trim = length $seq;
    }

    if ($qualCutoff == 0 or $percentage == 0){
        $pass = "True";
    }
    else{
        if ($qualTrimming eq "True" and $qualCutoff != 0){
            ($seq,$ascii) = &trimSeq4Qual($seq,$ascii);
        }
        my @ASCII = unpack("C*", $ascii);
        my $pass_qual_base = 0;
        foreach my $qbase (@ASCII){
            $qbase -= $subVal;
            $pass_qual_base++ if $qbase>=$qualCutoff;
        }
        if ($pass_qual_base/$length_trim*100 >= $percentage){
            $pass = "True";
        }
    }
    if ($pass eq "True"){
        $pass_reads++;
        print OUT "$id\n$seq\n$plus\n$ascii\n";
        ($base_gc,$base_n,$q20_base,$q30_base) = &statistic($seq,$ascii);
        $count_pass_base += $length_trim;
        $count_pass_gc += $base_gc;
        $count_pass_n += $base_n;
        $count_pass_q20 += $q20_base;
        $count_pass_q30 += $q30_base;
    }

    if ($showDetails){
        if (100*$total_reads/$totalLines > $temp_progress){
            print "Process approx $temp_progress% complete...\n";
            $temp_progress+=10;
        }
    }
}
if ($showDetails){
    print "Process 100% complete!\n";
}

my $pass_ratio = sprintf("%.2f",100*$pass_reads/$total_reads);
my $pass_base_ratio = sprintf("%.2f",100*$count_pass_base/$count_base);
print LOG "Input: $total_reads\t$count_base\nOutput: $pass_reads\t$count_pass_base\nCleanRatio: $pass_ratio\t$pass_base_ratio\n";

my $gc_rate = sprintf("%.2f",100*$count_gc/($count_base-$count_n));
my $q20_rate = sprintf("%.2f",100*$count_q20/$count_base);
my $q30_rate = sprintf("%.2f",100*$count_q30/$count_base);
my $gc_rate_pass = sprintf("%.2f",100*$count_pass_gc/($count_pass_base-$count_pass_n));
my $q20_rate_pass = sprintf("%.2f",100*$count_pass_q20/$count_pass_base);
my $q30_rate_pass = sprintf("%.2f",100*$count_pass_q30/$count_pass_base);
open XLS, ">$outfile.xls" or die "can't write statistics file\n";
print XLS "$sampleName\t$total_reads\t$count_base\t$gc_rate\t$q20_rate\t$q30_rate\t$pass_reads\t$pass_ratio\t$count_pass_base\t$pass_base_ratio\t$gc_rate_pass\t$q20_rate_pass\t$q30_rate_pass\n";

sub statistic{
    my ($seq,$qual) = @_;
    my $base_gc += $seq =~ tr/gcGC/gcGC/;
    my $base_n += $seq =~ tr/nN/nN/;
    my @ASCII = unpack("C*", $qual);
    my ($q20_base,$q30_base) = (0,0);
    foreach my $qbase (@ASCII){
        $qbase -= $subVal;
        $q20_base++ if $qbase>=20;
        $q30_base++ if $qbase>=30;
    }
    return $base_gc,$base_n,$q20_base,$q30_base;
}

sub trimSeq{
    my ($seq,$leftTrim,$rightTrim) = @_;
    $seq =~ /^.{$leftTrim}(.+).{$rightTrim}$/;
    return $1;
}

sub trimSeq4Qual{
    my ($seq,$qual) = @_;
    my @ASCII = unpack("C*", $qual);
    my $trimCount = 0;
    my $Count = 0;
    my $length = length $seq;
    for(my $i=@ASCII; $i>0; $i--) {
        $Count++;
        my $val = $ASCII[$i-1] - $subVal;
        if($val < $qualCutoff) {
            $trimCount++;
        }
        else {
            if($Count>=5 and ($trimCount/$Count)<0.5){
                last;
            }
            elsif($Count>=5 and $trimCount<=2){
                $trimCount=5;
                last;
            }
            else{
                last;
            }
        }
    }
    $seq = substr($seq,0,$length-$trimCount);
    $qual = substr($qual,0,$length-$trimCount);
    return $seq,$qual;
}
