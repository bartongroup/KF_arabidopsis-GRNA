#!/sw/bin/perl -w

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;



my $use = "perl $0 -c PREFIX -r REPLICATE [-p PAIRED] [-i FASTQDIR] [-o LEVFASTQDIR] [-n COUNTSFILE]";
{
    my $l = scalar(@ARGV);
    die "Usage: $use\n" if ($l < 4 || $l % 2 != 0);
}

=head1 NAME

B<$0>   

=head1 USAGE

$use

-c PREFIX       :   The replicate prefix for the _R1.fastq.gz and _R2.fastq.gz files, 
                    both input and output.
-r REPLICATE    :   Replicate index number.
-p PAIRED       :   0 = single-end, 1 = paired-end (0).
-i FASTQDIR     :   Directory where input files are located (./fastq).
-o LEVFASTQDIR  :   Directory where output files will be written (./fastqlev).
-n COUNTSFILE   :   Text file with read counts per replicate (./combined_counts/readcount.dat).

File name convention    :   PREFIX_REPLICATE_R!.fastq.gz and PREFIX_REPLICATE_R2.fastq.gz .
Counts file format      :   PREFIX_REPLICATE  readcount

=head1 DESCRIPTION

Equal count (levelling) normalization of fastq files.

Modified from Marek's original to remove excessive and unnecessary dependencies on other custom packages:
Modified from Marek's original to remove dependency on config file.
Modified from Marek's original to enable input from any specified pair of files, instead of entire directory.
Modified from Marek's original to enable output to any file destination.
Modified from Marek's original to enable paired-end levelling.
Modified from Marek's original to use a single readcount file.
Modified from Marek's original to use a modified readcount format, to accommodate multiple conditions in the same file.
Modified from MArek's original to simplify the internal hash of counts per replicate.

=cut

# Log.
my $logscript = (getpwuid $>)[7] . "/scripts/mylogs.py";
`python $logscript "$0 @ARGV"`;


$| = 1;

# Get parameters
my ($countsfile, $paired, $fastqdir, $levdir, $cond, $rep);
my %args = @ARGV;
$countsfile = defined($args{"-n"}) ? $args{"-n"} : "./combined_counts/readcount.dat";
$paired = defined($args{"-p"}) ? $args{"-p"} : 0;
$fastqdir = defined($args{"-i"}) ? $args{"-i"} : "./fastq";
$levdir = defined($args{"-o"}) ? $args{"-o"} : "./fastqlev";
if (defined($args{"-c"})) {   $cond = $args{"-c"};    } else  {   die "No PREFIX.";  }
if (defined($args{"-r"})) {   $rep = $args{"-r"};    } else  {   die "No REPLICATE.";    }

my ($cnt, $min) = ReadCountFile();

my $count = $cnt->{$cond}{$rep};
my $fastq1 = "${cond}_${rep}_R1.fastq.gz";
my $fastq2 = "${cond}_${rep}_R2.fastq.gz";
print "${cond}_${rep} paired n=$count  ";

Level("$fastqdir/$fastq1", "$levdir/$fastq1", $count, $min, "$fastqdir/$fastq2", "$levdir/$fastq2");

print "\n";



####################################################

sub Level
{
  my ($infile1, $outfile1, $nreads, $ntarget, $infile2, $outfile2) = @_;
  
  print "randomizing... ";
  my $r = random($nreads);    # random vector
  my $perm = qsorti $r;       # random permutation of $nreads indices
  my $sel = qsort $perm(0:$ntarget-1);   # $ntarget random indices
  
  my $tmpout1 = "./leveout$$.fastq1";
  my $tmpin1= "./levein$$.fastq1";
  my $tmpout2 = "./leveout$$.fastq2";
  my $tmpin2= "./levein$$.fastq2";
    
  print "gunzipping... ";
  `gunzip $infile1 -c > $tmpin1`;
  `gunzip $infile2 -c > $tmpin2` if $paired;
  chop $infile1; chop $infile1; chop $infile1;
  chop $infile2; chop $infile2; chop $infile2;
  
  print "sampling...";
  open IN1, "$tmpin1" or die "Problem with input R1.";
  open OUT1, ">$tmpout1" or die "Problem with output R1.";
  open IN2, "$tmpin2" or die "Problem with input R2."  if $paired;
  open OUT2, ">$tmpout2" or die "Problem with output R2." if $paired;
  
  my $i = 0;
  my $n = 0;
  while($i < $ntarget - 1)
  {
    my $read1 = GetRead(*IN1);
    my $read2 = GetRead(*IN2) if $paired;
    if(at($sel, $i) == $n)
    {
      #print "$i $n\n";
      print OUT1 "$read1";
      print OUT2 "$read2" if $paired;
      $i++
    }
    $n++
  }
  
  print "gzipping... ";
  `gzip -c $tmpout1 > $outfile1`;
  `gzip -c $tmpout2 > $outfile2` if $paired;
  
  unlink $tmpout1;
  unlink $tmpin1; 
  if ($paired) {
    unlink $tmpout2;
    unlink $tmpin2;
  } 
  
}

sub GetRead
{
  my $z = shift;
  my $r = <$z>;
  die "Unexpected end of file\n" unless defined $r;
  for my $i (2 .. 4)
    {$r .= <$z>}
  return $r;
}


sub ReadCountFile
{
  my %cnt = ();
  my $min = 1e32; # Start with a huge value
  {
    local *F;
    open F, $countsfile or die;
    while(my $line = <F>)
    {
      chomp $line;
      my ($cond, $rep, $count) = split /\s+/, $line;
        $cnt{$cond}{$rep} = $count;
      $min = $count if $count < $min
    }
  }
  return \%cnt, $min
}



