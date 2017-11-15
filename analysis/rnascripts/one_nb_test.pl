#!/sw/bin/perl

=head1 NAME

B<one_nb_test.pl>

=head1 DESCRIPTION

This scrip is a part of C<grid_launcher_nbtest.pl>.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;


use PDL;
use PDL::NiceSlice;

use RNASeq;
use RGetOpt;

use Tools;
use Distribution;

$| = 1;

my $created = "Created with $0 " . join(" ", @ARGV) . "\n";

ReadDefs();

my ($ncrit, $maxn);
my $batchsize = 30;
my $stat = 'm';
my $test = 'mi';
my $batch = 0;
my $outfile;
my ($help, $man);
GetMyOptions(
  'batch=i' => \$batch,
  'test=s' => \$test,
  'stat=s' => \$stat,
  'outfile=s' => \$outfile,
  'batchsize=i' => \$batchsize,
  'ncrit=i' => \$ncrit,
  'maxn=f' => \$maxn
);

my $hostname = $ENV{HOSTNAME};

print "$created\n";
print "Batch $batch\n";
print "Outfile $outfile\n";
print "Host $hostname\n" if defined $hostname;
print "ncrit $ncrit\n" if defined $ncrit;
print "maxn $maxn\n" if defined $maxn;

my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = GetNormalizedData(
  condition=>$cond,
  replicate=>$rep,
  exclude=>$exclude,
  nonzero=>1,
  norm=>$norm,
  type=>$type,
  clean=>$clean,
  spclean=>$spclean
);


####################################

#$batch++ if $batch == 0;

my $n1 = ($batch - 1) * $batchsize;
my $n2 = $n1 + $batchsize - 1;
$n2 = $ngen - 1 if $n2 > $ngen - 1;
exit if $n2 < $n1;

print "Processing batch $batch, $n1 .. $n2 \(n = $ngen)\n";

local *F;
open F, ">$outfile" or die;


for my $i ($n1 .. $n2)
{
  my $f = "h${hostname}_b${batch}_i${i}";
  my $gene = $genes->[at($sel, $i)];
  my @out = ($gene);
  my $x = $d(,$i;-);
  my ($m, $s) = stats($x);
  push @out, sprintf "%.3g", $m;
  push @out, sprintf "%.3g", $s;
  my $P;
  if($test eq 'meintanis'){
      print "Gene: $gene\tM: $m\tS: $s\tX: $x\tNcrit: $ncrit\tMaxN: $maxn\tStat: $stat\n";
      $P = NBBootstrapTest($stat, $x, $ncrit, $maxn)
  }
  elsif($test eq 'mi')
    {$P = NBMiTest($x, nsim=>$maxn, fname=>$f, debug=>1)}
  #print "$gene\t$m\t$s\t$P\n";
  printf F "%s\t%.4g\t%.4g\t%.4g\n", $gene, $m, $s, $P;
  Percent(($i-$n1)/($n2-$n1));
}
print "\n";

close F
