package GetOpt;

=head1 NAME

GetOpt 

=head1 SYNOPSIS

  use GetOpt;
  
=head1 DESCRIPTION
  
=head1 FUNCTIONS

=over 4

=cut

require Exporter;

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use GRNASeq;

our @ISA = ("Exporter");
our @EXPORT = qw(
  $psfile
  $cond
  $rep
  $nonzero
  $norm
  $type
  $exclude
  $include
  $clip
  $genes
  $genlist
  $gff
  $clean
  $spclean
  $randrep
  $randrep2
  $multicor
  $help
  $man
  GetMyOptions
);


our $psfile;
our $cond = 'WT';
our $rep;
our $nonzero = 1;
our $norm = 'deseq';
our $type = 'raw';
our $exclude;
our $include;
our $clip;
our $genes;
our $genlist;
our $clean;
our $spclean;
our $randrep;
our $randrep2;
our $multicor = 'bh';  # Benjamini-Hochberg
our $gff;
our ($help, $man);


sub GetMyOptions
{
  my %opt = @_;

  my $stat = GetOptions(
    'psfile=s' => \$psfile,
    'cond=s' => \$cond,
    'rep=i' => \$rep,
    'nonzero=i' => \$nonzero,
    'norm=s' => \$norm,
    'type=s' => \$type,
    'exclude=s' => \$exclude,
    'include=s' => \$include,
    'clip=i' => \$clip,
    'genes=s' => \$genes,
    'genlist=s' => \$genlist,
    'randrep=i' => \$randrep,
    'randrep2=i' => \$randrep2,
    'multicor=s' => \$multicor,
    'countsdir=s' => \$countsdir,
    'outlierdir=s' => \$outlierdir,
    'pileupdir=s' => \$pileupdir,
    'gff=s' => \$gff,
    clean => \$clean,
    spclean => \$spclean,
    help => \$help,
    man => \$man,
    %opt
  );
  die "Failed getting options. Check your command-line options and try again.\n" unless $stat;
  pod2usage(-verbose => 2) if $man;
  pod2usage(-verbose => 1) if $help;
}





1;
