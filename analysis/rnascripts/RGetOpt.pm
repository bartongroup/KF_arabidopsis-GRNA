package RGetOpt;

=head1 NAME

RGetOpt 

=head1 SYNOPSIS

  use RGetOpt;
  
=head1 DESCRIPTION
  
=head1 FUNCTIONS

=over 4

=cut

require Exporter;

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

our @ISA = ("Exporter");
our @EXPORT = qw(
  $starindex 
  $subreadindex
  $bwaindex
  $gtffile
  $genomefile
  $fastqdir
  $samdir
  $rawcountsdir
  $countsdir
  $cuffdir
  $tabmdir
  $balldir
  $misodir
  $outlierdir
  $pileupdir
  $misogffdir
  @experiment
  %nrep
  @condpairs
  @conds
  @condnames
  $rep
  $cond
  $nonzero
  $norm
  $type
  $exclude
  $include
  $CleanExclude
  $SpikeinCleanExclude
  $clip
  $genes
  $genlist
  $clean
  $spclean
  $randrep
  $randrep2
  $multicor
  $psfile
  $help
  $man
  $logdir
  $stardir
  $tophatdir
  $bowtiedir
  $samtoolsdir
  $picarddir
  $subreaddir
  $bwadir
  $cufftoolsdir
  $misobin
  $expfile
  $paired
  $stranded
  GetMyOptions
  ReadDefs
  ReadExperiment
);

our ($starindex, $gtffile, $genomefile, $subreadindex, $bwaindex);
our ($fastqdir, $samdir, $rawcountsdir, $countsdir, $outlierdir, $pileupdir, $cuffdir, $tabmdir, $balldir, $misodir);
our ($logdir, $stardir, $subreaddir, $tophatdir, $bowtiedir, $samtoolsdir, $picarddir, $cufftoolsdir, $misobin, $bwadir);

our $paired;
our $stranded;
our %nrep;
our @experiment;
our @conds;
our @condnames;
our @condpairs;
our $rep;
our $misogffdir;
our $cond = 'WT';
our $nonzero = 1;
our $norm = 'deseq';
our $type = 'raw';
our $exclude;
our $include;
our $clip;
our $genes;
our $genlist = 'genlist.txt';
our $clean;
our $spclean;
our $randrep;
our $randrep2;
our $multicor = 'none';  # bh: Benjamini-Hochberg, hb: Holm-Bonferroni 
our $genome;
our $psfile;
our ($help, $man);

my $defsfile = 'defs.dat';
our $expfile;

our $topdir;
our ($CleanExclude, $SpikeinCleanExclude); 

our %defs;

my @required = qw(topdir expfile countsdir cuffdir outlierdir pileupdir bootstrapdir gtffile clean spclean conditions condnames paired logdir stardir tophatdir bowtiedir samtoolsdir cufftoolsdir);


sub GetMyOptions
{
  my %opt = @_;

  my $conds;
  my $condnames;
  my $condpairs;
  my $spaired;
  my $stat = GetOptions(
    'psfile=s' => \$psfile,
    'cond=s' => \$cond,
    'rep=i' => \$rep,
    'conds=s' => \$conds,
    'condnames=s' => \$condnames,
    'condpairs=s' => \$condpairs,
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
    'fastqdir=s' => \$fastqdir,
    'samdir=s' => \$samdir,
    'rawcountsdir=s' => \$rawcountsdir,
    'countsdir=s' => \$countsdir,
    'cuffdir=s' => \$cuffdir,
    'tabmdir=s' => \$tabmdir,
    'balldir=s' => \$balldir,
    'misodir=s' => \$misodir,
    'misogffdir=s' => \$misogffdir,
    'outlierdir=s' => \$outlierdir,
    'pileupdir=s' => \$pileupdir,
    'starindex=s' => \$starindex,
    'bwaindex=s' => \$bwaindex,
    'gtffile=s' => \$gtffile,
    'genomefile=s' => \$genomefile,
    'logdir=s' => \$logdir,
    'stardir=s' => \$stardir,
    'tophatdir=s' => \$tophatdir,
    'subreaddir=s' => \$subreaddir,
    'bowtiedir=s' => \$bowtiedir,
    'samtoolsdir=s' => \$samtoolsdir,
    'picarddir=s' => \$picarddir,
    'cufftoolsdir=s' => \$cufftoolsdir,
    'misobin=s' => \$misobin,
    'genome=s' => \$genome,
    'paired=s' => \$spaired,  # yes or no
    clean => \$clean,
    spclean => \$spclean,
    help => \$help,
    man => \$man,
    %opt
  );
  die "Failed getting options. Check your command-line options and try again.\n" unless $stat;
  pod2usage(-verbose => 2) if $man;
  pod2usage(-verbose => 1) if $help;
  
  if(defined $conds)
  {
    my %iscond = map {$_ => 1} @conds;
    my @selconds = split(/\,/, $conds);
    for my $cond (@selconds)
      {die "Condition $cond not present in defs.dat, but specified in -conds\n" unless $iscond{$cond}}
    my %cond2i = map {$conds[$_] => $_} (0 .. @conds - 1);
    @condnames = map {$condnames[$cond2i{$_}]} @selconds;
    @conds = @selconds;
  }
  @condnames = split(/\,/, $condnames) if defined $condnames;
  @condpairs = split(/,/, $condpairs) if defined $condpairs;
  if(defined $spaired)
  {
    if($spaired eq 'yes' || $spaired eq 'y')
      {$paired = 1}
    elsif($spaired eq 'no' || $spaired eq 'n')
      {$paired = 0}
  }
}

################################################################

sub ReadDefs
{
  local *F;
  open F, $defsfile or die "Cannot open $defsfile\n";
  
  %defs = ();
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    if($line =~ /^\s*(\S+)\s*=\s*(\S+)\s*$/)
      {$defs{$1} = $2;}
  }
  
  for my $r (@required)
    {die "$r not defined in $defsfile\n" unless defined $defs{$r}}
  
  @conds = split(/,/, $defs{conditions});
  my %iscond = map {$_ => 1} @conds;
  $defs{condnames} =~ s/\_/ /g;
  @condnames = split(/,/, $defs{condnames});
  
  my @nrep = split(/,/, $defs{nrep});
  die "Replicates do not correspond to conditions\n" if $#conds != $#nrep;
  %nrep = map {$conds[$_] => $nrep[$_]} (0 .. @conds - 1); 
  
  if(defined $defs{condpairs})
  {
    @condpairs = split(/,/, $defs{condpairs});
    for my $c (@condpairs)
    {
      my ($cond1, $cond2) = split(/:/, $c);
      die "Incorrect format of -condpairs\n" unless defined $cond1 && defined $cond2;
      die "Condition $cond1 not recognized in -condpairs\n" unless $iscond{$cond1};
      die "Condition $cond2 not recognized in -condpairs\n" unless $iscond{$cond2};
    }
  }
  
  $topdir = $defs{topdir};
  $expfile = "$topdir/$defs{expfile}";
  $logdir = "$topdir/$defs{logdir}";
  $rawcountsdir = "$topdir/$defs{rawcountsdir}";
  $countsdir = "$topdir/$defs{countsdir}";
  $cuffdir = "$topdir/$defs{cuffdir}";
  $tabmdir = "$topdir/$defs{tabmdir}" if defined $defs{tabmdir};
  $balldir = "$topdir/$defs{balldir}" if defined $defs{balldir};
  $misodir = "$topdir/$defs{misodir}";
  $outlierdir = "$topdir/$defs{outlierdir}";
  $pileupdir = "$topdir/$defs{pileupdir}";
  $fastqdir = "$topdir/$defs{fastqdir}";
  $samdir = "$topdir/$defs{samdir}";
  $starindex = "$topdir/$defs{starindex}";
  $misogffdir = "$topdir/$defs{misogffdir}" if defined $defs{misogffdir};
  $subreadindex = "$topdir/$defs{subreadindex}" if defined $defs{subreadindex};
  $bwaindex = "$topdir/$defs{bwaindex}" if defined $defs{bwaindex};
  $stardir = $defs{stardir};
  $tophatdir = $defs{tophatdir};
  $bowtiedir = $defs{bowtiedir};
  $subreaddir = $defs{subreaddir};
  $samtoolsdir = $defs{samtoolsdir};
  $picarddir = $defs{picarddir};
  $cufftoolsdir = $defs{cufftoolsdir};
  $misobin = $defs{misobin};
  $bwadir = $defs{bwadir};
  
  $CleanExclude = $defs{clean};
  $SpikeinCleanExclude = $defs{spclean};
  
#  $genome = "$topdir/$defs{genome}";
  $gtffile = "$topdir/$defs{gtffile}";
  $genomefile = "$topdir/$defs{genomefile}";
  $paired = $defs{paired};
  $stranded = $defs{stranded};
}

sub ReadExperiment
{
  local *F;
  open F, $expfile or die "Cannot open $expfile\n";
  
  @experiment = ();
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^\s*$/;
    my ($file, $cond, $rep, $pair) = split " ", $line;
    push @experiment, [$file, $cond, $rep, $pair];
  }
}



1;
