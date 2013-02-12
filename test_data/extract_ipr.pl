#!/usr/bin/env perl

use warnings;
use strict;
use 5.010;
use List::MoreUtils qw/uniq/;
#use Data::Printer;
use Gonzales::Util::Cerial;

use Bio::Gonzales::Matrix::IO qw(mslurp mspew);

my $int = mslurp("topless_2lvl_int.tsv");

my %ids =  map { $_->[0] => {}, $_->[1] => {} } @$int;

my $ipr = mslurp("TAIR10_all.domains");

#p $ipr->[0];
#die;
#[
    #[0]  "AT1G08520.1",
    #[1]  "",
    #[2]  760,
    #[3]  "ProfileScan",
    #[4]  "PS50234",
    #[5]  "VWFA",
    #[6]  558,
    #[7]  754,
    #[8]  0.0,
    #[9]  "28-Oct-2010",
    #[10] "IPR002035",
    #[11] "von Willebrand factor, type "
#]


for my $i (@$ipr) {
  (my $id = $i->[0]) =~ s/\.\d+$//;
  next unless(exists($ids{$id}));
  my $dom = $i->[10];

  next if($dom eq 'NO_ASSIGNMENT');
  next if($dom eq 'NULL');

  $ids{$id}{$dom} = 1;
}

my @res;

while(my ($id, $dom) = each %ids) {
  for my $d (keys %$dom) {
    push @res, [$id, $d];
  }
}

mspew("ipr_labels.tsv", \@res);

