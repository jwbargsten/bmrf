#!/usr/bin/env perl

use warnings;
use strict;
use 5.010;
use List::MoreUtils qw/uniq/;
#use Data::Printer;
use Bio::Gonzales::Util::Cerial;

use Bio::Gonzales::Matrix::IO qw(mslurp mspew);

my $int = mslurp("topless_2lvl_int.tsv");

my %ids = map { $_->[0] => 1, $_->[1] => 1 } @$int;
my %ev = (
  IDA => 1,
  IGI => 1,
  IMP => 1,
);

my $go_parents = jslurp("/home/bargs001/projects/fun_pred/analysis/go_preprocess/2012-05-02/go_parents.json");

my $go = mslurp("ATH_GO_GOSLIM.txt.gz");

#[
#[0]  "AT1G01010",
#[1]  "gene:2200934",
#[2]  "AT1G01010.1",
#[3]  "located in",
#[4]  "nucleus",
#[5]  "GO:0005634",
#[6]  537,
#[7]  "C",
#[8]  "nucleus",
#[9]  "ISM",
#[10] "predicted protein features",
#[11] "",
#[12] "AnalysisReference:501750651",
#[13] "rkaundal",
#[14] "2012-08-31"
#]

my @go = grep { exists( $ids{ $_->[0] } ) && $_->[7] eq 'P' && exists( $ev{ $_->[9] } ) } @$go;
my %subset;
for my $g (@go) {
  my $got = $g->[5];
  my $id  = $g->[0];
  $subset{ $id . "__" . $got } = [ $id, $got, 'P' ];
  for my $gp ( @{ $go_parents->{$got} } ) {
    $subset{ $id . "__" . $gp } = [ $id, $gp, 'P' ];
  }
}

mspew( "go_labels.up.tsv", [ values %subset ] );
