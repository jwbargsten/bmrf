#!/usr/bin/env perl

use warnings;
use strict;
use 5.010;
use Data::Dumper;
use Bio::Gonzales::Matrix::IO qw(mslurp);
use Gonzales::Util::Cerial;



my $int = mslurp("topless_2lvl_int.tsv");
my %ids = map { $_->[0] => 0, $_->[1] => 0 } @$int;

my $go = mslurp("go_labels.up.tsv");

for my $g (@$go) {
  $ids{$g->[0]}++;
}
jspew("annotation_cnt.json", \%ids);
