#!/usr/bin/perl

use strict;
use warnings;
use Test::More tests => 5;

# ****************************************** #
# ~~~~~~~~~~ ALL INTERNAL MODULES ~~~~~~~~~~ #
# .......................................... #

BEGIN { use_ok('RNA'); }
BEGIN { use_ok('RNA::Design'); }
BEGIN { use_ok('RNA::Utils'); }

my $ViennaDesign = RNA::Design->new();
isa_ok($ViennaDesign, 'RNA::Design', 'object initialized');

# Check default Options
subtest 'Done with Default Options' => sub {
  plan tests => 3;
  is ($ViennaDesign->get_verb, 0, 'get_verb()');
  is ($ViennaDesign->get_constraint, '', 'get_constraint()');
  is (ref $ViennaDesign->get_structures, "ARRAY", 'get_structures()');
};


# ******************************************************** #
# ~~~~~~~~~~ Chemistry::SundialsWrapper::Solver ~~~~~~~~~~ #
# ........................................................ #

#subtest 'RNA::Design (Private Functions)' => sub {
#  plan tests => 4;
#   TRUE  == rewrite_neighbor('R', \'N');
#   FALSE == rewrite_neighbor('R', \'Y');
#}

# TODO: 
# @clist = make_dependency_graph(@s);
# ($constraint, $border, $nos) = explore_sequence_space(@clist, $constraint);
# $sequence = find_a_sequence(@clist, $constraint);
# $optseq = optimize_sequence($refseq, $m);
# $cost = eval_sequence($seq, $optfunc);
#
# mutate_seq()
# pake_pathseq()
# 
# base_avoid()
# base_prob()
# prob()
# prob_circ()
# eos()
# eos_circ()
# pfc()
# pfc_circ()
# gfe()
# gfe_circ()
#



