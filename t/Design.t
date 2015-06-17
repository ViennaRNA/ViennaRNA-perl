#!/usr/bin/perl

use strict;
use warnings;
use Test::More tests => 11;

# ****************************************** #
# ~~~~~~~~~~ ALL INTERNAL MODULES ~~~~~~~~~~ #
# .......................................... #

BEGIN { use_ok('RNA'); }
BEGIN { use_ok('RNA::Design'); }
BEGIN { use_ok('RNA::Utils'); }

my $ViennaDesign = RNA::Design->new(verb=>1);
isa_ok($ViennaDesign, 'RNA::Design', 'Object initialized');

# Check default Options
subtest 'RNA::Design -- Default Options' => sub {
  plan tests => 3;
  is ($ViennaDesign->get_verb, 0, 'get_verb()');
  is ($ViennaDesign->get_constraint, '', 'get_constraint()');
  is (ref $ViennaDesign->get_structures, "ARRAY", 'get_structures()');
};

subtest 'RNA::Design -- Internal Functions' => sub {
  plan tests => 13;
  my ($nuc1, $nuc2) = ('R', 'N');
  is ($ViennaDesign->rewrite_neighbor($nuc1, \$nuc2), 1, 'rewrite_neighbor() - change');
  is ($ViennaDesign->rewrite_neighbor($nuc1, \$nuc2), 0, 'rewrite_neighbor() - constant');

  my @path = ('N','N','N','N','C','N','N','N');
  my @outp = ('Y','R','Y','G','C','G','Y','R');
  my $cycle;
  is_deeply ([$ViennaDesign->update_constraint($cycle = 1, @path)], \@outp, 'update_constraint() - cycle');
  is_deeply ([$ViennaDesign->update_constraint($cycle = 0, @path)], \@outp, 'update_constraint() - path');

  # TODO: test these:
  is ($ViennaDesign->enumerate_pathways($cycle = 0, @path),
    $ViennaDesign->enumerate_pathways($cycle = 0, @outp), 'enumerate_pathways() - no overcounting');

  is ($ViennaDesign->enumerate_pathways($cycle = 0, ('A','N','N','N')), 3, 'enum_pathways() - path');
  is ($ViennaDesign->enumerate_pathways($cycle = 1, ('A','N','N','N')), 2, 'enum_pathways() - cycle');

  @path =     ('N','N','N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','R','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','B','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','C','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
  #@path = ('N','N','N','N','N','N','N','N','N');
  
  #my @b;# = ([0,0,0,0]) x @path;
  #push @b, [0,0,0,0] foreach @path;
  #my %bb = (
  #  'A' => 0,
  #  'C' => 1,
  #  'G' => 2,
  #  'U' => 3,
  #);

  #@path = $ViennaDesign->update_constraint($cycle = 0, @path);
  #print "@path\n";
  #for my $r (0..1000) {
  #  my @tmp = $ViennaDesign->make_pathseq($cycle = 0, @path);
  #  #print "@tmp\n";
  #  for my $t (0 .. $#tmp) {
  #    my $idx = $bb{$tmp[$t]};
  #    $b[$t]->[$idx]+=1;
  #  }
  #}
  #print "@$_\n" foreach @b;

  @path = @outp;

  @outp = ('C','G','C','G','C','G','U','A');
  srand(30); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  @outp = ('C','G','C','G','C','G','U','G');
  srand(30); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

  @outp = ('C','G','U','G','C','G','C','G');
  srand(1); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  srand(1); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

  @outp = ('U','A','U','G','C','G','U','G');
  srand(2); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  srand(2); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

};


my $plist;
my @structs;

@structs = ( '......' );
$plist = [[0],[1],[2],[3],[4],[5]];
is_deeply($ViennaDesign->find_dependency_paths(@structs), $plist, 'find_dependency_paths -- open chain');

@structs = ( '.(...)','.(...)');
$plist = [[0],[1,5],[2],[3],[4]];
is_deeply($ViennaDesign->find_dependency_paths(@structs), $plist, 'find_dependency_paths -- one bpair');

@structs = (
  #0....,....1....,....2....,....3',
  '.(....)(((...)))',
  '.(...).((....)).',
  '(..............)',
);

$plist = [[9,13,8,14,7,15,0],[5,1,6],[2],[3],[4],[10],[11],[12]];
is_deeply($ViennaDesign->find_dependency_paths(@structs), $plist, 'find_dependency_paths -- two pathways');

@structs = (
  #0....,....1....,....2....,....3',
  '.(...)((...)).(...)',
  '.(...((.(...)))...)',
  #NNNNNNNNNNNNNNNNNNN',
);

$plist = [[0],[18,1,5,14,18],[2],[3],[4],[13,6,12,8],[7,11],[9],[10],[15],[16],[17]];
is_deeply($ViennaDesign->find_dependency_paths(@structs), $plist, 'find_dependency_paths -- cycles,pathways,basepairs');

$ViennaDesign->set_constraint('GNNNNNNNNNNNNANNNUA');
my @explore = ('GUNNNRUNYNNNRAUNNUA', 39, 589824);
is_deeply([$ViennaDesign->explore_sequence_space()], \@explore, 'explore_sequence_space()');


print "done testing\n";

# TODO: 
# $sequence = find_a_sequence(@clist, $constraint);
# $optseq = optimize_sequence($refseq, $m);
# $cost = eval_sequence($seq, $optfunc);
#
# mutate_seq()
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

