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
  is ($ViennaDesign->get_verbosity, 0, 'get_verbosity()');
  is ($ViennaDesign->get_constraint, '', 'get_constraint()');
  is (ref $ViennaDesign->get_structures, "ARRAY", 'get_structures()');
};

subtest 'RNA::Design -- Internal Functions' => sub {
  plan tests => 16;
  my ($nuc1, $nuc2) = ('R', 'N');
  is ($ViennaDesign->rewrite_neighbor($nuc1, \$nuc2), 1, 'rewrite_neighbor() - change');
  is ($ViennaDesign->rewrite_neighbor($nuc1, \$nuc2), 0, 'rewrite_neighbor() - constant');

  my @path = ('N','N','N','N','C','N','N','N');
  my @outp = ('Y','R','Y','G','C','G','Y','R');
  my $cycle;
  is_deeply ([$ViennaDesign->update_constraint($cycle = 1, @path)], \@outp, 'update_constraint() - cycle');
  is_deeply ([$ViennaDesign->update_constraint($cycle = 0, @path)], \@outp, 'update_constraint() - path');

  is ($ViennaDesign->enumerate_pathways($cycle = 0, @outp), 15, 'enumerate_pathways() - path_long');
  is ($ViennaDesign->enumerate_pathways($cycle = 1, @outp), 13, 'enumerate_pathways() - cycle_long');

  @path = ('A','N','N','N');
  @outp = ('A','U','R','Y');
  is_deeply ([$ViennaDesign->update_constraint($cycle = 0, @path)], \@outp, 'update_constraint() - path');
  is ($ViennaDesign->enumerate_pathways($cycle = 0, @outp), 3, 'enum_pathways() - path');
  @outp = ('A','U','R','U');
  is_deeply ([$ViennaDesign->update_constraint($cycle = 1, @path)], \@outp, 'update_constraint() - cycle');
  is ($ViennaDesign->enumerate_pathways($cycle = 1, @outp), 2, 'enum_pathways() - cycle');

  @path = ('Y','R','Y','G','C','G','Y','R');
  @outp = ('C','G','C','G','C','G','U','A');
  srand(1); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  @outp = ('C','G','C','G','C','G','U','G');
  srand(6); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

  @outp = ('C','G','U','G','C','G','C','G');
  srand(20); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  srand(15); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

  @outp = ('U','A','U','G','C','G','U','G');
  srand(18); is_deeply ([$ViennaDesign->make_pathseq($cycle = 0, @path)], \@outp, 'make_pathseq() - path');
  srand(11); is_deeply ([$ViennaDesign->make_pathseq($cycle = 1, @path)], \@outp, 'make_pathseq() - cycle');

  # Playground to test pathway-math
  @path =     ('N','N','N','N','G','N');
  push @path, ('N','N','R','N','N','N');
  push @path, ('N','N','N','N','N','N');
  #push @path, ('N','N','N','N','N','N');
 
  @path = $ViennaDesign->update_constraint($cycle = 0, @path);
  my $count = $ViennaDesign->enumerate_pathways($cycle = 0, @path);
  # print "@path $count\n";
  
  # check if pathways are drawn with even probability
  # my %resp;
  # for my $r (0..10000) {
  #   my $tmp = join '', $ViennaDesign->make_pathseq($cycle = 0, @path);
  #   $resp{$tmp}++;
  # }
  # print "$_: $resp{$_}\n" foreach sort keys %resp;

  # for my $r (0..20) {
  #   srand ($r);
  #   my @tmp = $ViennaDesign->make_pathseq($cycle = 0, @path);
  #   print "$r: @tmp\n";
  # }
  
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

print "\n*Start with full example!*\n";

$ViennaDesign->set_verbosity(1);
my @tstructs = (
  '.(...)((...)).(...)',
  '.(...((.(...)))...)',
#  '(...(........)...).',
#  '(...)........(...).',
);
$ViennaDesign->add_structures(@tstructs);
$ViennaDesign->set_constraint('NNNNNNNNNNNNNNNNNNN');
$ViennaDesign->find_dependency_paths();
$ViennaDesign->explore_sequence_space();
my $seq = $ViennaDesign->find_a_sequence;

my %resp;
for my $r (0..10000) {
  my $tseq = $ViennaDesign->mutate_seq($seq);
  $resp{$tseq}++;
}
print "$seq ".scalar(keys %resp)."\n";
#print "$_: $resp{$_}\n" foreach sort {$resp{$b} <=> $resp{$a}} keys %resp;
#print "$_: $resp{$_}\n" foreach sort keys %resp;

$seq = $ViennaDesign->optimize_sequence($seq, 999999);
my $cost= $ViennaDesign->eval_sequence($seq);
print "$seq "." $cost\n";
 
# 
# # TODO: 
# # $sequence = find_a_sequence(@clist, $constraint);
# # $optseq = optimize_sequence($refseq, $m);
# # $cost = eval_sequence($seq, $optfunc);
# #
# # mutate_seq()
# # base_avoid()
# # base_prob()
# # prob()
# # prob_circ()
# # eos()
# # eos_circ()
# # pfc()
# # pfc_circ()
# # gfe()
# # gfe_circ()
# 
