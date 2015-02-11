#!/usr/bin/perl
#
package ViennaRNA::Reliability;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(rel_posent rel_probs rel_access);
%EXPORT_TAGS = ();



sub rel_posent{
  my @pp;
  my @sp;
  my $log2  = log(2);
  my $probs = shift;
  my $n     = shift;

  for my $i (1..$n){
    for my $j (($i+1)..$n){
      next unless exists $probs->{$i,$j};
      my $p  = $probs->{$i,$j};

      my $ss = ($p>0)?$p*log($p):0;
      $sp[$i] += $ss;
      $sp[$j] += $ss;
      $pp[$i] += $p;
      $pp[$j] += $p;
    }
  }

  for my $i (1..$n) {
    no warnings;  # $p[$i] may be undef
    $sp[$i] += ($pp[$i] < 1) ? (1 - $pp[$i]) * log(1 - $pp[$i]) : 0;
    $sp[$i] /= -$log2;
  }

  shift @sp; # get rid of [0] entry

  return \@sp;
}

sub rel_probs{
  my $probs = shift;
  my $mfe   = shift;
  my $n     = shift;
  my @pp;
  my @sp;
  for my $i (1..$n){
    for my $j (($i+1)..$n){
      next unless exists $probs->{$i,$j};
      my $p  = $probs->{$i,$j};

      $sp[$i] = $sp[$j] = $p if exists  $mfe->{$i,$j};

      $pp[$i] += $p;
      $pp[$j] += $p;
    }

    no warnings;  # $p[$i] may be undef
    $sp[$i] = 1-$pp[$i] if !defined $sp[$i];
  }

  shift @sp; # get rid of [0] entry

  return \@sp;
}

sub rel_access {
  my $probs = shift;
  my $mfe   = shift;
  my $n     = shift;
  my @pp;
  my @sp;
  for my $i (1..$n){
    for my $j (($i+1)..$n){
      next unless exists $probs->{$i,$j};
      my $p  = $probs->{$i,$j};

      $sp[$i] = $sp[$j] = 1-$p if exists $mfe->{$i,$j};

      $pp[$i] += $p;
      $pp[$j] += $p;
    }

    no warnings;  # $p[$i] may be undef
    $sp[$i] = 1-$pp[$i] if !defined $sp[$i];
  }

  shift @sp; # get rid of [0] entry

  return \@sp;
}


1;
