package ViennaRNA::Benchmark;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ();

=head1 AUTHOR

Ronny Lorenz (ronny@tbi.univie.ac.at)

=head1 NAME

ViennaRNA::Benchmark - A set of subroutines to benchmark RNA secondary structure
prediction performance

=head1 DESCRIPTION

This package provides various subroutines for benchmarking the prediction
performance of RNA secondary structure prediction tools.
The basic data structure for comparing secondary structures is a pair table,
which can be generated from dot-bracket annotation via the make_pair_table()
subroutine of the ViennaRNA::Utils package.

=head1 METHODS

=cut


=head2 CompareStructures($pt_gold, $pt_other, $fuzzy = 0, $verbose = 0)

Compare a secondary structure to a 'gold standard' to assess the number of
false positives (FP) and false negatives (FN)

Takes a trusted secondary structure ($pt_gold), and a predicted secondary
structure ($pt_other) provided as references to pair_tables to compute the
number of base pairs in the trusted structure (true positives, or TP), the
number of base pairs in the predicted structure that are not in the trusted
structure (false positives, or FP), and the number of base pairs in the
trusted structure but not in the predicted one (false negatives, or FN).
The maximal number of base pairs that are not part of the trusted structure,
the true negatives (TN), are estimates as n*(n-1)/2 for an RNA sequence of
length n.

The optional parameter $fuzzy allows for base pair slippage. Hence, for any
base pair (i,j) in the gold standard set of pairs ($pt_gold), a base pair
(p, q) in the second structure ($pt_other) is considered a true positive, if
i - $fuzzy <= p <= i + $fuzzy, and j - $fuzzy <= q <= j + $fuzzy.

The optional parameter $verbose may be used to produce progress output
of the comparison.

The subroutine returns a set of four values ($TP, $TN, $FN, $FP), representing
the true positives $TP, the true negatives $TN, the false negatives $FN, and
the false positives $FP.

=cut
sub compareStructures{
  my $pt_gold   = shift;
  my $pt_other  = shift;
  my $fuzzy     = shift;
  my $verbose   = shift;

  my ($TP, $TN, $FN, $FP);


  return ($TP, $TN, $FN, $FP);
}

=head2 TruePositiveRate($TP, $FP, $FN)

Compute the true positive rate (TPR), a.k.a. sensitivity

Takes the number of true positives ($TP), false positives ($FP),
and false negatives ($FN) to return the True Positive Rate.

=cut
sub TruePositiveRate{
  my ($TP, $FP, $FN) = @_;
  return $TP / ($TP + $FN) if ($TP + $FN > 0);
  return 0;
}

=head2 PositivePredictiveValue($TP, $FP, $FN)

Compute the positive predictive value (PPV), a.k.a. selectivity

Takes the number of true positives ($TP), false positives ($FP),
and false negatives ($FN) to return the Positve Predictive Value.

=cut
sub PositivePredictiveValue{
  my ($TP, $FP, $FN) = @_;
  return $TP / ($TP + $FP) if ($TP + $FP > 0);
  return 0;
}

=head2 MathewsCorrelationCoefficient($TP, $FP, $FN, $TN)

Compute the Mathews Correlation Coefficient (MCC)

Takes the number of true positives ($TP), false positives ($FP),
false negatives ($FN), and true negatives ($TN) to return the
Positve Predictive Value.

=cut
sub MathewsCorrelationCoefficient{
  my ($TP, $FP, $FN, $TN) = @_;
  return ($TP * $TN - $FP * $FN) / sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN));
}

=head2 FMeasure($PPV, $TPR)

Compute the F1-Measure

Takes the positive predictive value ($PPV), and the true positive
rate ($TPR) to compute the F1-Measure.

=cut
sub FMeasure{
  my ($ppv, $tpr) = @_;
  return 2 * ($ppv * $tpr) / ($ppv + $tpr) if ($ppv + $tpr > 0);
  return 0;
}

1;

=head1  SEE ALSO

ViennaRNA::Utils  for converting secondary structures to pair_table format
ViennaRNA::Files  for parsing secondary structures from different input file formats

=cut
