package ViennaRNA::Benchmark;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(compareStructures
                  TruePositiveRate
                  PositivePredictiveValue
                  MathewsCorrelationCoefficient
                  FMeasure);

%EXPORT_TAGS = (  ALL => [qw(&compareStructures &TruePositiveRate &PositivePredictiveValue &MathewsCorrelationCoefficient &FMeasure)],
                  MEASURES => [qw(&TruePositiveRate &PositivePredictiveValue &MathewsCorrelationCoefficient &FMeasure)],
                  MATCHERS => [qw(&compareStructures)]);

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

The optional parameter $verbose is mainly for debugging purposes, but may be
used to produce progress output of the comparison.

The subroutine returns a set of four values ($TP, $FP, $TN, $FN), representing
the true positives $TP, the true negatives $TN, the false negatives $FN, and
the false positives $FP.

=cut

sub compareStructures{
  my $pt_gold   = shift;
  my $pt_other  = shift;
  my $fuzzy     = shift;
  my $verbose   = shift;

  $fuzzy    = 0 if !defined($fuzzy);
  $verbose  = 0 if !defined($verbose);

  my ($TP, $TN, $FN, $FP) = (0, 0, 0, 0);

  # number of false positive, but compatible pairs
  my $compatible = 0;

  # count the number of base pairs in $pt_gold, and $pt_other
  my $bps_gold  = 0;
  my $bps_other = 0;

  for my $i (1 .. $pt_gold->[0] ){
    next if $pt_gold->[$i] < $i;
    $bps_gold++;
  }
  for my $i (1 .. $pt_other->[0] ){
    next if $pt_other->[$i] < $i;
    $bps_other++;
  }

  # got through $pt_other and compare its pairs to $pt_gold
  for my $i ( 1 .. $pt_other->[0] ) {
    next if ( $pt_other->[$i] < $i );

    my $j                = $pt_other->[$i];
    my $is_true_positive = 0;
    my $is_inconsistent  = 0;
    my $is_contradicting = 0;

    #############################################################
    # TRUE POSITIVES
    #############################################################

    # let's see if position i matches in gold position j +/- fuzzy
    for my $add ( 0 .. $fuzzy ) {
      if ( $pt_gold->[$i] == $j + $add and $j + $add <= $pt_other->[0] ) {
        $is_true_positive = 1;
        $bps_gold--;
        print STDERR "BP ($i, $j) in Reference matches BP ($i, ", ( $j + $add ), ") in Prediction.\n" if ($verbose);
        last;
      }
      if ( $pt_gold->[$i] == $j - $add and $j - $add >= 1 ) {
        $is_true_positive = 1;
        $bps_gold--;
        print STDERR "BP ($i, $j) in Reference matches BP ($i, ", ( $j - $add ), ") in Prediction.\n" if ($verbose);
        last;
      }
    }

    # let's see if position j matches in gold position i +/- fuzzy
    if ( $is_true_positive == 0 ) {
      for my $add ( 0 .. $fuzzy ) {
        if ( $pt_gold->[$j] == $i + $add and $i + $add <= $pt_other->[0] ) {
          $is_true_positive = 1;
          $bps_gold--;
          print STDERR "BP ($i, $j) in Reference matches BP ($j, ", ( $i + $add ), ") in Prediction.\n" if ($verbose);
          last;
        }
        if ( $pt_gold->[$j] == $i - $add and $i - $add >= 1 ) {
          $is_true_positive = 1;
          $bps_gold--;
          print STDERR "BP ($i, $j) in Reference matches BP ($j, ", ( $i - $add ), ") in Prediction.\n" if ($verbose);
          last;
        }
      }
    }

    #############################################################
    # FALSE POSITIVES
    #############################################################

    # let's check if base-pair i * j is inconsistent with reference structure
    # base i or base j in reference is paired to something else
    my $tmp_string1 = '';
    if ( $is_true_positive == 0 ) {
      if ( $pt_gold->[$i] != 0 ) {
        $is_inconsistent = 1;
        $tmp_string1     = "Reference: ($i, $pt_gold->[$i]) ";
      }
      if ( $pt_gold->[$j] != 0 ) {
        $is_inconsistent = 1;
        $tmp_string1 .= "Reference: ($j, $pt_gold->[$j]) ";
      }
    }

    # let's check if a base-pair is contradicting
    # there is a pair k * l in ther reference that k < i < l < j
    my $tmp_string2 = '';
    for my $k ( 1 .. $pt_gold->[0] ) {
      next if ( $pt_gold->[$k] < $k);
      my $l = $pt_gold->[$k];
      if ( $k < $i and $i < $l and $l < $j ) {
        $tmp_string2      = "crossing pairs: $k < $i < $l < $j";
        $is_contradicting = 1;
        last;
      }
    }

    $TP++ if ( $is_true_positive == 1 );

    if ($is_true_positive == 0){
      if ( $is_inconsistent == 1 or $is_contradicting == 1 ) {
        print STDERR "BP ($i, $j) is a false positive." if ( $verbose );
        print STDERR "\t$tmp_string1" if ( $verbose and $tmp_string1 ne '');
        print STDERR "\t$tmp_string2"  if ( $verbose and $tmp_string2 ne '');
        print STDERR "\n" if ( $verbose );
      } else {
        print STDERR "BP ($i, $j) is a false positive, but compatible with Reference.\n" if ( $verbose );
        $compatible++;
      }
      $FP++;
    }
  }

  print "Reference contains ", $bps_gold, " BPs that were not matched to Prediction\n" if ( $verbose );

  $FN = $bps_gold;

  # make an educated guess for the number of true negatives
  $TN = 0.5 * ($pt_gold->[0] * ($pt_gold->[0] - 1));

  return ($TP, $FP, $TN, $FN, $compatible);
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

=over

=item * ViennaRNA::Utils  for converting secondary structures to pair_table format

=item * ViennaRNA::Files  for parsing secondary structures from different input file formats

=back

=cut
