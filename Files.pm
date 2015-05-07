#!/usr/bin/perl
#
package ViennaRNA::Files;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(read_dotplot);
%EXPORT_TAGS = ();




sub read_dotplot{
  my $dotplot = shift;
  my %mfe;
  my %probs;
  my $seq;
  my $seq_start;

  if( -e $dotplot ){
    open(DOT, "<".$dotplot);
    while(<DOT>){
      if(defined($seq_start)){
        undef($seq_start), next  unless /^([AUCGTNaucgtn]+)/;
        $seq .= $1;
      }

      $seq_start = 1 if /^\/sequence/;

      next unless /(\d+)\s+(\d+)\s+([0-9.Ee-]+)\s+(ubox|lbox)/;
      my ($i, $j, $p) = ($1,$2,$3);
      if($4 eq "ubox"){
        $p *= $p;
        $probs{$i,$j} = $p;
      } else {
        $mfe{$i,$j} = 1;
      }
    }
    close(DOT);
    return ($seq, \%mfe, \%probs);
  } else {
    return undef;
  }
}

1;
