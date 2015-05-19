package RNA::Design;
use strict;
use warnings;
use RNA;
use Carp;

use Exporter;

our $VERSION    = 1.00;
our @ISA        = qw(Exporter);
our @EXPORT     = ();
our @EXPORT_OK  = ();

=head1 AUTHOR

Stefan Badelt (stef@tbi.univie.ac.at)

=head1 NAME

RNA::Design -- plug-and-play design of nucleic acid sequences

=head1 DESCRIPTION

This package provides various subroutines to design RNA molecules optimized for
single or multiple conformations. The properties of the RNA sequence have to be
specified in the Main Object: $Design = RNA::Design->new(). Currently, this package
supports every structure constraint that can be specified by two separate
well-formed dot-bracket strings. Sequence optimization functions can be composed of
the energies of input structures 'eos(s)', the ensemble-free-energy 'gfe()',
the free energy of a constrained ensemble 'pfc(s)' and conditional probabilites
of certain structures 'prob(s1,s2)'. All of these functions exist for linear
and circular sequences and they allow to specify a temperature.

Additionally, cost-functions are multiplied with a term that corrects for specified 
base-probabilities (see set_prob()) and penalties for particular subsequences that 
shall be avoided (see set_avoid()).

=head1 METHODS

=cut

=head2 new()

Initialize the Global Object for Nucleic Acid Design. It contains various
public elements, such as a list of structures specified in the cost-
function, the cost-function itself, a probability distribution of bases, 
and a set of penalized nucleic acid sequences. See get/set routines for 
description of public functions.

=cut


sub new {
  my $class = shift;
  my $self = {
    #cutp  => -1,
    fibo  => [0,1],
    border=> 0,
    nos   => undef,
    clist => [],
    prob  => ({A => 0.25, C => 0.25, G => 0.25, U => 0.25}),
    avoid => ['AAAAA','CCCCC','GGGGG','UUUUU'],
    verb  => 0,

    structures  => [],
    constraint  => '',
    optfunc     => 'eos(1)+eos(2) - 2*gfe() + 0.3*(eos(1)-eos(2)+0.00)**2',

    iupack => ({
      'A' => 'A',
      'C' => 'C',
      'G' => 'G',
      'T' => 'U',
      'U' => 'U',
      'R' => 'AG',  # purine
      'Y' => 'CU',  # pyrimidine
      'S' => 'CG',
      'M' => 'AC',
      'W' => 'AU',
      'K' => 'GU',
      'V' => 'ACG', # not T
      'H' => 'ACU', # not G
      'D' => 'AGU', # not C
      'B' => 'CGU', # not A
      'N' => 'ACGU'
    }),
    neighbor => ({  # ACGU => ACGU
      'A' => 'U',   # 1000 => 0001 
      'C' => 'G',   # 0100 => 0010 
      'G' => 'Y',   # 0010 => 0101 
      'U' => 'R',   # 0001 => 1010 
      'R' => 'Y',   # 1010 => 0101 
      'Y' => 'R',   # 0101 => 1010 
      'S' => 'B',   # 0110 => 0111 
      'M' => 'K',   # 1100 => 0011 
      'W' => 'D',   # 1001 => 1011 
      'K' => 'N',   # 0011 => 1111 
      'V' => 'B',   # 1110 => 0111 
      'H' => 'D',   # 1101 => 1011 
      'D' => 'N',   # 1011 => 1111 
      'B' => 'N',   # 0111 => 1111 
      'N' => 'N',   # 1111 => 1111 
    }),
    iupack_bin => ({# ACGU
      'A' => 8,     # 1000,
      'C' => 4,     # 0100,
      'G' => 2,     # 0010,
      'T' => 1,     # 0001,
      'U' => 1,     # 0001,
      'R' => 10,    # 1010,  # purine
      'Y' => 5,     # 0101,  # pyrimidine
      'S' => 6,     # 0110, 
      'M' => 12,    # 1100, 
      'W' => 9,     # 1001, 
      'K' => 3,     # 0011, 
      'V' => 14,    # 1110,  # not T
      'H' => 13,    # 1101,  # not G
      'D' => 11,    # 1011,  # not C
      'B' => 7,     # 0111,  # not A
      'N' => 15,    # 1111,
    }),
    bin_iupack => [ # ACGU
      '',           # 0000  0 
      'U',          # 0001  1 
      'G',          # 0010  2
      'K',          # 0011  3
      'C',          # 0100  4
      'Y',          # 0101  5
      'S',          # 0110  6
      'B',          # 0111  7
      'A',          # 1000  8
      'W',          # 1001  9 
      'R',          # 1010 10
      'D',          # 1011 11
      'M',          # 1100 12
      'H',          # 1101 13
      'V',          # 1110 14
      'N'           # 1111 15
    ],
  };

  bless $self, $class;
  return $self;
}


=head2 get/set parameters

Get and Set stuff

=cut


{ # get/set routines
  sub set_verb {
    my ($self, $var) = @_;
    $self->{verb}=$var;
    return $self->{verb};
  }

  sub get_verb {
    my $self = shift;
    return $self->{verb};
  }

  sub set_optfunc {
    my ($self, $var) = @_;
    $self->{optfunc} = $var;
    return $self->{optfunc};
  }

  sub get_optfunc {
    my $self = shift;
    return $self->{optfunc};
  }

  sub set_constraint {
    my ($self, $var) = @_;
    croak "overwriting old constraint" if $self->{constraint};
    $self->{constraint} = $var;
    return $self->{constraint};
  }

  sub get_constraint {
    my $self = shift;
    return $self->{constraint};
  }

  sub set_avoid {
    my $self = shift;
    $self->{avoid} = [@_];
    return $self->{avoid};
  }

  sub get_avoid {
    my $self = shift;
    return $self->{avoid};
  }

  sub add_structures {
    my $self = shift;
    push @{$self->{structures}}, @_;
    return $self->{structures};
  }

  sub get_structures {
    my $self = shift;
    return $self->{structures};
  }

  sub fill_fibo {
    my ($self, $var) = @_;
    while ($#{$self->{fibo}} < $var) {
      push @{$self->{fibo}}, $self->{fibo}->[-1] + $self->{fibo}->[-2];
    }
    return $self->{fibo};
  }

  sub get_fibo {
    my ($self, $var) = @_;
    $self->fill_fibo($var) unless defined $self->{fibo}[$var];
    return $self->{fibo}[$var];
  }
}

=head2 make_dependency_graph()

TODO: make_dependency_graph(@structures). If @structures is empty, it will make
the dependency graph from all structures added with the add_strucutres()
routine. In some cases the user may not want all the structures in the cost
function to be part of the dependency graph, therfore @structures can be 
specified to select for those that shall define the dependencies.

TODO: rewrite and use RNA::Utils::make_pair_table();

=cut

sub make_dependency_graph {
  my $self = shift;
  my @seen;

  my @fist = make_pair_table($self->{structures}[0]);
  my @sest = make_pair_table($self->{structures}[1]);

  for my $i (0..$#fist) {
    my @cycle1 = ();
    next if ($seen[$i]);
    push @cycle1, $i;
    $seen[$i]=1;
    my $j = $fist[$i];
    while ($j>=0 && !$seen[$j]) {
      push @cycle1, $j;
      $seen[$j] = 1;
      $j = ($#cycle1 % 2) ? $sest[$j] : $fist[$j];
    }
    $j = $sest[$i];
    my @cycle2 = ();
    while ($j>=0 && !$seen[$j]) {
      unshift @cycle2, $j;
      $seen[$j] = 1;
      $j = ($#cycle2 % 2) ? $sest[$j] : $fist[$j];
    }
    # duplicate first element if closed cycle
    unshift @cycle2, $j if ($j!=-1);
    push @{$self->{clist}}, [@cycle2, @cycle1];
  }

  ## foreach position store its path in @plist
  #foreach my $pp (@{$self->{clist}}) {
  #  my @path = @$pp;
  #  for my $pos (0..$#path) {
  #    $self->{plist}->[$path[$pos]] = [$pp, $pos];
  #  }
  #}

  return $self->{clist};
}

=head2 explore_sequence_space()

This function needs to be called to initialize the subsequent sequence design. It reads
the previously computed dependency graph to (i) update the sequence constraint, (ii)
caclulate the total number of sequences able to fulfill sequence and structure 
constraints, (iii) set a stop-condition for sequence design dependent on the number
of possible mutation-moves.

TODO: The routine does not include sequence constraints for step (ii) and (iii) and 
therefore overestimates the sequence space. Fixing this should greatly reduce the runtime.

=cut

sub explore_sequence_space {
  my $self = shift;
  my $verb = 0;

  my @clist=@{$self->{clist}};
  my $con  = $self->{constraint};

  my ($border, $max_len, $nos) = (0,0,1);
  foreach my $path (@clist) {
    my @pseq=();
    my $l = @$path;
    my $cycle=0;

    # Identify cycles
    do{--$l; $cycle=1} if ($path->[0] eq $path->[-1] && (@$path>1));

    # Statistics
    $max_len = $l if $l > $max_len;
    $border += ($cycle) ? 
    2*($self->get_fibo($l+1)+$self->get_fibo($l-1)) : 
    2*($self->get_fibo($l+1)+$self->get_fibo($l));
    $nos    *= ($cycle) ? 
    2*($self->get_fibo($l+1)+$self->get_fibo($l-1)) : 
    2*($self->get_fibo($l+1)+$self->get_fibo($l));

    # Translate path into nucleotide constraints
    push @pseq, substr($con, $_, 1) foreach @$path;

    print "@$path => @pseq\n" if $l>1 && 0;

    # Update Constraint
    my $jailbreak=0;
    while (1) {
      my $mutated = 0;
      my $i;
      if ($cycle) {
        for $i (-$#pseq .. $#pseq-1) {
          next if $pseq[$i] eq 'N';
          $mutated = 1 if $self->rewrite_neighbor($pseq[$i], \$pseq[$i+1]);
        }
        for $i (reverse(-$#pseq .. $#pseq)) {
          next if $pseq[$i] eq 'N';
          $mutated = 1 if $self->rewrite_neighbor($pseq[$i], \$pseq[$i-1]);
        }
      } else {
        for $i (0 .. $#pseq-1) {
          next if $pseq[$i] eq 'N';
          $mutated = 1 if $self->rewrite_neighbor($pseq[$i], \$pseq[$i+1]);
        }
        for $i (reverse(1 .. $#pseq)) {
          next if $pseq[$i] eq 'N';
          $mutated = 1 if $self->rewrite_neighbor($pseq[$i], \$pseq[$i-1]);
        }
      }
      last unless $mutated;
      die "escaping constraint updates" if $jailbreak++ > 100;
    }

    substr($con, $path->[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  }
  $self->{constraint}=$con;
  $self->{border}=$border unless $self->{border};
  $self->{nos}=$nos;
  printf STDERR "Total number of sequences w.o. constraints: %g; Border: %g\n",
  $self->{nos}, $self->{border} if $verb;
  #TODO: return (constraint, border, nos);
}

sub rewrite_neighbor {
  my ($self, $c1, $c2) = @_;
  print "$c1 -> $$c2\n" if 0;

  my $nb = $self->{neighbor}->{$c1};
  if ($nb eq $$c2) {
    return 0;
  } else {
    my $bin_c2 = $self->{iupack_bin}{$$c2};
    my $bin_nb = $self->{iupack_bin}{$nb};

    my $new_nb = $self->{bin_iupack}[($bin_c2 & $bin_nb)];
    croak "sequence constraints cannot be fulfilled\n" unless $new_nb;

    return 0 if $new_nb eq $$c2;

    $$c2 = $new_nb;
    return 1;
  }
}

=head2 find_a_sequence()

Uses the depencency graph and the sequence constraint to make a random sequence. 

TODO: The random-sequence should already include the information of get_probs().

=cut


sub find_a_sequence {
  # make random sequence or randomly mutate sequence
  # by replacing all @clists by random paths
  my $self = shift;

  my @clist=@{$self->{clist}};
  my $con  = $self->{constraint};
  my $seq  = $con;

  foreach my $path (@clist) {
    my @pseq = ();
    my $l = @$path;

    # Identify cycles
    do{--$l; $l=-$l} if ($$path[0] eq $$path[-1] && (@$path>1));

    # Translate path into nucleotide constraints
    push @pseq, substr($con, $_, 1) foreach @$path;

    @pseq = $self->make_pathseq($l, @pseq);

    substr($seq, $$path[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  }
  return $seq;
}

=head2 optimize_sequence(sequence, maximum_number_of_mutations)

The standard optimization function. Whenever a sequence mutation results in a
better score, it replaces the current solution. In case there are too many 
useless mutations, (see explore_sequence_space()), the current sequence is returned.

=cut

sub optimize_sequence {
  my $self  = shift;
  my $refseq= shift;
  my $m     = shift;

  my $verb    = $self->{verb};
  my $border  = $self->{border};
  my $refcost = $self->eval_sequence($refseq);
  my ($mutseq,$newcost);

  my $reject = 0;
  my %seen = ();
  for my $d (1..$m) {
    $mutseq = $self->mutate_seq($refseq);
    if (!exists $seen{$mutseq}) {
      $seen{$mutseq}=1;
      $newcost = $self->eval_sequence($mutseq);
      if ($newcost < $refcost) {
        $refseq = $mutseq;
        $refcost= $newcost;
        $reject = 0;
        printf STDERR "%4d %s %6.3f\n", $d, $mutseq, $newcost if $verb;
      } 
    }
    else {
      $seen{$mutseq}++;
      last if (++$reject >= $border);
    }
  }
  return $refseq;
}

=head2 mutate_seq()

Choose a random cycle, mutate it using make_pathseq of the sequence constraint.

TODO: that could result in the same sequence as before, which is inefficient!

=cut

sub mutate_seq {
  my $self  = shift;
  my $seq   = shift;

  my $con   = $self->{constraint};
  my @clist = @{$self->{clist}};

  my ($path, @pseq);

  $path = $clist[int rand @clist];
  my $l = @$path;
  do{--$l; $l=-$l} if ($$path[0] eq $$path[-1] && (@$path>1));
  push @pseq, substr($con, $_, 1) foreach @$path;

  @pseq = $self->make_pathseq($l, @pseq);

  substr($seq, $$path[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  return $seq;
}

=head2 make_pathseq()

Takes an Array of IUPACK code and randomly rewrites it into a valid array of Nucleotides

=cut

sub make_pathseq {
  my $self = shift;
  my $len  = shift;
  my @pseq = @_; # ('N', 'N', 'R')
  my @path;

  my %iupack   = %{$self->{iupack}};

  my $cycle=0;
  if ($len < 0) {
    $len = -$len;
    $cycle=1;
  }

  for (my $i=0; $i<=$#pseq; ++$i) {
    my $c = $pseq[$i];
    my @i = (split '', $iupack{$c});
    $path[$i] = $i[int rand @i];

    if ($i==0 && $cycle) {
      my ($a, $j) = ($path[$i], -1);
      do { $a=$pseq[$j--] } while $self->rewrite_neighbor($a, \$pseq[$j]);
    }

    $self->rewrite_neighbor($path[$i], \$pseq[$i+1]) if $i < $#pseq;
  }

  return @path;
}

=head2 eval_sequence()

Evaluate a given sequence according to your cost function.

=cut

sub eval_sequence {
  my $self = shift;
  my $seq  = shift;

  my $verb = $self->{verb};
  my $ofun = $self->{optfunc};
  print STDERR $ofun."\n" if $verb > 1;

  my %results;
  $RNA::fold_constrained=1;
  foreach my $func (qw/eos eos_circ pfc pfc_circ gfe gfe_circ prob prob_circ/) {
    while ($ofun =~ m/$func\(([0-9\,\s]*)\)/) {
      print STDERR "next: $&\n" if $verb > 1;
      $results{$&} = eval "\$self->$func(\$seq, $1)" unless exists $results{$&};
      $ofun =~ s/$func\(([0-9\,\s]*)\)/ $results{$&} /;
    }
  }
  $RNA::fold_constrained=0;

  print STDERR $ofun."\n" if $verb > 1;
  my $r = eval $ofun;
  warn $@ if $@;

  my $p = $self->base_prob($seq);
  my $a = $self->base_avoid($seq, 5);

  return $r*$p*$a;
}

sub base_avoid {
  my ($self, $seq, $pen) = @_;
  my $penalty=1;

  my @avoid = @{$self->{avoid}};
  foreach my $string (@avoid) {
    $penalty *= $pen while ($seq =~ m/$string/g);
  }
  return $penalty;
}

sub base_prob {
  my ($self, $seq) = @_;

  my %prob = %{$self->{prob}};
  my %base;
  $base{$_}++ foreach (split //, $seq);

  my $cost=1;
  foreach my $k (keys %base) {
    $base{$k} /= length $seq;
    #print "Dist $k $base{$k} vs $prob{$k} => \n";
    if ($base{$k} > $prob{$k}) {
      $cost *= $base{$k}/$prob{$k};
    } else {
      $cost *= $prob{$k}/$base{$k};
    }
  }
  return $cost;
}

sub prob {
  my ($self, $seq, $i, $j, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;

  my $kT=0.6163207755;
  my $s_i = ($i) ? $self->{structures}[$i-1] : undef;
  my $s_j = ($j) ? $self->{structures}[$j-1] : undef;
  my $tmp;

  $tmp=$s_i; my $dGi = RNA::pf_fold($seq, $tmp);
  $tmp=$s_j; my $dGj = RNA::pf_fold($seq, $tmp);

  return exp(($dGj-$dGi)/$kT);
}

sub prob_circ {
  my ($self, $seq, $i, $j, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;

  my $kT=0.6163207755;
  my $s_i = ($i) ? $self->{structures}[$i-1] : undef;
  my $s_j = ($j) ? $self->{structures}[$j-1] : undef;
  my $tmp;

  $tmp=$s_i; my $dGi = RNA::pf_circ_fold($seq, $tmp);
  $tmp=$s_j; my $dGj = RNA::pf_circ_fold($seq, $tmp);

  return exp(($dGj-$dGi)/$kT);
}

sub eos {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  croak "Could not find ".$i."th structure\n" unless $str;
  return RNA::energy_of_struct($seq, $str);
}

sub eos_circ {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  croak "Could not find ".$i."th structure\n" unless $str;
  return RNA::energy_of_circ_struct($seq, $str);
}

sub pfc {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  # seemingly useless string modification
  # to avoid ViennaRNA SWIG interface bugs
  $str.='.'; $str=substr($str,0,-1);
  croak "Could not find ".$i."th structure\n" unless $str;
  return RNA::pf_fold($seq, $str);
}

sub pfc_circ {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  croak "Could not find ".$i."th structure\n" unless $str;
  # seemingly useless string modification
  # to avoid ViennaRNA SWIG interface bugs
  $str.='.'; $str=substr($str,0,-1);
  return RNA::pf_circ_fold($seq, $str);
}

sub gfe {
  my ($self, $seq, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  return RNA::pf_fold($seq, undef);
}

sub gfe_circ {
  my ($self, $seq, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  return RNA::pf_circ_fold($seq, undef);
}

sub make_pair_table {
  my $struc = shift;
  my ($i, $j, @stack, @ptable);

  $i=0;
  foreach my $c (split(//,$struc)) {
    if ($c eq '(') {
      push @stack, $i;
    } elsif ($c eq ')') {
      $j = pop @stack;
      croak "ERROR [".((caller(0))[3])."] unbalanced brackets!" unless defined $j;
      $ptable[$i]=$j;
      $ptable[$j]=$i;
    } else {
      $ptable[$i]=-1;
    }
    $i++;
  }
  croak "ERROR [".((caller(0))[3])."] unbalanced brackets!" if @stack != 0;
  return @ptable;
}

1;
