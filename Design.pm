package RNA::Design;
use strict;
use warnings;
use RNA;
use RNA::Utils;
use Carp;
use Data::Dumper;

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
    plist => [],
    rlist => [],
    avoid => ['AAAAA','CCCCC','GGGGG','UUUUU'],
    avoid_penalty => 5,
    base_probs  => ({A => 0.25, C => 0.25, G => 0.25, U => 0.25}),
    verb  => 0,

    structures  => [],
    constraint  => '',
    optfunc     => 'eos(1)+eos(2) - 2*gfe() + 0.3*(eos(1)-eos(2)+0.00)**2',

    base => ({
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'U' => 3,
      }),

    pair => [
      #A C G U
      [0,0,0,1], # A 
      [0,0,1,0], # C
      [0,1,0,1], # G
      [1,0,1,0]  # U
    ],

    max_const_plen => 20,
    solution_space => ({

      }),

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
    carp "overwriting old constraint" if $self->{constraint};
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

  sub set_avoid_penalty {
    my ($self, $var) = @_;
    $self->{avoid_penalty} = $var;
    return $self->{avoid_penalty};
  }

  sub get_avoid_penalty {
    my $self = shift;
    return $self->{avoid_penalty};
  }

  sub set_base_probs {
    my ($self, $var) = @_;
    $self->{base_probs} = $var;
    return $self->{base_probs};
  }

  sub get_base_probs {
    my $self = shift;
    return $self->{base_probs};
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

=head2 find_dependency_paths(@s)

If @s is empty, it will make the dependency graph from all structures added
with the add_structures() routine. In some cases the user may not want all the
structures in the cost function to be part of the dependency graph, therfore
@structures can be specified to select for those that shall define the
dependencies. 

=cut

sub find_dependency_paths {
  my $self = shift;
  my @structures = @_;
  my (@pt1, @pt2, @lx1, @lx2, @seen);

  if ($self->{plist}) {
    #carp "overwriting old dependency-path!\n";
    $self->{plist}=[];
  }

  @structures = @{$self->{structures}} unless @structures;
  croak "no structure input found" unless @structures;

  # Merge all strucutes into two pairtables and complain if not possible!
  for (my $s=0; $s<@structures; ++$s) {
    my @pt = @{RNA::Utils::make_pair_table($structures[$s])};
    if ($s == 0) {
      @pt1 = @pt;
      @pt2 = (0) x @pt1;
    } else {
      @lx1 = RNA::Utils::make_loop_index_from_pt(@pt1);
      @lx2 = RNA::Utils::make_loop_index_from_pt(@pt2);
      for my $i (1..$#pt) {
        next unless $i < $pt[$i]; # '('
        my $j = $pt[$i];
        # the base-pairs are in there already
        next if ($pt1[$i]==$j && $pt1[$j]==$i);
        next if ($pt2[$i]==$j && $pt2[$j]==$i);
        # new base-pairs:
        if ($lx1[$i]==$lx1[$j] && !$pt1[$i] && !$pt1[$j]) {
          $pt1[$i]=$pt[$i];
          $pt1[$j]=$pt[$j];
        } elsif ($lx2[$i]==$lx2[$j] && !$pt2[$i] && !$pt2[$j]) {
          $pt2[$i]=$pt[$i];
          $pt2[$j]=$pt[$j];
        } else {
          carp "ignoring basepair $i - $j from structure ".($s+1);
        }
      }
    }
  }

  # Construct a list of depencency-pathways (ZERO-Based)
  for my $i (1..$#pt1) {
    next if $seen[$i];

    # initialize path with $i
    my @path  = ($i-1);
    $seen[$i] = 1;

    # expand @path to the right (push)
    my ($pt, $j) = (1, $pt1[$i]);
    while ($j && !$seen[$j]) {
      push @path, $j-1;
      $seen[$j]=1;
      $j  = ($pt == 1) ? $pt2[$j] : $pt1[$j];
      $pt = ($pt == 1) ? 2 : 1;
    }

    # expand @path to the left (unshift)
    ($pt, $j) = (2, $pt2[$i]);
    while ($j && !$seen[$j]) {
      unshift @path, $j-1;
      $seen[$j]=1;
      $j  = ($pt == 1) ? $pt2[$j] : $pt1[$j];
      $pt = ($pt == 1) ? 2 : 1;
    }

    # duplicate element if we have a cycle 
    if ($j) {
      croak "messed up cycle!" unless $j-1 == $path[-1];
      unshift @path, $j-1 unless @path == 2; # 2 => duplicate basepair
    }

    push @{$self->{plist}}, [@path];
  }
  return $self->{plist};
}

=head2 explore_sequence_space()

This function needs to be called to initialize the subsequent sequence design. It reads
the previously computed dependency graph to (i) update the sequence constraint, (ii)
caclulate the total number of sequences able to fulfill sequence and structure 
constraints, (iii) set a stop-condition for sequence design dependent on the number
of possible mutation-moves, (iv) reduces the pathlist to a *slim* pathlist that only
contains mutateable elements.

=cut

sub explore_sequence_space {
  my $self = shift;
  my $verb = 0;

  my @plist =@{$self->{plist}};
  my $con   =  $self->{constraint};
  my %iupack=%{$self->{iupack}};

  croak "need to comupute depencendy pathways first!" unless @plist;
  croak "need constraint infromation!" unless $con;

  my ($border, $max_len, $nos) = (0,0,1);
  my @slim_plist;
  my @rand_plist = (0);
  foreach my $path (@plist) {
    my @pseq=();
    my $cycle=0;
    my $num;

    # Let's do the simple thingy separate
    if (@$path == 1) {
      # Get constraint
      $pseq[0] = substr($con, $path->[0], 1);

      # Statistics (bos, border)
      $num  = length $iupack{$pseq[0]};
    } else {
      $cycle=1 if ($path->[0] eq $path->[-1]);

      # Translate path into nucleotide constraints
      my $l = ($cycle) ? $#$path-1 : $#$path;
      push @pseq, substr($con, $$path[$_], 1) for (0 .. $l);

      # Update Constraint
      @pseq = $self->update_constraint($cycle, @pseq);

      # Enumerate all possible pathways 
      $num = $self->enumerate_pathways($cycle, @pseq);
    }

    if ($num > 1) {
      push @slim_plist, $path;# for (1 .. $num);
      push @rand_plist, $rand_plist[-1]+$num;
      $border += $num;
      $nos    *= $num;
    }

    #print "@$path => @pseq | Num: $num\n" if 1;
    substr($con, $path->[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  }
  shift @rand_plist;
  #print "@rand_plist\n";
  #print "SLP: ".@slim_plist." => $border\n";
  #print "$rand_plist[$_] => @{$slim_plist[$_]}\n" foreach (0 .. $#slim_plist);
  #print "@{$slim_plist[$_]}\n" foreach (0 .. $#slim_plist);

  $self->{plist}=\@slim_plist;
  $self->{rlist}=\@rand_plist;
  $self->{constraint}=$con;
  $self->{border}=$border unless $self->{border};
  $self->{nos}=$nos;
  return ($con, $border, $nos);
}

=head2 update_constraint()

Rewrite the constraint for a particular dependency path, such that starting an
arbitrary position would result in a correct path: "NNNRNNN" => "YRYRYRY";

=cut

sub update_constraint {
  my $self = shift;
  my $cycle = shift;
  my @pseq  = @_;

  croak "cycle of uneven length!" if $cycle && @pseq % 2;

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
    croak "escaping from seemingly endless constraint updates" if $jailbreak++ > 100;
  }
  return @pseq;
}

=head2 enumerate_pathways()

For a given constrained depencendy path, calculate the number of sequences
fulfilling that constraint. If it is a constrained path, exhaustivley enumerate
and store the solution-tree in a hash. If it is constrained and too long for 
exhaustive enumeration, estimate the solutions with the fibronacci numbers and
then (in make_pathseq()) shuffle them with a greedy heuristic.

=cut

sub enumerate_pathways {
  my $self = shift;
  my $cycle= shift;
  my @pseq = @_;

  my $pstr = join '', @pseq;
  my $plen = length $pstr;

  croak "cycles must have even length" if $cycle && $plen%2;
  croak "first call update_constraint(), then enumerate_pathways()!" if ($pstr =~ m/[N]/g) && ($pstr =~ m/[^N]/g);

  my %solutions = %{$self->{solution_space}};
  my $max_plen  = $self->{max_const_plen};
  my %iupack    = %{$self->{iupack}};
  my %base      = %{$self->{base}};
  my @pair      = @{$self->{pair}};

  if (exists $solutions{$pstr.$cycle}) {
    return scalar($solutions{$pstr.$cycle}->get_leaves);
  } elsif ($pstr =~ m/[N]/g) { 
    # its a path with all N's => (fibronacci: switch.pl)
    my $l = ($cycle) ? $plen-1 : $plen;
    return 2*($self->get_fibo($plen+1) + $self->get_fibo($l));
  } 
  # enable this as soon as fair sampling in this case works!
  # elsif ($pstr !~ m/[^RY]/g) {
  #   # its a path with all RY's
  #   my $l = ($cycle) ? $plen-1 : $plen;
  #   return ($self->get_fibo($plen+1) + $self->get_fibo($l));
  # } 
  #
  elsif ($plen > $max_plen) { 
    # path too long to enumerate with constraints
    carp "constrained path too long, fallback to greedy heuristic!";
    my $l = ($cycle) ? $plen-1 : $plen;
    return ($self->get_fibo($plen+1) + $self->get_fibo($l));
    # calculate RYRYR as NNNN/2 and use greedy-shuffle

    # compute the theoretical sequence length (n) to construct an 
    # example that breaks plen: n = (plen-1)*4 + 1 
    # so for plen=26 => n=101
    # if we assume a helix has 4 base-pairs, it is:
    # n = (plen-1)*12 + 1 = 301
  }

  my $tree = RNA::Design::Tree->new();
  my @le = $tree->get_leaves;
  $tree->init_new_leaves;
  foreach (split //, $iupack{$pseq[0]}) {
    $tree->push_to_current_leaves($_, $le[0]);
  }
  
  for my $i (1 .. $#pseq) {
    my @choices = split //, $iupack{$pseq[$i]};
    my @leaves  = $tree->get_leaves;
    $tree->init_new_leaves;

    if ($i == $#pseq && $cycle) {
      foreach my $l (@leaves) {
        foreach my $c (@choices) {
          if ($pair[$base{$$l[0]}][$base{$c}]) {
            # check if $c also pairs with anything in the first row!
            my $string = $tree->get_full_path($l);
            my $f = substr $string, 0, 1;
            if ($pair[$base{$c}][$base{$f}]) {
              $tree->push_to_current_leaves($c, $l);
            }

          }
        }
      }
    } else {
      foreach my $c (@choices) {
        foreach my $l (@leaves) {
          if ($pair[$base{$$l[0]}][$base{$c}]) {
            $tree->push_to_current_leaves($c, $l);
          }
        }
      }
    }
  }

  #DEBUG:
  #print $pstr.$cycle."\n";
  #foreach ($tree->get_leaves) {
  #  print $tree->get_full_path($_)."\n";
  #}
  
  $self->{solution_space}->{$pstr.$cycle} = $tree;
  return scalar($tree->get_leaves);
}

sub fill_solution_space {

}
 
=head2 rewrite_neighbor()

Takes a letter and a reference to the neighbor. Rewrites the neighbor and 
returns 0 or 1 if there exists a neighbor or not.

=cut

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

Uses the dependency graph and the sequence constraint to make a random sequence. 

=cut

sub find_a_sequence {
  # make random sequence or randomly mutate sequence
  # by replacing all @plists by random paths
  my $self = shift;

  my @plist=@{$self->{plist}};
  my $con  = $self->{constraint};
  my $seq  = $con;

  foreach my $path (@plist) {
    my @pseq = ();
    my $cycle=0;

    $cycle = 1 if ($$path[0] eq $$path[-1] && (@$path>1));

    # Translate path into nucleotide constraints
    my $l = ($cycle) ? $#$path-1 : $#$path;
    push @pseq, substr($con, $$path[$_], 1) for (0 .. $l);

    @pseq = $self->make_pathseq($cycle, @pseq);

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
      last if (++$reject >= 2*$border);
    }
  }
  return $refseq;
}

=head2 mutate_seq()

Choose a random cycle, mutate it using make_pathseq of the sequence constraint.
Choses randomly according to the number of solutions, but can be done more efficiently
in n*log(n) time!

=cut

sub mutate_seq {
  my $self  = shift;
  my $seq   = shift;

  my $con   = $self->{constraint};
  my @plist = @{$self->{plist}};
  my @rlist = @{$self->{rlist}};

  my $reject_same_solution=1;

  my ($path, @pseq, @ref_pseq);
  my $cycle = 0;

  # chose a random path => weight by number of solutions to make it fair!
  my $rand = (int rand $rlist[-1]) + 1;
  for (my $i=0; $i<=$#plist; ++$i) {
    if ($rand <= $rlist[$i]) {
      $path = \@{$plist[$i]};
      last;
    }
  }

  # this would be random independent of the pathlength!
  #$path = $plist[int rand @plist];
  #print "@$path\n";

  $cycle = 1 if ($$path[0] eq $$path[-1] && (@$path>1));

  my $l = ($cycle) ? $#$path-1 : $#$path;
  push @pseq, substr($con, $$path[$_], 1) for (0 .. $l);

  if ($reject_same_solution) {
    push @ref_pseq, substr($seq, $$path[$_], 1) for (0 .. $l);
    my @npseq;
    my $jailbreak = 0;
    do {
      @npseq = $self->make_pathseq($cycle, @pseq);
      die if ++$jailbreak > 100;
    } while (join('', @ref_pseq) eq join('', @npseq));
    @pseq = @npseq;
  } else {
    @pseq = $self->make_pathseq($cycle, @pseq);
  }

  substr($seq, $$path[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  return $seq;
}

=head2 make_pathseq()

Takes an Array of IUPACK code and rewrites it into a valid array of Nucleotides. There
are four cases: (i) path of length 1 is a randomly shuffled nucleotide accoriding to 
iupack-letter, (ii) The solution-tree has been built before and so we save some work,
(iii) the sequence is only N's, so fibronacci is used (i.e. switch.pl method), and
(iv) the path is constrained and longer than *max_const_plen* so the solution space
is counted by fibronacci and paths are shuffled with a greedy heuristic.

=cut

sub make_pathseq {
  my $self = shift;
  my $cycle= shift;
  my @pseq = @_; 

  my $pstr = join '', @pseq;
  my $plen = length $pstr;

  croak "cycles must have even length" if $cycle && $plen%2;

  my %iupack    = %{$self->{iupack}};
  my %iupack_bin= %{$self->{iupack_bin}};
  my %solutions = %{$self->{solution_space}};
  my $max_plen  =   $self->{max_const_plen};

  if ($plen == 1) {
    # if path length 1 => return random iupack
    my @i = split '', $iupack{$pstr};
    return ($i[int rand @i]);
  } elsif (exists $solutions{$pstr.$cycle}) {
    my $tree = $solutions{$pstr.$cycle};
    my @i = $tree->get_leaves;
    return (split '', $tree->get_full_path($i[int rand @i]));
  } elsif ($pstr !~ m/[^N]/g) { 
    # its a path with all N's => do it like switch.pl
    
    #TODO: need to implement this method for RY-pathways as well
    
    # This is original switch.pl code that I don't understand!
    my $l = $plen;
    my $ll = ($cycle) ? $l-1 : $l;
    my $n = 2*($self->get_fibo($l+1) + $self->get_fibo($ll));

    my $rand = rand($n);

    my @seq = ();
    # set first base in sequence
    if ($rand < $self->get_fibo($ll)) {
      push @seq, 'A', 'U';
    } elsif ($rand < 2*$self->get_fibo($ll)) {
      push @seq, 'C', 'G'; 
      $rand -= $self->get_fibo($ll);
    } else {
      $rand -= 2*$self->get_fibo($ll); 
      $ll=$l;
      push @seq, ($rand>=$self->get_fibo($l+1))?'U':'G';
      $rand -= $self->get_fibo($l+1) if $rand >= $self->get_fibo($l+1);
    }

    # grow sequence to length $l
    # if we have a cycle starting with A or C $ll=$l-1, else $ll=$l
    while (@seq < $l) {
      if ($rand < $self->get_fibo($ll-@seq)) {
        push @seq, 'C','G' if ($seq[-1] eq 'G');
        push @seq, 'A','U' if ($seq[-1] eq 'U');
      } else {
        $rand -= $self->get_fibo($ll-@seq);
        push @seq, ($seq[-1] eq 'G') ? 'U' : 'G';
      }
    }
    pop @seq if (@seq > $l); # in case we've added one base too many
    return @seq;
  } else {
    # do a greedy constraint shuffling!
    croak "some special case not covered yet! (@pseq)" unless $plen > $max_plen;

    my @path = @pseq;
    for (my $i=0; $i<=$#path; ++$i) {
      my $c = $path[$i];
      my @i = split '', $iupack{$c};

      $path[$i] = $i[int rand @i];

      if ($i==0 && $cycle) {
        my ($a, $j) = ($path[$i], -1);
        $a=$path[$j--] while $self->rewrite_neighbor($a, \$path[$j]);
      }

      $self->rewrite_neighbor($path[$i], \$path[$i+1]) if $i < $#path;
    }
    return @path;
  }
  return;
}

=head2 eval_sequence()

Evaluate a given sequence according to your cost function.

=cut

sub eval_sequence {
  my $self = shift;
  my $seq  = shift;

  my $verb = $self->{verb};
  my $ofun = $self->{optfunc};
  my $apen = $self->{avoid_penalty};
  warn "$seq\n".$ofun."\n" if $verb > 1;
  croak "no structures specified!\n" unless @{$self->{structures}};

  my %results;
  $RNA::fold_constrained=1;
  foreach my $func (qw/eos eos_circ pfc pfc_circ gfe gfe_circ prob prob_circ/) {
    while ($ofun =~ m/$func\(([0-9\,\s]*)\)/) {
      warn "next: $&\n" if $verb > 1;
      $results{$&} = eval "\$self->$func(\$seq, $1)" unless exists $results{$&};
      croak $@ if $@;
      $ofun =~ s/$func\(([0-9\,\s]*)\)/ $results{$&} /;
    }
  }
  $RNA::fold_constrained=0;

  warn $ofun."\n" if $verb > 1;
  my $r = eval $ofun;
  croak $@ if $@;

  my $p = $self->base_prob($seq);
  my $a = $self->base_avoid($seq, $apen);

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

  my %prob = %{$self->{base_probs}};
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
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  croak "cannot find structure number $j" unless $self->{structures}[$j-1];
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
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  croak "cannot find structure number $j" unless $self->{structures}[$j-1];
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
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  return RNA::energy_of_struct($seq, $str);
}

sub eos_circ {
  my ($self, $seq, $i, $t) = @_;
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  return RNA::energy_of_circ_struct($seq, $str);
}

sub pfc {
  my ($self, $seq, $i, $t) = @_;
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
  # seemingly useless string modification
  # to avoid ViennaRNA SWIG interface bugs
  $str.='.'; $str=substr($str,0,-1);
  return RNA::pf_fold($seq, $str);
}

sub pfc_circ {
  my ($self, $seq, $i, $t) = @_;
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $str = $self->{structures}[$i-1];
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

1;

package RNA::Design::Tree;

use strict;
use warnings;
use Data::Dumper;

sub new {
  my $class = shift;
  my $self  = {
    tree => [
      [
        ['root', undef]
      ]
    ],
  };

  bless $self, $class;
  return $self;
}

sub init_new_leaves {
  my $self = shift;
  push @{$self->{tree}}, [];
  return $self->{tree};
}

sub push_to_current_leaves {
  my $self = shift;
  my ($string, $parent) = @_;

  my $tree = $self->{tree};
  my $leaves = $$tree[-1];

  push @$leaves, [$string, $parent];
  
  return $self->{tree};
}

sub get_leaves {
  my $self = shift;
  return @{${$self->{tree}}[-1]};
}

sub get_first {
  my $self = shift;
  return @{${$self->{tree}}[1]};
}

sub get_full_path {
  my $self  = shift;
  my $leave = shift;

  my $string;
  my ($char, $ref) = @$leave;
  while ($char ne 'root') {
    $string .= $char;
    ($char, $ref) = @$ref;
  }
  return reverse $string;
}

1;

