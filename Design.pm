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

B<RNA::Design> -- plug-and-play design of nucleic acid sequences

=head1 DESCRIPTION

This package provides various subroutines to design RNA molecules optimized for
single or multiple conformations. The properties of the target molecules have
to be specified in the Main Object: C<$Design = RNA::Design-E<gt>new()>.
B<RNA::Design> supports structure constraints that can be specified by two
well-formed dot-bracket strings. Additional base-pairs specified as input will
be ignored and produce a warning. Sequence optimization functions can be
composed of the sub-functions: C<eos(i,t)> (energy of structure), C<efe(i,t)>
(ensemble free energy), and C<prob(i,j,t)> (probability of structure), where
C<i> and C<j> stand for the number of the input-structure and C<t> allows the 
user to specify a temperature. All of these functions exist for circular sequences 
by attatching a C<_circ> flag (for example C<eos_circ(i,t)>.

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
    fibo  => [0,1],
    border=> 0,
    nos   => undef,
    plist => [],
    rlist => [],
    avoid => ['AAAAA','CCCCC','GGGGG','UUUUU'],
    avoid_penalty => 5,
    base_probs    => ({A => 0.25, C => 0.25, G => 0.25, U => 0.25}),
    verb  => 0,

    structures  => [],
    cut_point   => -1,
    constraint  => '',
    optfunc     => 'eos(1)+eos(2) - 2*efe() + 0.3*(eos(1)-eos(2)+0.00)**2',

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

    iupac => ({
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
    iupac_bin => ({ # ACGU
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
    bin_iupac => [  # ACGU
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

Get and Set parameters for sequence design. For every C<set> routine, there
exists the equivalent C<get> routine. 

=cut

=head3 set_verbosity(<INT>)

set the verboisty of warnings to STDERR.

=cut 

sub set_verbosity {
  my ($self, $var) = @_;
  $self->{verb}=$var;
  return $self->{verb};
}

sub get_verbosity {
  my $self = shift;
  return $self->{verb};
}

=head3 set_optfunc(<STRING>)

set the cost-function for sequence optimization.

=cut

sub set_optfunc {
  my ($self, $var) = @_;
  $self->{optfunc} = $var;
  return $self->{optfunc};
}

sub get_optfunc {
  my $self = shift;
  return $self->{optfunc};
}

=head3 set_constraint(<STRING>)

set the sequence constraints for sequence optimization.

=cut

sub set_constraint {
  my ($self, $var) = @_;
  $self->{constraint} = $var;
  return $self->{constraint};
}

sub get_constraint {
  my $self = shift;
  return $self->{constraint};
}

=head3 set_avoid_motifs(<ARRAY>)

set a list of sequence motifs that will cause a penalty during optimization.

=cut

sub set_avoid_motifs {
  my $self = shift;
  $self->{avoid} = [@_];
  return $self->{avoid};
}

sub get_avoid_motifs {
  my $self = shift;
  return $self->{avoid};
}

=head3 set_avoid_penalty(<INT>)

set the penalty for sequence motifs specified with C<set_avoid_motifs()>

=cut

sub set_avoid_penalty {
  my ($self, $var) = @_;
  $self->{avoid_penalty} = $var;
  return $self->{avoid_penalty};
}

sub get_avoid_penalty {
  my $self = shift;
  return $self->{avoid_penalty};
}

=head3 set_base_probs(<HASH>)

set the penalty for sequence motifs specified with C<set_avoid_motifs()>

=cut

sub set_base_probs {
  my ($self, $var) = @_;
  $self->{base_probs} = $var;
  return $self->{base_probs};
}

sub get_base_probs {
  my $self = shift;
  return $self->{base_probs};
}

=head3 set_structures(<ARRAY>)

set the list of structures for sequence optimization

=cut

sub set_structures {
  my $self = shift;
  my @structs = @_;
  if ($self->{verb}) {
    foreach (@structs) {
      croak "only dot-bracket strings allowed in add_structures()" if m/[^\(\)\.]/g; 
    }
  }
  $self->{structures} = [@structs];
  return $self->{structures};
}

sub add_structures {
  my $self = shift;
  foreach (@_) {
    if ($self->{verb} && m/[^\(\)\.]/g) { 
      croak "only dot-bracket strings allowed in add_structures()";
    }
    push @{$self->{structures}}, $_;
  }
  return $self->{structures};
}

sub get_structures {
  my $self = shift;
  return $self->{structures};
}

=head3 set_cut_point(<INT>)

set a cut_point when designing two interacting molecules. The cut point
indicates the first nucleotide of the second sequence.

=cut

sub set_cut_point {
  my ($self, $var) = @_;
  if ($self->{verb} && $self->{cut_point} != -1 && $self->{cut_point} != $var) {
    carp "overwriting old cut_point";
  }
  $self->{cut_point} = $var;
  return $self->{cut_point};
}

sub get_cut_point {
  my $self = shift;
  return $self->{cut_point};
}

=head3 get_num_of_seqs()

get the number of sequences after calling explore_sequence_space()

=cut

sub get_num_of_seqs {
  my $self = shift;
  return $self->{nos};
}

=head3 get_num_of_nbors()

get the number of neighbors of a sequence after calling explore_sequence_space()

=cut

sub get_num_of_nbors {
  my $self = shift;
  return $self->{border};
}


=head3 get_fibo(<INT>)

get the fibronaccy number at position <INT>. The list starts with [0,1] at 
positions 0,1. Note that there is no C<set> routine for fibonacci.

=cut

sub get_fibo {
  my ($self, $var) = @_;
  $self->fill_fibo($var) unless defined $self->{fibo}[$var];
  return $self->{fibo}[$var];
}

sub fill_fibo {
  my ($self, $var) = @_;
  while ($#{$self->{fibo}} < $var) {
    push @{$self->{fibo}}, $self->{fibo}->[-1] + $self->{fibo}->[-2];
  }
  return $self->{fibo};
}

=head2 find_dependency_paths(@structures)

If @structures is empty, the mutation dependency graph is constructed from all
structures added with the set_structures() routine. Specify @structures only if
you do not want all of your structures for the cost-function to be part of the
depencency graph.

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

This function initializes internal paramters for subsequent sequence
optimization. It reads the previously computed dependency graph to 

=over 

=item 1. 

update and check the sequence constraint 
  
=item 2.

caclulate the total number of sequences compatible with all constraints
  
=item 3.

set a stop-condition for optimization from the number of possible mutations
  
=item 4.

store the number of possible solutions per depencendy-path for subsequent fair
sampling.

=back

=cut

sub explore_sequence_space {
  my $self = shift;
  my $verb = 0;

  my @plist =@{$self->{plist}};
  my $con   =  $self->{constraint};
  my %iupac =%{$self->{iupac}};

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
      $num  = length $iupac{$pseq[0]};
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
      push @rand_plist, $rand_plist[-1]+($num-1);
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
  
  print "$nos Number of sequences, $border neighbors for each sequence.\n" if $self->{verb};

  $self->{plist}=\@slim_plist;
  $self->{rlist}=\@rand_plist;
  $self->{constraint}=$con;
  $self->{border}=$border unless $self->{border};
  $self->{nos}=$nos;
  return ($con, $border, $nos);
}

=head3 update_constraint($is_cycle, @path)

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

=head3 enumerate_pathways($is_cycle, @path)

For a given depencendy path, calculate the number of compatible sequences, and
initialze data-structures for C<make_pathseq()>. If it is an unconstrained
path, use the fibronacci numbers to enumerate pathways. If it is a constrained
path, exhaustivley count and store the solutions in a tree structure. If the
path is constrained and longer than C<max_const_plen>, estimate the number of
pathways with the fibronacci numbers.

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
  my %iupac     = %{$self->{iupac}};
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
  foreach (split //, $iupac{$pseq[0]}) {
    $tree->push_to_current_leaves($_, $le[0]);
  }
  
  for my $i (1 .. $#pseq) {
    my @choices = split //, $iupac{$pseq[$i]};
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
 
=head3 rewrite_neighbor($iupac, \$iupac)

Takes a IUPAC letter and a reference to a base-pairing IUPAC letter. Updates
the IUPAC code if necessary and returns 0 or 1 if there exists a neighbor or
not.

=cut

sub rewrite_neighbor {
  my ($self, $c1, $c2) = @_;
  print "$c1 -> $$c2\n" if 0;

  my $nb = $self->{neighbor}->{$c1};
  if ($nb eq $$c2) {
    return 0;
  } else {
    my $bin_c2 = $self->{iupac_bin}{$$c2};
    my $bin_nb = $self->{iupac_bin}{$nb};

    my $new_nb = $self->{bin_iupac}[($bin_c2 & $bin_nb)];
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
better score from the cost-function, replaces the current solution. In case
there are too many useless mutations, (see explore_sequence_space()), the
current sequence is considered local-optimal and returned.

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

=head3 mutate_seq($sequence)

Choose a random depencency path, mutate it using make_pathseq() of the sequence
constraint. Choses randomly according to the number of solutions, but can be
implemented more efficiently in C<log(n)> time!

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
  return $seq unless defined $rlist[-1];
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
      last if ++$jailbreak > 100;
    } while (join('', @ref_pseq) eq join('', @npseq));
    @pseq = @npseq;
  } else {
    @pseq = $self->make_pathseq($cycle, @pseq);
  }

  substr($seq, $$path[$_], 1, $pseq[$_]) for (0 .. $#pseq);
  return $seq;
}

=head3 make_pathseq($is_cycle, @path)

Takes an Array of IUPAC code and rewrites it into a valid array of Nucleotides.
There are four cases: (i) path of length 1 is a randomly shuffled nucleotide
accoriding to IUPAC-letter, (ii) The solution-tree has been built during
explore_sequence_space() and so we use it to sample, (iii) the sequence is only
N's, so the fibronacci numbers are used for sampling (i.e. the switch.pl
method), and (iv) the path is constrained and longer than C<max_const_plen>, so
the paths are sampled using a greedy heuristic.

=cut

sub make_pathseq {
  my $self = shift;
  my $cycle= shift;
  my @pseq = @_; 

  my $pstr = join '', @pseq;
  my $plen = length $pstr;

  croak "cycles must have even length" if $cycle && $plen%2;

  my %iupac     = %{$self->{iupac}};
  my %iupac_bin = %{$self->{iupac_bin}};
  my %solutions = %{$self->{solution_space}};
  my $max_plen  =   $self->{max_const_plen};

  if ($plen == 1) {
    # if path length 1 => return random iupac
    my @i = split '', $iupac{$pstr};
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
      my @i = split '', $iupac{$c};

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

=head2 eval_sequence($sequence)

Evaluate a given sequence according to the cost function.

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
  $RNA::cut_point=$self->{cut_point};
  foreach my $func (qw/eos eos_circ efe efe_circ prob prob_circ/) {
    while ($ofun =~ m/$func\(([0-9\,\s]*)\)/) {
      warn "next: $&\n" if $verb > 1;

      if (exists $results{$&}) {
        #print "OR: $& => $results{$&}\n";
      } else {
        my $res = eval("\$self->$func(\$seq, $1)");
        $res = sprintf "%.2f", $res;
        croak $@ if $@;
        #print "NR: $& => $res\n";
        $results{$&}=$res;
      }
      $ofun =~ s/$func\(([0-9\,\s]*)\)/ $results{$&} /;
    }
  }
  $RNA::cut_point=-1;
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
  my $cpnt = $self->{cut_point};

  if ($cpnt == -1) {
    foreach my $string (@avoid) {
      $penalty *= $pen while ($seq =~ m/$string/g);
    }
  } else {
    my $left  = substr $seq, 0, $cpnt-1;
    my $right = substr $seq, $cpnt-1;
    #print "$left&$right\n";
    foreach my $string (@avoid) {
      $penalty *= $pen while ($left   =~ m/$string/g);
      $penalty *= $pen while ($right  =~ m/$string/g);
    }
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
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $kT=0.6163207755;

  croak "cannot find structure number $i" if $i && !$self->{structures}[$i-1];
  croak "cannot find structure number $j" if $i && !$self->{structures}[$j-1];

  my $s_i = ($i) ? $self->{structures}[$i-1] : undef;
  my $s_j = ($j) ? $self->{structures}[$j-1] : undef;

  my ($dGi, $dGj, $tmp);
  if ($RNA::cut_point == -1) {
    $tmp=$s_i; $dGi = RNA::pf_fold($seq, $tmp);
    $tmp=$s_j; $dGj = RNA::pf_fold($seq, $tmp);
  } else {
    $tmp=$s_i; $dGi = RNA::co_pf_fold($seq, $tmp);
    $tmp=$s_j; $dGj = RNA::co_pf_fold($seq, $tmp);
  }

  return exp(($dGj-$dGi)/$kT);
}

sub prob_circ {
  my ($self, $seq, $i, $j, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  my $kT=0.6163207755;

  croak "cannot find structure number $i" if $i && !$self->{structures}[$i-1];
  croak "cannot find structure number $j" if $j && !$self->{structures}[$j-1];
  carp "ignoring cut_point for circular costfunction" if ($RNA::cut_point != -1);

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
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  my $str = $self->{structures}[$i-1];
  return RNA::energy_of_struct($seq, $str);
}

sub eos_circ {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  croak "cannot find structure number $i" unless $self->{structures}[$i-1];
  my $str = $self->{structures}[$i-1];
  return RNA::energy_of_circ_struct($seq, $str);
}

sub efe {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  croak "cannot find structure number $i" if $i && !$self->{structures}[$i-1];

  my $str = ($i) ? $self->{structures}[$i-1] : undef;
  # seemingly useless string modification
  # to avoid ViennaRNA SWIG interface bugs
  if ($str) {
    $str.='.'; $str=substr($str,0,-1);
  }
  return ($RNA::cut_point == -1) ? RNA::pf_fold($seq, $str) : RNA::co_pf_fold($seq, $str);
}

sub efe_circ {
  my ($self, $seq, $i, $t) = @_;
  $t = 37 unless defined $t;
  $RNA::temperature=$t;
  croak "cannot find structure number $i" if $i && !$self->{structures}[$i-1];
  carp "ignoring cut_point for circular costfunction" if ($RNA::cut_point != -1);

  my $str = ($i) ? $self->{structures}[$i-1] : undef;
  # seemingly useless string modification
  # to avoid ViennaRNA SWIG interface bugs
  if ($str) {
    $str.='.'; $str=substr($str,0,-1);
  }
  return RNA::pf_circ_fold($seq, $str);
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

