package Range;

use Exporter;

our @ISA=qw(Exporter);
our @EXPORT=qw(cmpRange);

=head1 NAME

Range - Range with direction, used for sequence analysis

=head1 DESCRIPTION

start <= end, use direction to indiction direction

=head2 new($start, $end, $direction)

   Alternative input ($start,end), where if $start>$end
   direction of - is implied.

this class is a reference to an array
this can be called either from a instance or 
as a class method

$direction is either + or -

The underlying data structure is an array

   0     1       2
[start, end, direction]

start <= end, direction (+/-) tells the 5'->3' direction
of the feature.

   If direction is not given and $start,$end is directional
   range. For example [3,5] indicate + and [50,10] indicate
   -, then this constructor will also work to fill in the
   direction parameter.

=cut

sub new {
    my $invocant = shift;
    my $class = ref $invocant ? ref $invocant : $invocant;
    my $self = [@_];
    if (@$self == 2) {
        if ($self->[0] < $self->[1]) { $self->[2]='+'; }
        else { 
            my $tmp=$self->[0];
            $self->[0]=$self->[1];
            $self->[1]=$tmp;
            $self->[2]='-'; 
        }
    }

    bless $self, $class;
    return $self;
}

=head2 copy()

Make a new copy of the original object

=cut

sub copy {
    my $self=shift;
    return Range->new($self->[0], $self->[1], $self->[2]);
}

=head2 direction or getDirection

return the direction of the range + or - or ' '

=cut

sub direction {
    $self=shift;
    unless (defined $self->[2]) {
        die "direction not defined!\n";
    }
    return $self->[2];
}
sub getDirection {
    $self=shift;
    unless (defined $self->[2]) {
        die "direction not defined!\n";
    }
    return $self->[2];
}

=head2 length

return the length of the range

=cut

sub length {
    $self=shift;
    return $self->[1]-$self->[0]+1;
}

sub empty {
    my $self=shift;
    return scalar(@$self) < 1;
}

=head2 overlap($another_range)

return the overlap size of the two ranges only if they
are in the same direction, otherwise return negative overlap

return 0 of no overlap.

=cut

sub overlap {
    my $r1=shift; # is self
    my $r2=shift;
# smaller end - larger begin
    my $olp=($r1->[1] < $r2->[1] ? $r1->[1] : $r2->[1]) - 
            ($r1->[0] > $r2->[0] ? $r1->[0] : $r2->[0]);
    if ($r1->direction() eq $r2->direction()) {
        if ($olp>0) { return $olp; }
        else { return 0; }
    }
    else {
      if ($olp>0) { return -$olp; }
      else { return 0; }
    }
}

sub containPoint {
   my $self=shift;
   my $p=shift;
   if ($p >= $self->[0] && $p <= $self->[1]) {
      return 1;
   }
   return 0;
}

=head2 overlapFraction($another_range)

return an array of two elements containing the 
fraction of the overlap relative to the length of 
this range and the other range.

=cut

sub overlapFraction {
    my $self=shift;
    my $r2=shift;
    my $olp=$self->overlap($r2);
    return ($olp/$self->length, $olp/$r2->length);
}

=head2 overlapFractionMoreThan($anotherRange, $cutoff)

    Return the length of the overlap if the overlap as computed with respect to
       either one of the sequences is more than the $cutoff.
       else return 0.

test the overlap is greater than the cutoff value
The cutoff value should be number smaller than 1.

=cut

sub overlapFractionMoreThan {
    my $self=shift;
    my $r2=shift;
    my $cutoff=shift;
    my $olp=$self->overlap($r2);
    if ($olp/$self->length > $cutoff || $olp/$r2->length > $cutoff) {
        return $olp;
    }
    else { return 0; }
}

=head2 merge

combine two ranges so that this range contains both
of the two ranges before the merge.

Note: this object is enlarged.

The caller is responsible to check for the directions.

=cut

sub merge {
    my $self=shift;
    my $r2=shift;
    $self->[0] = $self->[0] < $r2->[0] ? $self->[0] : $r2->[0];
    $self->[1] = $self->[1] > $r2->[1] ? $self->[1] : $r2->[1];
}

=head2 getBegin

=cut
sub getBegin {
    my $self=shift;
    return $self->[0];
}


=head2 getEnd

=cut
sub getEnd {
    my $self=shift;
    return $self->[1];
}

=head2 setBounds($b,$e)

set new bounds

=cut
sub setBounds {
    my $self=shift;
    $self->[0]=$_[0]; $self->[1]=$_[1];
}

=head2 less($other_range)

compare this range to $other_range
regardless of direction. 

return true if left of the other range

=cut
sub less {
    my $self=shift;
    my $other=shift;
    if ($self->[0] < $other->[0]) { return 1; }
    elsif ($self->[0] == $other->[0]) {
        if ($self->[1] < $other->[1]) { return 1; }
        else { return 0; }
    }
    else { return 0; }
}

=head2 lessDirection

    + before -

=cut
sub lessDirection {
    my $self=shift;
    my $other=shift;
    if ($self->direction eq '+' && $other->direction eq '-') {
        return 1;
    }
    elsif ($self->direction eq '-' && $other->direction eq '+') {
        return 0;
    }

    if ($self->[0] < $other->[0]) { return 1; }
    elsif ($self->[0] == $other->[0]) {
        if ($self->[1] < $other->[1]) { return 1; }
        else { return 0; }
    }
    else { return 0; }
}

=head2 cmpRange

This is not an object method, needs to be exported

    return 1 Range1 > Range2
           0 equal
           -1 Range1 < Range2

=cut
sub cmpRange {
    my $r1=shift;
    my $r2=shift;
    if ($r1->lessDirection($r2)) { return -1; }
    elsif ($r2->lessDirection($r1)) {  return 1; }
    else {  return 0; }
}

=head2 show

    @param $fh output file stream

    show in a human readable format: [B,E]+/-
    print to output stream.

=cut
sub show {
    my $self=shift;
    my $fh=shift;
    if (!$fh) { $fh=\*STDOUT; }
    print $fh '[', $self->[0], '-', $self->[1], ']', $self->[2];
}

1;
