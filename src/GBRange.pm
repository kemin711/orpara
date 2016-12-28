package GBRange;

use base("Range");

=head1 NAME 

GBRange - a derived class of Range, with property
attribute to remember whether it has start or ends

The underlying data structure is an array
that extended from the parnet class Range

[start, end, direction, has_start, has_end ] 

has_start referes to start,
has_end referes to end, always from small to big

This class is mainly used for GenBank Submission

this is an derived  class, it can have extra
arguments, it will use the base class 
new method.

GBRange->new(10, 100, '+', 1, 1)
if not given element 3,4 will be undef
   which will be interpreted as true.

=cut


=head2 nobegin($with_begin)

$with_begin is a boolean variable 0 or 1.

nobegin(0) means there is nobegin.
This methods is misleading. But old code has been using it this way.

=cut

sub nobegin {
    my $self=shift;
    my $nb=shift;
    $self->[3]=$nb;
}

sub noend {
    my $self=shift;
    my $ne=shift;
    $self->[4]=$ne;
}

=head2 setStart($has_start)

should rename this method setHasStart
I was too lazy in typin the extra letters.

=cut
sub setStart {
    my $self=shift;
    $self->[3]=$_[0];
}

=head2 setHasStart

=cut
sub setHasStart {
    my $self=shift;
    $self->[3]=$_[0];
}
=head2 setFinish($has_stop)

=cut
sub setFinish {
    my $self=shift;
    $self->[4]=$_[0];
}
=head2 setStop

alias for setFinish

=cut
sub setStop {
    my $self=shift;
    $self->[4]=$_[0];
}

=head2 setHasStop

This should be the properly named method. Use this one in the
future.
This method replaces setStop()

=cut
sub setHasStop {
    my $self=shift;
    $self->[4]=$_[0];
}

=head2 hasStart

return the gene has Start, either it 
has the start codon for CDS or has 5'-UTR for mRNA

This methods can also be used to set the internal vale
    hasStart(0) will mean set not start.

=cut

sub hasStart {
    my $self=shift;
    if (@_>0) { $self->[3]=$_[0]; }
    if (!defined $self->[3]) { return 1; }
    return $self->[3];
}

=head2 hasStop

    Can be used to get the states and set the internal state
    hasStop(1) will set has stop to true.

=cut
sub hasStop {
    my $self=shift;
    if (@_>0) { $self->[4]=$_[0]; }
    if (!defined $self->[4]) { return 1; }
    return $self->[4];
}

=head2 output($fh)

write the range to $fh

the start,end value is always from small to large.
the direction is indicated by strand.
completeness is stored without regard to direction.

=cut

sub output {
    my $self=shift;
    my $fh=shift;
    if ($self->direction eq '+') {
        unless ($self->hasStart) { print $fh '<'; }
        print $fh $self->[0], "\t";
        unless ($self->hasStop) { print $fh '>'; }
        print $fh $self->[1];
    }
    elsif ($self->direction eq '-') {
        unless ($self->hasStart) { print $fh '<'; }
        print $fh $self->[1], "\t";
        unless ($self->hasStop) { print $fh '>'; }
        print $fh $self->[0];
    }
    else {
        die "direction unknown for GBRange\n";
    }
}

1;
