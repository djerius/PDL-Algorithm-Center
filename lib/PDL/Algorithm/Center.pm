package PDL::Algorithm::Center;

# ABSTRACT: Various methods of finding the center of a sample

use strict;
use warnings;

require 5.010000;

use feature 'state';

our $VERSION = '0.05';

use Carp;

use Safe::Isa;
use Ref::Util qw< is_arrayref is_ref is_coderef  >;

use custom::failures;
use Package::Stash;
use Hash::Wrap;

use PDL::Algorithm::Center::Failure;

use PDL::Algorithm::Center::Types -all;
use Types::Standard -types;
use Types::Common::Numeric -types;
use Type::Params qw[ compile_named ];


use PDL::Lite ();

use Exporter 'import';

our @EXPORT_OK = qw[ sigma_clip iterate ];


sub _weighted_mean_center {
    my ( $coords, $wmask, $weight ) = @_;

    $weight //= $wmask->dsum;

    return ( $coords * $wmask->dummy( 0 ) )->xchg( 0, 1 )->dsumover
      / $weight;
}

sub _distance {

    my ( $last, $current ) = @_;

    return sqrt( ( ( $last->center - $current->center )**2 )->dsum );
}

sub _sigma_clip_initialize {

    my ( $init_clip, $dtol, $coords, $wmask, $current, $work ) = @_;

    $current->{clip} = $init_clip;

    my $r2 = $work->{r2} = PDL->null;
    $r2 .= ( ( $coords - $current->center )**2 )->dsumover;


    if ( defined $init_clip ) {
        $wmask = $wmask->copy;
        $wmask *= $r2 <= $init_clip**2;
    }

    $current->weight( $wmask->dsum );
    $current->nelem( ( $wmask > 0 )->dsum );
    $current->{sigma} = sqrt( ( $wmask * $r2 )->dsum / $current->weight );

    return;
}

sub _sigma_clip_wmask {

    my ( $nsigma, $coords, $wmask, $iter, $work ) = @_;

    my $r2 = $work->{r2};
    $r2 .= ( ( $coords - $iter->center )**2 )->dsumover;

    $iter->clip( $nsigma * $iter->sigma );
    $wmask *= $r2 < $iter->clip**2;

    $iter->weight( $wmask->dsum );
    $iter->nelem( ( $wmask > 0 )->dsum );

    $iter->sigma( sqrt( ( $wmask * $r2 )->dsum / $iter->weight ) );

    return;
}


sub _sigma_clip_is_converged {

    my ( $init_clip, $dtol, $coords, $wmask, $last, $current ) = @_;

    $current->{dist} = undef;

    # stop if standard deviations and centers haven't changed

    if ( $current->sigma == $last->sigma
         && PDL::all( $current->center == $last->center) ) {

        $current->dist( _distance( $last, $current ) )
          if defined $dtol;

        return 1;
    }

    # or, if a tolerance was defined, stop if distance from old
    # to new centers is less than the tolerance.
    if ( defined $dtol ) {
        $current->dist( _distance( $last, $current ) );
        return  1 if $current->dist <= $dtol;
    }

    return;
}


## no critic (ProhibitAccessOfPrivateData)

=pod

=sub sigma_clip

  $results = sigma_clip(
    center    => Optional [ Center | CodeRef ],
    clip      => Optional [PositiveNum],
    coords    => Optional [Coords],
    dtol      => PositiveNum,
    iterlim   => Optional [PositiveInt],
    log       => Optional [CodeRef],
    mask      => Optional [Undef | Piddle_min1D_ne],
    nsigma    => PositiveNum,
    weight    => Optional [Undef | Piddle_min1D_ne],
 );

Center a dataset by iteratively excluding data outside of a radius
equal to a specified number of standard deviations. The dataset may be
specified as a list of coordinates and optional weights, or as a
weight piddle of shape I<NxM> (e.g., an image).  If only the weight
piddle is provided, it is converted internally into a list of
coordinates with associated weights.

The center of a data set is determined by:

=over

=item 1

clipping (ignoring) the data whose distance to the current center is
greater than a specified number of standard deviations

=item 2

calculating a new center by performing a (weighted) centroid of the
remaining data

=item 3

calculating the standard deviation of the distance from the remaining
data to the center

=item 4

repeat step 1 until either a convergence tolerance has been met or
the iteration limit has been exceeded

=back

The initial center may be explicitly specified,  or may be calculated
by performing a (weighted) centroid of the data.

The initial standard deviation is calculated using the initial center and either
the entire dataset, or from a clipped region about the initial center.

=head3 Options

The following options are available:

=over

=item C<center> => I<ArrayRef | Piddle1D_ne | coderef >

The initial center.  It may be

=over

=item *

An array of length I<N>

The array may contain undefined values for each dimension for which the center should
be determined by finding the mean of the values in that dimension.

=item *

a piddle with shape I<N>  (or
something that can be coerced into one, see L</TYPES>),

=item *

a coderef which will return the center as a piddle with shape I<N>.
The subroutine is called as

  &$center( $coords, $wmask, $weight );

with

=over

=item C<$coords>

a piddle with shape I<NxM> containing I<M> coordinates with dimension I<N>

=item C<$wmask>

a piddle of shape I<M> containing weights

=item C<$weight>

a scalar which is the sum of the weights in C<$wmask>

=back

=back

=item C<clip> => I<positive number>

I<Optional>.  The clipping radius used to determine the initial standard deviation.

=item C<coords> => I<Coords>

I<Optional>.  The coordinates to center.  C<coords> is a piddle of
shape I<NxM> (or anything which can be coerced into it, see
L</TYPES>) where I<N> is the number of dimensions in the data and
I<M> is the number of data elements.

C<weight> may be specified with coords to indicate weighted data.

C<coords> is useful if the data cube is not fully populated; for dense
data, use C<weight> I<instead>.

=item C<dtol> => I<positive number>

I<Optional>.  If specified iteration will cease when successive centers are closer
than the specified distance.

=item C<iterlim> => I<positive integer>

I<Optional>. The maximum number of iterations to run.  Defaults to 10.

=item C<log> => I<coderef>

I<Optional>. A subroutine which will be called before the first iteration and at
the end of each iteration. It is passed a copy of the current
iteration's results object; see L</Sigma Clip Iteration Results>.

=item C<mask> => I<piddle>

I<Optional>. This specifies data elements to ignore completely.
True values indicate elements to be used, false those to be ignored.

When used with C<coords>, C<mask> must be a piddle of shape I<M>,
where I<M> is the number of data elements in C<coords>.

If C<coords> is not specified, C<mask> should have the same shape as
C<weight>.

=item C<nsigma> => I<scalar>

The size of the clipping radius, in units of the standard deviation.

=item C<weight> => I<piddle>

I<Optional>. Data weights. When used with C<coords>, C<weight> must be
a piddle of shape I<M>, where I<M> is the number of data elements in
C<coords>. If C<coords> is not specified, C<weight> is a piddle of
shape I<NxM>, where I<N> is the number of dimensions in the data and
I<M> is the number of data elements.

=back

=head3 Sigma Clip Results

B<sigma_clip> returns an object which includes all of the attributes
from the final iteration object (See L</Sigma Clip Iterations> ), with
the following additional attributes/methods:

=over

=item C<iterations> => I<arrayref>

An array of results objects for each iteration.

=item C<success> => I<boolean>

True if the iteration converged, false otherwise.

=item C<error> => I<error object>

If convergence has failed, this will contain an error object
describing the failure.  See L</Errors>.

=back

=head4 Sigma Clip Iterations

The results for each iteration are stored in an object with the
following attributes/methods:

=over

=item C<center> => I<piddle|undef>

A 1D piddle containing the derived center.  The value for the last
iteration will be undefined if all of the elements have been clipped.

=item C<iter> => I<integer>

The iteration index.  An index of C<0> indicates the values determined
before the iterative loop was entered, and reflects the initial
clipping and mask exclusion.

=item C<nelem> => I<integer>

The number of data elements used in the center.

=item C<weight> => I<number>

The combined weight of the data elements used to determine the center.

=item C<sigma> => I<number|undef>

The standard deviation of the clipped data.  The value for the last
iteration will be undefined if all of the elements have been clipped.

=item C<clip> => I<number|undef>

The clipping radius.  This will be undefined for the first iteration
if the C<clip> option was not specified.

=item C<dist> => I<number>

I<Optional>. The distance between the previous and current centers. This is defined
only if the C<dtol> option was passed.

=back


=cut

use Hash::Wrap {
    -as     => 'new_iteration',
    -create => 1,
    -class  => 'PDL::Algorithm::Center::Iteration',
    -clone  => sub {
        my $hash = shift;

        return {
            map {
                my $value = $hash->{$_};
                $value = $value->copy if $value->$_isa( 'PDL' );
                ( $_, $value )
            } keys %$hash
        };
    },
  },
  {
    -as     => 'return_iterate_results',
    -class  => 'PDL::Algorithm::Center::Iterate::Results',
    -create => 1,
  };


sub sigma_clip {

    state $check = compile_named(
        center    => Optional [ ArrayRef[Num|Undef] | Center | CodeRef ],
        clip      => Optional [PositiveNum],
        coords    => Optional [Coords],
        dtol      => PositiveNum,
        iterlim   => Optional [PositiveInt],
        log       => Optional [CodeRef],
        mask      => Optional [Undef | Piddle_min1D_ne],
        nsigma    => PositiveNum,
        weight    => Optional [Undef | Piddle_min1D_ne],
    );

    my $opt = do {
        local $@;
        my %opt = eval { %{ $check->( @_ ) } };
        parameter_failure->throw( $@ ) if $@;
        wrap_hash( \%opt );
    };

    $opt->{iterlim} //= 10;

    #---------------------------------------------------------------

    # now, see what kind of data we have, and ensure that all dimensions
    # are consistent

    if ( defined $opt->{coords} ) {

        for my $name ( 'mask', 'weight' ) {

            my $value = $opt->{$name};
            next unless defined $value;

            parameter_failure->throw(
                "<$name> must be a 1D piddle if <coords> is specified" )
              if $value->ndims != 1;

            my $nelem_c = $value->getdim( -1 );
            my $nelem_p = $opt->coords->getdim( -1 );

            parameter_failure->throw(
                "number of elements in <$name> ($nelem_p) ) must be the same as in <coords> ($nelem_c)"
            ) if $nelem_c != $nelem_p;
        }

    }

    elsif ( defined $opt->{weight} ) {

        $opt->{coords} = $opt->weight->ndcoords( PDL::indx );

        if ( defined $opt->{mask} ) {
            parameter_failure->throw( "mask must have same shape as weight\n" )
              if  $opt->mask->shape != $opt->weight->shape;

            $opt->mask( $opt->mask->flat );
        }

        $opt->weight( $opt->weight->flat );
    }

    else {

        parameter_failure->throw( "must specify one of <coords> or <weight>" );
    }


    my ( $ndims ) = $opt->coords->dims;


    if ( defined $opt->{center} && is_arrayref( $opt->center) ) {

        my $icenter = pdl( @{ $opt->center } );

        parameter_failure->throw( "<center> must have $ndims elements" )
          unless  $icenter->nelem == $ndims ;

        my $defined = pdl( map { defined } @{ $opt->center } );

        if ( $defined->not->any ) {

            $icenter = $icenter->where( $defined )->sever;

            $opt->{center} = sub {

                my ( $coords, $wmask, $weight ) = @_;
                my $center = _weighted_mean_center( $coords, $wmask, $weight );
                $center->where( $defined) .= $icenter;;

                return $center;
            };
        }
    }
    else {

        $opt->{center} //= \&_weighted_mean_center;
    }


    my $nsigma = delete $opt->{nsigma};
    $opt->{calc_wmask} //= sub {
        _sigma_clip_wmask( $nsigma, @_ );
        return;
    };

    $opt->{calc_center} //= sub {
        my ( $coords, $wmask, $iter ) = @_;

        _weighted_mean_center( $coords, $wmask, $iter->weight );
    };

    my ( $clip, $dtol ) = delete @{$opt}{ 'clip', 'dtol' };
    $opt->{initialize} = sub {
        _sigma_clip_initialize( $clip, $dtol, @_ );
    };

    $opt->{is_converged} = sub {
        _sigma_clip_is_converged( $clip, $dtol, @_ );
    };

    delete @{$opt}{grep { ! defined $opt->{$_} } keys %$opt};


    iterate( %$opt );
}


=sub iterate

  $result = iterate(
    center       => Center | CodeRef,
    initialize   => CodeRef,
    calc_center  => CodeRef,
    calc_wmask   => CodeRef,
    is_converged => CodeRef,
    coords       => Coords,
    iterlim      => PositiveInt,
    log          => Optional [CodeRef],
    mask         => Optional [Piddle1D_ne],
    weight       => Optional [Piddle1D_ne],
  );

A generic iteration loop for centering data using callbacks for
calculating centers, weight masks, and iteration completion.

=head3 Options

The following options are accepted:

=over

=item C<center> => I<Piddle1D_ne | coderef >

The initial center.  It may either be a piddle with shape I<N> (or
something that can be coerced into one, see L</TYPES>) or a coderef
which will return the center as a piddle with shape I<N>.  The coderef
is called as

  $initial_center = &$center( $coords, $wmask, $weight );

with

=over

=item C<$coords>

a piddle with shape I<NxM> containing I<M> coordinates with dimension I<N>

=item C<$wmask>

a piddle of shape I<M> containing weights

=item C<$weight>

a scalar which is the sum of the weights in C<$wmask>

=back

=item C<initialize> => I<coderef>

This subroutine provides initialization prior to entering the
iteration loop.  It should initialize the passed iteration object and
work storage.

It is invoked as:

  &$initialize( $coords, $wmask, $current, $work );

with

=over

=item C<$coords>

a piddle of shape I<NxM> with the coordinates of each element

=item C<$wmask>

a piddle of shape I<M> with weights for each element

=item C<$current>

a reference to a L<Hash::Wrap> based object containing data for the
current iteration.  C<initialize> may augment the underlying hash with
its own data (but see L</Workspace>). The following attributes
are provided by C<iterate>:

=over

=item C<nelem>

the number of elements

=item C<weight>

the sum of the weights of the elements

=back

=item C<$work>

a hashref which  may use to store temporary data (e.g. work
piddles) which will be available to all of the callback routines.

=back

=item C<calc_center> => I<coderef>

This subroutine should return a piddle of shape I<N> with the
calculated center.

It will be called as:

  $center = &$calc_center( $coords, $wmask, $current, $work );

with

=over

=item C<$coords>

a piddle of shape I<NxM> with the coordinates of each element

=item C<$wmask>

a piddle of shape I<M> with weights for each element

=item C<$current>

a reference to a L<Hash::Wrap> based object containing
data for the current iteration.

C<calc_center> may augment the underlying hash with its own data (but
see L</Iteration Objects>). The following attributes are provided by
C<iterate>:

=over

=item C<nelem>

the number of elements

=item C<weight>

the sum of the weights of the elements

=back

=item C<$work>

a hashref which  may use to store temporary data (e.g. work
piddles) which will be available to all of the callback routines.

=back

=item C<calc_wmask> => I<coderef>

This subroutine should determine the current weight for each element.
To ignore an element set its weight to zero.

It will be called as:

  &$calc_mask( $coords, $wmask, $current, $work );

with

=over

=item C<$coords>

a piddle of shape I<NxM> with the coordinates of each element

=item C<$wmask>

a piddle of shape I<M> with the initial weights for each element( as
passed via the C<weight> option). C<calc_mask> should update it for
the current iteration.

=item C<$current>

a reference to a L<Hash::Wrap> based object containing data for the
current iteration.

C<calc_center> may augment the underlying hash with its own data (but
see L</Workspace>). The following attributes are provided by
C<iterate>:

=over

=item C<nelem>

the number of elements with non-zero weight.  If C<calc_mask> changes
C<$wmask>, this must either be updated or set to the undefined value.

=item C<weight>

the sum of the weights of the elements.  If C<calc_mask> changes
C<$wmask>, this must either be updated or set to the undefined value.

=back

=item C<$work>

a hashref which  may use to store temporary data (e.g. work
piddles) which will be available to all of the callback routines.

=back

=item C<is_converged> => I<coderef>

This subroutine should return a boolean value indicating whether the
iteration has converged.

It is invoked as:

  $bool = &$is_converged( $coords, $wmask, $last, $current, $work );

with

=over

=item C<$coords>

a piddle of shape I<NxM> with the coordinates of each element

=item C<$wmask>

a piddle of shape I<M> with weights for each element

=item C<$last>

a reference to a L<Hash::Wrap> based object containing data for the
previous iteration.  C<is_converged> may augment the underlying hash
with its own data (but see L</Workspace>). The following
attributes are provided by C<iterate>:

=over

=item C<nelem>

the number of elements

=item C<weight>

the sum of the weights of the elements

=back

=item C<$current>

a reference to a L<Hash::Wrap> based object containing data for the
current iteration, with attributes as described above for C<$last>

=item C<$work>

a hashref which  may use to store temporary data (e.g. work
piddles) which will be available to all of the callback routines.

=back

The C<is_converged> routine is passed references to the B<actual>
objects used by B<sigma_clip> to keep track of the iterations.  This
means that the C<is_converged> routine may manipulate the starting
point for the next iteration by altering its C<$current> parameter.

C<is_converged> is called prior to entering the iteration loop with
C<$last> set to C<undef>.  This allows priming the C<$current> structure,
which will be used as C<$last> in the first iteration.

=item C<coords> => I<Coords>

The coordinates to center.  C<coords> is a piddle of
shape I<NxM> (or anything which can be coerced into it, see
L</TYPES>) where I<N> is the number of dimensions in the data and
I<M> is the number of data elements.

=item C<iterlim>

a positive integer specifying the maximum number of iterations.

=item C<log> => I<coderef>

I<Optional>. A subroutine which will be called

=over

=item between the call to C<initialize> and the start of the first iteration

=item at the end of each iteration

=back

It is invoked as

  &$log( $iteration );

where C<$iteration> is a I<copy> of the current iteration object.  The object will
have at least the following fields:

=over

=item C<center> => I<piddle|undef>

A piddle of shape I<N> containing the derived center.  The value for
the last iteration will be undefined if all of the elements have been
clipped.

=item C<iter>

The iteration index

=item C<nelem>

The number of elements in the current selected set.

=item C<weight>

The summed weight of the selected elements.

=back

There may be other attributes added by the various callbacks
(C<calc_wmask>, C<calc_center>, C<is_converged>). See for example,
L</Sigma Clip Iterations>.

=item C<mask> => I<piddle>

I<Optional>. This specifies data elements to ignore completely.
True values indicate elements to be used, false those to be ignored.

When used with C<coords>, C<mask> must be a piddle of shape I<M>,
where I<M> is the number of data elements in C<coords>.

If C<coords> is not specified, C<mask> should have the same shape as
C<weight>.


=item C<weight> => I<piddle>

I<Optional>. Data weights.  When used with C<coords>, C<weight> must
be a piddle of shape I<M>, where I<M> is the number of data elements
in C<coords>. If C<coords> is not specified, C<weight> is a piddle of
shape I<NxM>, where I<N> is the number of dimensions in the data and
I<M> is the number of data elements.

=back

Callbacks are provided with L<Hash::Wrap> based objects which contain
the data for the current iteration.  They should add data to the
objects underlying hash which records particulars about their specific
operation,

=head3 Workspace

Callbacks are passed L<Hash::Wrap> based iteration objects and a
reference to a C<$work> hash.  The iteration objects may have additional
elements added to them (which will be available to the caller),
but should refrain from storing unnecessary data there, as each
new iteration's object is I<copied> from that for the previous iteration.

Instead, use the passed C<$work> hash.  It is shared amongst the
callbacks, so use it to store data which will not be returned to
the caller.

=head3 Results

B<iterate> returns an object which includes all of the attributes
from the final iteration object (See L</Iteration Object> ), with
the following additional attributes/methods:

=over

=item C<iterations> => I<arrayref>

An array of result objects for each iteration.

=item C<success> => I<boolean>

True if the iteration converged, false otherwise.

=item C<error> => I<error object>

If convergence has failed, this will contain an error object
describing the failure.  See L</Errors>.

=back

The value of the C<center> attribute in the last iteration will be
undefined if all of the elements have been clipped.

=head4 Iteration Object

The results for each iteration are stored in an object with the
following attributes/methods (in addition to those added by the
callbacks).

=over

=item C<center> => I<piddle|undef>

A 1D piddle containing the derived center.  The value for the last
iteration will be undefined if all of the elements have been clipped.

=item C<iter> => I<integer>

The iteration index.  An index of C<0> indicates the values determined
before the iterative loop was entered, and reflects the initial
clipping and mask exclusion.

=item C<nelem> => I<integer>

The number of data elements used in the center.

=item C<weight> => I<number>

The combined weight of the data elements used to determine the center.

=back

=head3 Iteration Steps

Before the first iteration:

=over

=item 1

Extract an initial center from C<center>.

=item 2

Create a new iteration object.

=item 3

Call C<initialize>.

=item 4

Call C<log>

=back

For each iteration:

=over

=item 1

Creat a new iteration object by B<copying> the old one.

=item 2

Call C<calc_wmask>

=item 3

Update summed weight and number of elements if C<calc_wmask> sets them to C<undef>.

=item 4

Call C<calc_center>

=item 5

Call C<is_converged>

=item 6

Call C<log>

=item 7

Goto step 1 if not converged and iteration limit has not been reached.

=back

=cut

sub iterate {

    state $check = compile_named(
        center       => Center | CodeRef,
        initialize   => CodeRef,
        calc_center  => CodeRef,
        calc_wmask   => CodeRef,
        is_converged => CodeRef,
        coords       => Coords,
        iterlim      => PositiveInt,
        log          => Optional [CodeRef],
        mask         => Optional [Piddle1D_ne],
        weight       => Optional [Piddle1D_ne],
    );

    my $opt = wrap_hash( $check->( @_ ) );

    $opt->{log} //= undef;

    my ( $ndims, $nelem ) = $opt->coords->dims;

    parameter_failure->throw( "<$_> must have $nelem elements" )
      for grep { defined $opt->{$_} && $opt->{$_}->nelem != $nelem }
      qw[ mask weight ];

    # use mask to remove unwanted elements
    if ( defined $opt->{mask} ) {

        my $select = $opt->mask != 0;
        $opt->coords( $opt->coords->mv(-1,0)->whereND( $select )->mv(0,-1)->sever );

        ( $ndims, $nelem ) = $opt->coords->dims;

        $opt->weight( $opt->weight->where( $select )->sever )
          if defined $opt->{weight};
    }

    my $wmask_base = defined $opt->{weight} ? PDL::convert( $opt->weight, PDL::double ) : PDL->ones( PDL::double, $opt->coords->getdim( -1 ) );
    my $wmask_base_weight = $wmask_base->dsum;
    my $wmask_base_nelem  = ( $wmask_base > 0 )->dsum;

    $opt->center( $opt->center->( $opt->coords, $wmask_base, $wmask_base_weight ) )
      if is_coderef( $opt->center );

    parameter_failure->throw( "<center> must be a 1D piddle with $ndims elements" )
      unless is_Piddle1D( $opt->center ) && $opt->center->nelem == $ndims;



    #############################################################################

    # Iterate until convergence.

    my @iteration;

    my $wmask = $wmask_base->copy;
    my $work  = {};

    # Set up initial state

    push @iteration,
      new_iteration( {
          center => $opt->center,
          weight => $wmask_base_weight,
          nelem => $wmask_base_nelem,
          iter  => 0,
      } );

    $opt->initialize->( $opt->coords, $wmask, $iteration[-1], $work );

    $opt->log && $opt->log->( new_iteration( $iteration[-1] ) );

    my $iteration;
    my $converged;

    eval {

        while ( ! $converged && ++$iteration <= $opt->iterlim ) {

            my $last = $iteration[-1];

            my $current = new_iteration( $last );
            push @iteration, $current;

            ++$current->{iter};

            $current->weight( $wmask_base_weight );
            $current->nelem( $wmask_base_nelem );
            $wmask .= $wmask_base;

            $opt->calc_wmask->( $opt->coords, $wmask, $current, $work );

            $current->weight( $wmask->dsum ) unless defined $current->weight;
            $current->nelem( ($wmask > 0)->dsum ) unless defined $current->nelem;

            iteration_empty_failure->throw( msg => "no elements left after clip" )
              if $current->nelem == 0;

            $current->center(
                $opt->calc_center->( $opt->coords, $wmask, $current, $work ) );

            $converged = $opt->is_converged->( $opt->coords, $wmask, $last, $current, $work );

            $opt->log && $opt->log->( new_iteration( $current ) );
        }

    };

    my $error = $@;

    $error
      = iteration_limit_reached_failure->new(
        msg => "iteration limit (@{[ $opt->iterlim ]}) reached" )
      if $iteration > $opt->iterlim;

    return_iterate_results( {
        %{ $iteration[-1] },
        iterations => \@iteration,
        success    => !$error,
        error      => $error
    } );
}




1;

# COPYRIGHT

__END__

=pod

=head1 DESCRIPTION

C<PDL::Algorithm::Center> is a collection of algorithms which
specialize in centering datasets.


=head1 SUBROUTINES

See L</TYPES> for information on the types used in the subroutine descriptions.


=head1 TYPES

In the L<description of the subroutines|/Subroutines>, the following
types are specified:

=over

=item Center

This accepts a non-null, non-empty 1D piddle, or anything that can be converted
into one (for example, a scalar, a scalar piddle, or an array of numbers );

=item CodeRef

A code reference.

=item PositiveNum

A positive real number.

=item PositiveInt

A positive integer.

=item Coords

This accepts a non-null, non-empty 2D piddle, or anything that can be converted or
up-converted to it.

=item Piddle_min1D_ne

This accepts a non-null, non-empty piddle with a minimum of 1 dimension.

=item Piddle1D_ne

This accepts a non-null, non-empty 1D piddle.

=back

=head1 ERRORS

Errors are represented as objects in the following classes:

=over

=item Parameter Validation

These are unconditionally thrown.

  PDL::Algorithm::Center::parameter

=item Iteration

These are stored in the result object's C<error> attribute.

  PDL::Algorithm::Center::iteration::limit_reached
  PDL::Algorithm::Center::iteration::empty

=back

The objects stringify to a failure message.

=head1 SEE ALSO
