package PDL::Algorithm::Center;

# ABSTRACT: Various methods of finding the center of a sample

use strict;
use warnings;

require 5.010000;

use feature 'state';
use Carp;

use Safe::Isa;
use Ref::Util qw< is_arrayref is_ref is_coderef  >;

use custom::failures;
use Package::Stash;
use Hash::Wrap;

use PDL::Algorithm::Center::Types -all;
use Types::Standard -types;
use Types::Common::Numeric -types;
use Type::Params qw[ compile_named ];

our $VERSION = '0.01';

use PDL::Lite ();

use Exporter 'import';

our @EXPORT_OK = qw[ sigma_clip iterate ];


my @failures;

BEGIN {

    my @failures = qw<
      parameter
      iteration::limit_reached
      iteration::empty
    >;

    custom::failures->import( __PACKAGE__, @failures );

    my $stash = Package::Stash->new( __PACKAGE__ );

    for my $failure ( @failures ) {

        ( my $name = $failure ) =~ s/::/_/g;
        $stash->add_symbol( "&${name}_error",
            sub (@) { ( __PACKAGE__ . "::$failure" )->throw( @_ ) } );
        $stash->add_symbol( "&${name}",
            sub (@) { ( __PACKAGE__ . "::$failure" )->new( @_ ) } );
    }
}

sub _error {

    my ( $type, $msg );

    if ( is_arrayref( $_[0] ) ) {

        $type = shift @{ $_[0] };
        $msg  = $_[0];
    }

    else {

        $type = shift;
        $msg  = \@_;

    }

    ( __PACKAGE__ . "::$type" )->throw( @$msg );

}


sub _weighted_mean_center {
    my ( $coords, $wmask, $weight ) = @_;

    $weight //= $wmask->dsum;

    return ( $coords * $wmask->dummy( 0 ) )->xchg( 0, 1 )->dsumover
      / $weight;
}

sub _sigma_clip_wmask {

    my ( $nsigma, $coords, $wmask, $iter, $work ) = @_;

    my $r2 = ( $work->{r2} //= PDL->null );
    $r2 .= ( ( $coords - $iter->center )**2 )->dsumover;

    $iter->clip( $nsigma * $iter->sigma );
    $wmask *= $r2 < $iter->clip**2;

    $iter->weight( $wmask->dsum );
    $iter->nelem( ( $wmask > 0 )->dsum );

    $iter->sigma( sqrt( ( $wmask * $r2 )->dsum / $iter->weight ) );

    return;
}


sub _distance {

    my ( $last, $current ) = @_;

    return sqrt( ( ( $last->center - $current->center )**2 )->dsum ) 
}

sub _sigma_clip_is_converged {

    my ( $init_clip, $dtol, $coords, $wmask, $last, $current, $work ) = @_;

    if ( ! defined $last ) {

        $current->{clip} = $init_clip;

        my $r2 = ( $work->{r2} //= PDL->null );
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

    $current->{dist} = undef;

    # stop if standard deviations and centers haven't changed

    if ( $current->sigma == $last->sigma
         && PDL::all( $current->center == $last->center) ) {
        $current->dist( _distance( $last, $current ) );
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


=sub iterate

=for ref

Center a dataset by iteratively excluding data outside of a radius equal to a specified number of standard deviations

Elements may be individually weighted.

=for usage

Usage:


  $results = sigma_clip( coords => $coords,
                         weight => $weight,
                         mask   => $mask,
                         %opts);

  $results = sigma_clip( weight => $img, %opts);


B<sigma_clip> finds the center of a data set by:

=over

=item 1

ignoring the data whose distance to the current center is greater than
a specified number of standard deviations

=item 2

calculating a new center by performing a (weighted) centroid of the
remaining data

=item 3

calculating the standard deviation of the distance from the data to
the center

=item 4

repeat at step 1 until either a convergence tolerance has been met or
the iteration limit has been exceeded

=back

The initial center may be explicitly specified, or may be calculated
by performing a (weighted) centroid of the data.

The initial standard deviation is calculated using the initial center and either
the entire dataset, or from a clipped region about the initial center.

=head3 Options

The following (case-insensitive) options are available:

=over 8

=item C<coords> => I<arrayref|piddle>

I<Optional> (see C<weight>)

The coordinates to center.  C<coords> may be either a I<N> element
list of piddles of shape I<M> or a single piddle of shape I<NxM>, where

=over

=item I<N> is the number of dimensions in the data

=item I<M> is the number of data elements

=back

C<coords> is useful if the data are sparse; for dense data, use C<weight> instead.

C<weight> may be specified with coords for sparse, weighted data.

=item C<weight> => I<piddle>

I<Optional> (see C<coords>)

Data weights.

For sparse data (i.e., when used with C<coords>)
C<weight> must be a piddle of shape I<M>, where I<M> is the number of
data elements in C<coords>.

For densely packed data, C<weight> is a piddle of shape I<NxM>, where

=over

=item I<N> is the number of dimensions in the data

=item I<M> is the number of data elements

=back

=item C<mask> = I<piddle>

I<Optional>

This specifies data elements to ignore completely.
True values indicate elements to be used, false those to be ignored.

For sparse data (i.e., when used with C<coords>) C<mask> must be a
piddle of shape I<M>, where I<M> is the number of data elements in
C<coords>.

For densely packed data, C<mask> should have the same shape as C<weight>.


=item C<center|centre> => I<arrayref>|I<1D piddle>|I<coderef>

I<Optional>

The initial center.  Defaults to the (weighted) average of the data.

=item C<is_converged> => I<subroutine reference>

I<Optional>

A subroutine which determines whether the iterations have converged.
It is called with two iteration objects which contain information about the
previous and current iterations, i.e.

    $stop_iteration = is_converged( $last, $current );

The structure of the objects is described in L</Iteration Results>.

It should return true if convergence has been achieved, false
otherwise.

The C<is_converged> routine is passed references to the B<actual>
objects used by B<sigma_clip> to keep track of the iterations.  This
means that the C<is_converged> routine may manipulate the starting
point for the next iteration by altering its C<$current> parameter.

The default behavior is to stop if both the standard deviation and
center have not changed between iterations, or if the C<dtol> option
was specified, the centers are closer than C<dtol>.

=item C<iterlim> => I<integer>

I<Optional>

The maximum number of iterations to run.  Defaults to 10.

=item C<dtol> => I<float>

I<Optional>

If specified, and the default convergence behavior is in use (see the
C<is_converged> option) iteration will cease when successive centers are
closer than the specified distance.

=item C<log> => I<subroutine reference>

I<Optional>

A subroutine which will be called at the end of each iteration. It is passed
a copy of the current iteration's results object see L</Iteration Results>).


=item C<transform> => PDL::Transform object

This is applied to the coordinates prior to centroiding.

=back

=head3 Iteration Results

The results for each iteration are stored in object of class
C<PDL::Algorithm::Center::sigma_clip::Iteration> with the
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

=item C<weight> => I<float>

The combined weight of the data elements used to determine the center.

=item C<sigma> => I<float|undef>

The standard deviation of the data.  The value for the last iteration
will be undefined if all of the elements have been clipped.

=item C<variance> => I<float|undef>

The calculated variance (i.e., C<sqrt( sigma )>.  The value for the
last iteration will be undefined if all of the elements have been
clipped.

=item C<clip> => I<float|undef>

The clipping radius.  This will be undefined for the first iteration
if the C<clip> option was not specified.

=item C<dist> => I<float>

I<Optional>

The distance between the previous and current centers. This is present
only if the default convergence routine is in use.

=back

=head3 Returned Results

B<sigma_clip> returns an object of class
C<PDL::Algorithm::Center::sigma_clip::Result>.  It is a subclass of
C<PDL::Algorithm::Center::sigma_clip::Iteration> (the common
attributes refer to the results of the final iteration) with these
additional attributes/methods:

=over

=item C<iterations> => I<arrayref>

A list of the iteration result objects.

=item C<success> => I<boolean>

True if the iteration converged, false otherwise.

=item C<error> => I<error object>

If convergence has failed, this will contain an error object
describing the failure.  See L</Errors>.

=back


=head3 Errors

Errors are represented as objects in the following classes:

=over

=item Parameter Validation

These are unconditionally thrown.

  PDL::Algorithm::Center::parameter
  PDL::Algorithm::Center::parameter::type
  PDL::Algorithm::Center::parameter::dimension
  PDL::Algorithm::Center::parameter::missing
  PDL::Algorithm::Center::parameter::value

=item Iteration

These are stored in the result object's C<error> attribute.

  PDL::Algorithm::Center::iteration::limit_reached
  PDL::Algorithm::Center::iteration::empty

=back

The objects stringify to a failure message.

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
    -as     => 'return_iterate',
    -class  => 'PDL::Algorithm::Center::Iterate',
    -create => 1,
  };


sub _upgrade_dims {

    return
         $_[0]->$_isa( 'PDL' )
      && !$_[0]->isempty
      && $_[0]->ndims == 0 ? $_[0]->dummy( 0 ) : $_[0];
}

sub sigma_clip {

    state $check = compile_named(
        center    => Optional [ Center | CodeRef ],
        clip      => Optional [PositiveNum],
        converged => Optional [CodeRef],
        coords    => Optional [Coords],
        dtol      => PositiveNum,
        iterlim   => Optional [PositiveInt],
        log       => Optional [CodeRef],
        mask      => Optional [Piddle_min1D_ne],
        nsigma    => PositiveNum,
        transform => Optional [ InstanceOf ['PDL::Transform'] ],
        weight    => Optional [Piddle_min1D_ne],
    );

    my %opt;
    do {
        local $@;
        %opt = eval { %{ $check->( @_ ) } };
        parameter_error( $@ ) if $@;
    };

    $opt{iterlim} //= 10;

    my $opt = wrap_hash( \%opt );

    #---------------------------------------------------------------

    # now, see what kind of data we have, and ensure that all dimensions
    # are consistent

    if ( defined $opt{coords} ) {

        for my $name ( 'mask', 'weight' ) {

            my $value = $opt{$name};
            next unless defined $value;

            parameter_error(
                "<$name> must be a 1D piddle if <coords> is specified" )
              if $value->ndims != 1;

            my $nelem_c = $value->getdim( -1 );
            my $nelem_p = $opt->coords->getdim( -1 );

            parameter_error(
                "number of elements in <$name> ($nelem_p) ) must be the same as in <coords> ($nelem_c)"
            ) if $nelem_c != $nelem_p;
        }

    }

    elsif ( defined $opt{weight} ) {

        $opt{coords} = $opt->weight->ndcoords( PDL::indx );

        parameter_error( "mask must have same shape as weight\n" )
          if defined $opt{mask}
          && $opt->mask->shape != $opt->weight->shape;


        $opt->weight( $opt->weight->flat );
    }

    else {

        parameter_error( "must specify one of <coords> or <weight>" );
    }

    $opt->coords( $opt->transform->apply( $opt->coords ) ) if $opt{transform};

    $opt{center} //= \&_weighted_mean_center;
    my $nsigma = delete $opt{nsigma};
    $opt{calc_wmask} //= sub {
        _sigma_clip_wmask( $nsigma, @_ );
        return;
    };

    $opt{calc_center} //= sub {
        my ( $coords, $wmask, $iter, $work ) = @_;

        _weighted_mean_center( $coords, $wmask, $iter->weight );
    };

    my ( $clip, $dtol ) = delete @{opt}{ 'clip', 'dtol' };
    $opt{is_converged} //= sub {
        _sigma_clip_is_converged( $clip, $dtol, @_ );
    };


    $DB::single=1;

    iterate( %opt );

}


sub iterate {

    state $check = compile_named(
        calc_center  => CodeRef,
        calc_wmask   => CodeRef,
        center       => Center | CodeRef,
        is_converged => CodeRef,
        coords       => Coords,
        iterlim      => PositiveInt,
        log          => Optional [CodeRef],
        mask         => Optional [Piddle1D_ne],
        weight       => Optional [Piddle1D_ne],
    );

    my %opt = %{ $check->( @_ ) };
    my $opt = wrap_hash( \%opt );

    $DB::single=1;

    my ( $ndims, $nelem ) = $opt->coords->dims;

    # mask for inclusion of data elements
    #
    # * Use the non-zero elements in the product of the inputmask and
    #   weight as a selection mask.
    # * Make this a double, as it will be multipled by $weight again.
    # * It will be 1D, to match the coords.


    my $wmask_base = do {

        parameter_error( "<$_> must have $nelem elements" )
          for grep { defined $opt{$_} && $opt{$_}->nelem != $nelem }
          qw[ mask weight ];

        # weighted mask. not a boolean mask
        my $wmask
          = defined $opt{weight} && defined $opt{mask}
          ? $opt{weight} * $opt{mask}
          : $opt{weight} // $opt{mask};

        defined $wmask
          ? PDL::convert( $wmask, PDL::double )
          : PDL->ones( PDL::double, $opt{coords}->getdim( -1 ) );
    };

    my $wmask_base_weight = $wmask_base->dsum;
    my $wmask_base_nelem  = ( $wmask_base > 0 )->dsum;


    $opt->center( $opt->center->( $opt->coords, $wmask_base, $wmask_base_weight ) )
      if is_coderef( $opt->center );

    parameter_error( "<center> must be a 1D piddle with $ndims elements" )
      unless is_Piddle1D( $opt->center ) && $opt->center->nelem == $ndims;



    #############################################################################

    # Iterate until convergence.

    my @iteration;

    # storage to avoid more memory thrashing. will get allocated upon
    # first use
    my $wmask = $wmask_base->copy;
    my $work  = {};

    # Set up initial state

    push @iteration,
      new_iteration( {
          center => $opt->center,
          weight => $wmask_base_weight,
          nelem => $wmask_base_nelem,
      } );

    $opt{is_converged}->( $opt->coords, $wmask, undef, $iteration[-1], $work );

    defined $opt{log} && $opt->log->( new_iteration( $iteration[-1] ) );

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

            iteration_empty_error( msg => "no elements left after clip" )
              if $current->nelem == 0;

            $current->center(
                $opt->calc_center->( $opt->coords, $wmask, $current, $work ) );

            $converged = $opt->is_converged->( $opt->coords, $wmask, $last, $current, $work );

            defined $opt{log} && $opt->log->( new_iteration( $current ) );
        }

    };

    my $error = $@;

    $error
      = iteration_limit_reached(
        msg => "iteration limit (@{[ $opt->iterlim ]}) reached" )
      if $iteration > $opt->iterlim;

    return_iterate( {
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

=head1 SYNOPSIS


=head1 SEE ALSO

