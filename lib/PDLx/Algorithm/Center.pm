# --8<--8<--8<--8<--
#
# Copyright (C) 2011 Smithsonian Astrophysical Observatory
#
# This file is part of PDLx::Algorithm::Center
#
# PDLx::Algorithm::Center is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -->8-->8-->8-->8--

package PDLx::Algorithm::Center;

use strict;
use warnings;

require 5.010000;

use feature 'state';
use Carp;

use Safe::Isa;
use Ref::Util qw< is_arrayref is_ref  >;

use custom::failures;
use Package::Stash;

use Validate::Tiny qw< validate is_required >;
use PDLx::Algorithm::Center::Validate;

our $VERSION = '0.01';

use Exporter 'import';

our @EXPORT_OK = qw[ sigma_clip ];

use PDL;
use PDL::Options;

my @failures;

BEGIN {

    my @failures = qw< parameter
		       parameter::type
		       parameter::dimension
		       parameter::missing
		       parameter::value
		       iteration::limit_reached
		       iteration::empty
		     >;

    custom::failures->import( __PACKAGE__, @failures );

    my $stash = Package::Stash->new( __PACKAGE__ );

    for my $failure ( @failures ) {

	( my $name = $failure ) =~ s/::/_/g;
	$stash->add_symbol( "&${name}_error", sub (@) { (__PACKAGE__ . "::$failure")->throw( @_ ) } );
	$stash->add_symbol( "&${name}", sub (@) { (__PACKAGE__ . "::$failure")->new( @_ ) } );
    }
}

sub _error {

    my ( $type, $msg );

    if ( is_arrayref( $_[0] ) ) {

	$type = shift @{$_[0]};
	$msg = $_[0];
    }

    else {

	$type = shift;
	$msg = \@_;

    }

    ( __PACKAGE__ . "::$type")->throw( @$msg );

}


## no critic (ProhibitAccessOfPrivateData)

=pod


=head1 PDLx::Algorithm::Center

Various methods of finding the center of a sample

=head1 FUNCTIONS

=head2 sigma_clip

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

ignoring the data whose distance to the current center is a specified
number of standard deviations

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


=item C<nsigma> => I<float>

I<Required>

The number of standard deviations at which to clip.

=item clip => I<scalar|arrayref>

I<Optional>

The clipping radius used to determine the initial standard deviation.

=item C<center|centre> => I<arrayref>|I<1D piddle>

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


=back

=head3 Iteration Results

The results for each iteration are stored in object of class
C<PDLx::Algorithm::Center::sigma_clip::Iteration> with the
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

The standard deviation of the data.
The value for the last iteration will be
undefined if all of the elements have been clipped.

=item C<variance> => I<float|undef>

The calculated variance (i.e., C<sqrt( sigma )>.
The value for the last iteration will be
undefined if all of the elements have been clipped.

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
C<PDLx::Algorithm::Center::sigma_clip::Result>.  It is a subclass
of C<PDLx::Algorithm::Center::sigma_clip::Iteration> (the common
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

  PDLx::Algorithm::Center::parameter
  PDLx::Algorithm::Center::parameter::type
  PDLx::Algorithm::Center::parameter::dimension
  PDLx::Algorithm::Center::parameter::missing
  PDLx::Algorithm::Center::parameter::value

=item Iteration

These are stored in the result object's C<error> attribute.

  PDLx::Algorithm::Center::iteration::limit_reached
  PDLx::Algorithm::Center::iteration::empty

=back

The objects stringify to a failure message.

=cut

{
    package PDLx::Algorithm::Center::sigma_clip::Iteration;

    use Safe::Isa;

    use Class::Tiny qw( center
      iter
      nelem
      weight
      sigma
      variance
      clip
      dist
    );

    sub copy {

        my $self = shift;

        return __PACKAGE__->new(
            map {
                my $key   = $_;
                my $value = $self->{$_};
                $value = $value->copy if $value->$_isa( 'PDL' );
                ( $key, $value )
            } Class::Tiny->get_all_attributes_for( __PACKAGE__ ) );
    }
}

{
    package PDLx::Algorithm::Center::sigma_clip::Result; 

    use parent -norequire, 'PDLx::Algorithm::Center::sigma_clip::Iteration';

    use Class::Tiny qw(
			iterations
			success
			error
    );

}

sub sigma_clip
{
    my %uopts = @_;


    #---------------------------------------------------------------

    # Options, first
    my $opts = PDL::Options->new (
		      {
		       nsigma 	     => undef,
		       clip   	     => undef,
		       iterlim       => 10,
		       center 	     => undef,
		       dtol   	     => undef,
		       mask   	     => undef,
		       log 	     => undef,
		       is_converged     => undef,
		       coords 	     => undef,
		       weight 	     => undef,
		       transform     => undef,
		      } );
    $opts->add_synonym( { 'CENTRE' => 'CENTER' } );

    my $opt = $opts->options( \%uopts );

    # vaildate the easy stuff
    state $spec = {

        fields => [],

        filters => [

            # convert from an array of piddles into a
            # single piddle of coordinate tuples
            coords => sub {

                if (
                    is_arrayref( $_[0] )
                    && !grep {
                        !(     $_->$_isa( 'PDL' )
                            && $_->ndims <= 1
                            && $_->nelem == $_[0][0]->nelem )
                    } @{ $_[0] } )
                {

                    return PDL->glue( map { $_->squeeze->dummy( 0 ) }
                          @{ $_[0] } );
                }


                # make sure there's a real dimension
                return $_[0]->dummy( 0 )
                  if $_[0]->$_isa( 'PDL' )
                  && !$_[0]->isempty
                  && $_[0]->ndims == 0;


                return $_[0];
            },


            center => sub {

                return if !defined $_[0];

                # try and convert. if it doesn't it'll get picked up
                # late in the check suite
                eval {
                    my $pdl = PDL::Core::topdl( $_[0] );
                    return $pdl if $pdl->isempty;

                    # make sure there's a real dimension
                    return $pdl->ndims == 0 ? $pdl->dummy( 0 ) : $pdl;
                };

                return $_[0];
            },

            # upgrade a non-empty piddle of 0 dimensions to 1 dimension
            [ 'mask', 'weight' ] => sub {
                return $_[0]->dummy( 0 )
                  if $_[0]->$_isa( 'PDL' )
                  && !$_[0]->isempty
                  && $_[0]->ndims == 0;

                $_[0];
            }
          ],

        checks => [

            [ 'clip', 'dtol' ] => is_a_positive_number,

	    nsigma => [ is_required, is_a_positive_number ],

            iterlim => [ is_required, is_a_positive_integer ],

	    # accept point PDL (ndims == 0 ), e.g. pdl(1)
            center => is_a_ndim_PDL( 0, 1, "must be a 1D piddle or an array of numbers" ),

            transform => is_a( 'PDL::Transform' ),

            # more is done later if these are both present or present
            # with coords
            [ 'mask', 'weight' ] => is_a_nonempty_PDL,

            # coords might have been converted from an array, so
            # need to emit an error message which encompasses that
            # possibility
	    # accept point PDL (ndims == 0 ), e.g. pdl(1)
            coords => is_a_ndim_PDL(
                0,
                2,
                "must be an array of 1D piddles or a non-empty 1D or 2D piddle"
            ),

            [ 'log', 'is_converged' ] => is_a_subroutine_reference,

        ],
    };

    my $result = validate($opt, $spec);


    # grab the first error that we get and throw on it
    _error( ( values %{ $result->{error} })[0] )
      unless $result->{success};


    # this is redundant at the moment.
    $opt = $result->{data};

    #---------------------------------------------------------------

    # now, see what kind of data we have, and ensure that all dimensions
    # are consistent

    my $ndim;			# dimensions of the data


    # This starts out as a multiplicative combination of the input
    # mask and the input weight.  Later, elements which are excluded
    # from the calculation are assigned a mask of zero.


    my $coords = $opt->{coords};

    if ( defined $coords ) {

	$ndim = $coords->dims(0);

	for my $name ( 'mask', 'weight' ) {

	    my $value = $opt->{$name};
	    next unless defined $value;

	    parameter_dimension_error( "<$name> must be a 1D piddle if <coords> is specified" )
	      if $value->ndims != 1;

	    my $nelem_c = $value->dim( -1 );
	    my $nelem_p = $coords->dim( -1 );

	    parameter_dimension_error( "number of elements in <$name> ($nelem_p) ) must be the same as in <coords> ($nelem_c)" )
	      if $nelem_c != $nelem_p;
	}

    }

    elsif ( defined $opt->{weight} ) {

	$coords = $opt->{weight}->ndcoords( indx );

	parameter_dimension_error( "mask must have same shape as weight\n" )
	  if defined $opt->{mask} && $opt->{mask}->shape != $opt->{weight}->shape;

	$ndim = $opt->{weight}->ndims;

    }

    else {

	parameter_missing_error( "must specify one of <coords> or <weight>" )
    }


    $coords = $opt->{transform}->apply( $coords ) if $opt->{transform};

    # mask for inclusion of data elements
    #
    # * Use the non-zero elements in the product of the inputmask and
    #   weight as a selection mask.
    # * Make this a double, as it will be multipled by $weight again.
    # * It will be 1D, to match the coords.


    my $mask_base = do {

        # weighted mask. not a boolean mask
        my $wmask
          = defined $opt->{weight} && defined $opt->{mask}
          ? $opt->{weight} * $opt->{mask}
          : $opt->{weight} // $opt->{mask};

        defined $wmask
          ? convert( $wmask != 0, double )
          : ones( double, $coords->dim( -1 ) );
      };

    parameter_dimension_error( "<center> must have $ndim elements" )
      if defined $opt->{center} && $opt->{center}->nelem != $ndim;

    my $nvar = $opt->{nsigma} ** 2;

    unless( defined $opt->{is_converged} ) {

	$opt->{is_converged} = sub {

	    my ( $last, $current ) = @_;

	    # stop if standard deviations and centers haven't changed
	    return 1 if $current->{sigma} == $last->{sigma}
	      && all( $current->{center} == $last->{center} );

	    # or, if a tolerance was defined, stop if distance from old
	    # to new centers is less than the tolerance.
	    if ( defined $opt->{dtol} )
	    {
		$current->{dist} = sqrt( dsum( ( $last->{center} - $current->{center} )**2 ) );
		return 1 if $current->{dist} <= $opt->{dtol};
	    }

	    return 0;
	}
    }


    #############################################################################

    # Iterate until convergence.

    my @iteration;

    # storage to avoid more memory thrashing. will get allocated upon
    # first use
    my $wmask = PDL->null;
    my $r2 = PDL->null;
    my $clipped = PDL->null;

    # Set up initial state

    push @iteration,
      _sigma_clip_iter0( $mask_base, $coords, $opt->{center}, $opt->{weight}, $opt->{clip}, $wmask, $clipped );

    $opt->{log} && $opt->{log}->( $iteration[-1]->copy );

    my $done = 0;
    my $error;
    while ( $opt->{iterlim}-- && !$done ) {

        my $last = $iteration[-1];

        # determine distance**2 to current center
        $r2 .= dsumover( ( $coords - $last->{center} )**2 );

        # derive new clipping radius based upon variance within old radius
        my $clip2 = $nvar * $last->variance;

        # reset mask and mask out those too far out
        $clipped .= $r2 < $clip2;
        $wmask   .= $mask_base * $clipped;
        my $nelem = my $total_weight = $wmask->sum;

        if ( defined $opt->{weight} ) {
            $wmask *= $opt->{weight};
            $total_weight = $wmask->dsum;
        }

        my ( $center, $variance ) = do {

            if ( $total_weight == 0 ) {

		$error = iteration_empty( msg => "no elements left after clip" );
                ( undef, undef );
            }

            else {

                (
                    _centroid( $coords, $wmask, $total_weight ),
                    ( $wmask * $r2 )->dsum / $total_weight
                );

            }
        };

        push @iteration,
          PDLx::Algorithm::Center::sigma_clip::Iteration->new(
            center   => $center,
            iter     => @iteration + 0,
            nelem    => $nelem,
            weight   => $total_weight,
            sigma    => defined $variance ? sqrt( $variance ) : undef,
            variance => $variance,
            clip     => sqrt( $clip2 ),
          );

	# do we stop because things can't continue?
	my $valid = defined $variance && defined $center;

	$done = $opt->{is_converged}->( $iteration[-2], $iteration[-1] )
	  if $valid;

        $opt->{log} && $opt->{log}->( $iteration[-1]->copy );

        last if ! $valid;
    }

    $error = iteration_limit_reached( msg => "iteration limit ($opt->{iterlim}) reached" )
      if ! $done && ! $error;

    return PDLx::Algorithm::Center::sigma_clip::Result->new(
        %{ $iteration[-1] },
        iterations => \@iteration,
        success    => $done,
        error      => $error
    );

}

sub _sigma_clip_iter0 {

    my ( $mask_base, $coords, $center, $weight, $clip, $wmask_storage, $clipped ) = @_;

    my $nelem = $mask_base->sum;
    my $wmask;
    my $total_weight;

    # get a weighted centroid if no center is specified
    if ( ! defined $center ) {

	if ( defined $weight ) {

	    $wmask = $wmask_storage;

	    $wmask .= $mask_base * $weight;
	    $total_weight = $wmask->dsum;
	}

	else {

	    $wmask = $mask_base;
	    $total_weight = $nelem;
	}

	$center = _centroid( $coords, $wmask, $total_weight );
    }

    # initial distances
    my $r2 = dsumover( ( $coords - $center)**2 );


    # initial clip radius. this defines $clip2, which will get
    # recorded in the log
    my $clip2;

    if ( defined $clip ) {

	$wmask = $wmask_storage;

	$clip2 = $clip ** 2;

	$clipped .= $r2 < $clip2;
	$wmask .= $mask_base & $clipped;

	$nelem = $total_weight = $wmask->sum;

	if ( defined $weight ) {
	    $wmask *= $weight;
	    $total_weight = $wmask->dsum;
	}
    }
    else {

	if ( defined $weight ) {

	    $wmask = $wmask_storage;

	    $wmask .= $mask_base * $weight;
	    $total_weight = $wmask->dsum;

	}

	else {

	    $wmask = $mask_base;
	    $total_weight = $nelem;

	}
    }

    # now the variance, which is wt*dist^2 / sum(wt);
    my $variance = ($wmask * $r2 )->dsum / $total_weight;

    return PDLx::Algorithm::Center::sigma_clip::Iteration->new( 
        center => $center->copy,
        iter   =>  0,
        nelem  => $nelem,
        weight => $total_weight,
        variance => $variance,
        sigma  => sqrt( $variance ),
        clip => defined $clip2 ? sqrt( $clip2 ) : undef,
    );

}


sub _centroid
{
  my ( $coords, $weight, $total_weight ) = @_;

  return ( $coords * $weight->dummy(0) )->xchg(0,1)->dsumover / $total_weight ;
}

1;

__END__

=pod

=head1 AUTHOR

Diab Jerius, E<lt>djerius@cpan.orgE<gt>


=head1 COPYRIGHT AND LICENSE

Copyright 2016 Smithsonian Astrophysical Observatory

This software is released under the GNU General Public License.  You
may find a copy at

          http://www.gnu.org/licenses



=cut
