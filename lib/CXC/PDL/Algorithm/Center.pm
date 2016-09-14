# --8<--8<--8<--8<--
#
# Copyright (C) 2011 Smithsonian Astrophysical Observatory
#
# This file is part of CXC::PDL::Algorithm::Center
#
# CXC::PDL::Algorithm::Center is free software: you can redistribute
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

package CXC::PDL::Algorithm::Center;

use strict;
use warnings;
use Carp;

our $VERSION = '0.01';

use parent 'Exporter';

our @EXPORT_OK = qw[ maxpix sigma_clip ];

use PDL;
use PDL::Options;

## no critic (ProhibitAccessOfPrivateData)

sub _sclr_or_arr {

    my ( $param, $value, $ndims ) = @_;

    my @value;

    if ( ! ref $value )
    {
	@value = ( $value ) x $ndims;
    }

    elsif ( 'ARRAY' eq ref $value )
    {
	barf( 'array dimensions ( ',
	      scalar @$value,
	      " ) must be the same as the input data ( $ndims )\n" )
	  unless $ndims == @$value;

	@value = @$value;
    }
    else
    {
	barf( "parameter $param must be scalar or array of $ndims dimensions\n" );
    }

    return @value;
}


# ============================================================

=pod

=head2 maxpix

=for ref

determine coordinates of maximum pixel of binned data

The input piddles are treated as sets coordinates and are binned into
an image. The coordinates of the pixel with the largest value are
returned.  The data may have weights.  It is currently limited
to two-dimensional data sets.

=for usage

@coords = maxpix( @coords, \%options )

=for options

=over

=item shape

The shape of the image, in pixels.  If it is a scalar it is used for
all dimensions.  Otherwise it should be a Perl array with the extents
of each dimension.

=item center

The center of the image, in input coordinates.  If it is a scalar it
is used for all dimensions.  Otherwise it should be a Perl array with
the center of each dimension.

=item weight

If present, this is a piddle of weights for each coordinate tuple.

=back

=cut

sub maxpix
{
    my %uopts = @_;

    my $opts = PDL::Options->new( { shape   => 512,   # output dimensions
				    center => undef,  # initial center
				    coords => undef,  # coordinate data
				    weight => undef,  # data weight
				    shrink => 0.5,    # shrinkage factor
				    minmax => undef,  # optional, minimum value for max pixel
				  }

				);

    my $opt = $opts->options( \%uopts );

    barf( "must specify one of coords or weight\n" )
      unless defined $opt->{coords} or defined $opt->{weight};

    my $data;
    my $ndim;

    if ( defined $opt->{coords} )
    {
	$data = $opt->{coords};

	# if we're passed an array, there is one 1D piddle per coordinate
	if ( 'ARRAY' eq ref $data )
	{
	    # these are data indices; must be 1D
	    barf( "coordinates must be 1D!\n" )
	      if grep { $_->dims > 1 } @$data;

	    barf( "inconsistent extent of indices\n" )
	      if grep { $_->nelem != $data->[0]->nelem } @$data;

	    barf( "weight must have same shape as data\n" )
	      if defined $opt->{weight} &&
		(    $opt->{weight}->dims != 1
		  || $opt->{weight}->nelem != $data->[0]->nelem );

	    $ndim = @$data;
	}

	# eventually this will be a piddle of coordinate tuples
	else
	{
	    barf( "single piddle of coordinates is not yet supported\n" );
	}
    }

    else {

	barf( "dense data is not yet supported\n" );

    }


    barf( "cannot handle input data other than two dimensions\n" )
      if $ndim != 2;

    my $shape = pdl( _sclr_or_arr( 'shape', $opt->{shape}, $ndim ) );

    my @center;
    @center = _sclr_or_arr( 'center', $opt->{center}, $ndim )
      if defined $opt->{center};

    my $shrink = pdl( _sclr_or_arr( 'shrink', $opt->{shrink}, $ndim ) );

    require Img2D::Events;
    require PDL::Image2D;
    require PDL::Transform;

    my ( $image, $max, $ix, $iy );

    while( 1 ) {

	my %opts;
	@opts{'xpix','ypix'} = $shape->list;

	@opts{'xc','yc'} = @center
	  if @center;

	$image = Img2D::Events->new( @$data, $opts->{weight} ? $opts->{weight} : (), \%opts )->img;

	$image->hdr->{NAXIS} = 2;

	( $max, $ix, $iy ) = $image->max2d_ind;

	last if ! defined $opt->{minmax} || $max > $opt->{minmax};

	$shape = ( $shape * $shrink)->floor;

	barf( __PACKAGE__ . "::maxpix: image shrunk too far\n" )
	  if $shape->max <= 1;
    }

    my $xfrm = PDL::Transform::t_fits( $image );

    my $idx = double( cat( $ix, $iy ) );

    my %results = ( image  => $image,
		    pixpos => $idx,
		    max    => $max,
		    trans  => $xfrm,
		    center => [$xfrm->apply( $idx )->list],
		  );

    return wantarray ? @{ $results{center} } : \%results;
}

# ============================================================

=pod

=head2 sigma_clip

=for ref

Center a dataset by performing interative sigma clipping.

Elements may be individually weighted.

It returns the center, the standard deviation, the number of iterations,
the number of remaining elements, the final clipping distance, and the
specified number of sigma at which to clip.

=for usage

Usage:


  # @idx is a list of 1D piddles of coordinate values
  @center = sigma_clip( coords => \@idx, %opts);
  @center = sigma_clip( coords => \@idx, weight => $wt, %opts);

  # $img is an nD piddle of weights
  @center = sigma_clip( weight => $img, %opts);

  # $coords is an NxM piddle, with M tuples of N-dimensional coordinates
  @center = sigma_clip( coords => $coords, %opts);
  @center = sigma_clip( coords => $coords, %opts);

  # results are returned as a list, a hashref, or as a hash
  ( $center, $sigma, $iter, $n, $r, $nsigma, [$mask] ) = sigma_clip( ... );

  $results = sigma_clip( ... );
  $results = { center => $center, sigma => $sigma, iter => $iter,
	       n => $n, r => $r, nsigma => $nsigma };

  %results = ( center => $center, sigma => $sigma, iter => $iter,
	       n => $n, r => $r, nsigma => $nsigma, [ mask => $mask ] );


If called in array mode, the indicated values are returned as a list.
If called in scalar mode, a reference to a
hash with the keys C<center>, C<sigma>, C<iter>, C<n>, C<r>, C<nsigma> is
returned.

The values returned are:

=over 8

=item center

The derived center.  This is an arrayref, one element per coordinate axis.

=item sigma

The standard deviation of the data surviving the clip in the last iteration.

=item iter

The number of iterations required to converge.

=item n

The number of elements remaining in the last iteration.

=item r

The final clipping radius.

=item nsigma

The number of standard deviations used in clipping.

=back

The following (case-insensitive) options are available:

=over 8

=item NSigma

The number of sigma at which to clip.  Defaults to 3.

=item Clip

The initial clipping radius. This defaults to infinity.

=item Center

The initial center.  This is an arrayref containing a first
guess at the center of the data.  If the data are one dimensional,
this may be a scalar.

If not specified, the (weighted) average of
the data is used.

=item Iterlim

The maximum number of iterations to run.  Defaults to 10.

=item dtol

Successive centers must be no farther apart than this for the centering
to succeed.  If this is not defined, then centering will stop
when the centroids and variances of two successive iterations are
exactly the same (or if the iteration limit has been hit).

=item Mask

This may be set to a piddle which indicates elements in the data
to be ignored.  The mask should be set to true for those
elements to use, false for those to be ignored.

=back

=cut

sub sigma_clip
{
    my %uopts = @_;

    #---------------------------------------------------------------

    # Options, first
    my $opts = PDL::Options->new (
		      {
		       nsigma => 3,
		       clip   => undef,
		       iterlim  => 10,
		       center => undef,
		       dtol   => undef,
		       mask   => undef,
		       log => sub { 1 },
		       converged => undef,
		       coords => undef,
		       weight => undef,
		      } );

    $opts->add_synonym( { 'CENTRE' => 'CENTER' } );

    my $opt = $opts->options( \%uopts );

    #---------------------------------------------------------------

    # now, see what kind of data we have, and ensure that all dimensions
    # are consistent

    my $ndim;			# dimensions of the data

    barf( "must specify one of coords or weight\n" )
      unless defined $opt->{coords} or defined $opt->{weight};

    # sparse data means coords or coords + weight
    # dense  data means weight only
    my $is_sparse;

    my $data;

    if ( defined $opt->{coords} )
    {
	$is_sparse = 1;
	$data = $opt->{coords};

	# if we're passed an array, there is one 1D piddles per coordinate
	if ( 'ARRAY' eq ref $data )
	{
	    # these are data indices; must be 1D
	    barf( "coordinates must be 1D!\n" )
	      if grep { $_->dims > 1 } @$data;

	    barf( "inconsistent extent of indices\n" )
	      if grep { $_->nelem != $data->[0]->nelem } @$data;

	    barf( "weight must have same shape as data\n" )
	      if defined $opt->{weight} &&
		(    $opt->{weight}->dims != 1
		  || $opt->{weight}->nelem != $data->[0]->nelem );

	    $ndim = @$data;
	}

	# eventually this will be a piddle of coordinate tuples
	else
	{
	    barf( "single piddle of coordinates is not yet supported\n" );
	}
    }

    else {

	$is_sparse = 0;
	$ndim = $opt->{weight}->dims;
	$opt->{coords} = $opt->{weight};
    }

    my $wt = $opt->{weight};

    # was a center defined?  if so, ensure that it has the correct dimensions
    my @center;
    if ( defined $opt->{center} )
    {
	if ( 'ARRAY' eq ref $opt->{center} )
	{
	    @center = @{$opt->{center}};
	}

	# we'll accept a scalar as the center for a 1D data set
	else
	{
	    @center = ( $opt->{center} );
	}

	barf( 'dimensionality of center (' . scalar @center .
	      ") inconsistent with that of data ($ndim)\n")
	  if  @center != $ndim;
    }

    # no center provided; make room for it.
    else
    {
	@center = (undef) x $ndim;
    }

    my $nvar = $opt->{nsigma} ** 2;

    # initial clip radius
    my $clip2 = $opt->{clip} ? $opt->{clip} ** 2 : undef;

    # mask for inclusion of data elements; this may be multiplied
    # by a floating point weight, so make it a double.
    my $mask = ones( double, $is_sparse ? $data->[0]->dims : $data->dims );

    # storage for squared distances
    my $r2 = PDL->zeroes( double, $mask->dims );

    unless( defined $opt->{converged} ) {

	$opt->{converged} = sub {

	    my ( $last, $current ) = @_;

	    # stop if standard deviations and centers haven't changed
	    my $done = $current->{sigma} == $last->{sigma};

	    $done += $current->{center}[$_] == $last->{center}[$_]
	      foreach 0..$ndim-1;

	    return 1
	      if $done == $ndim + 1;

	    # or, if a tolerance was defined, stop if distance from old
	    # to new centers is less than the tolerance.
	    if ( defined $opt->{dtol} )
	    {
		my $dist2 = 0;

		$dist2 += ( $last->{center}[$_] - $current->{center}[$_] ) ** 2
		  foreach 0..$ndim-1;

		return 1 if $dist2 <= $opt->{dtol} **2;
	    }

	    return 0;
	}
    }


    #---------------------------------------------------------------

    # Iterate until convergence.

    my $iter = 0;

    my @iterations;

    # set up initial state

    # initial center
    _centroid( $is_sparse, $data, $wt, \@center, $mask );

    # generate mask
    my $nelem = _mask( $clip2, $data, $mask, $opt->{mask}, $r2, $is_sparse, $wt, \@center );

    barf( "too few elements ($nelem) after initial clip\n" )
      if $nelem  < 2;

    # now the variance, which is wt*dist^2 / sum(wt);
    my $variance = ($mask * $r2 )->sum / $mask->sum;

    push @iterations, {
		       center => [@center ],
		       iter => $iter,
		       nelem => $mask->sum,
		       sigma => sqrt( $variance ),
		       ( defined $clip2 ? (clip => sqrt($clip2)) : () ),

		      };

    $opt->{log} && $opt->{log}->( $iterations[-1] );

    my $done = 0;
    while ( $opt->{iterlim}-- && ! $done )
    {
	$iter++;

	# derive new clipping radius based upon variance within old radius
	my $last_variance = $variance;
	$clip2 = $nvar * $last_variance;

	# generate mask
	my $nelem = _mask( $clip2, $data, $mask, $opt->{mask}, $r2, $is_sparse, $wt, \@center );

	$nelem = $mask->sum;
	last if $nelem == 0;

	@center = _centroid( $is_sparse, $data, $wt, [ ( undef ) x $ndim ], $mask );

	# now the variance, which is wt*dist^2 / sum(wt);
	$variance = ($mask * $r2 )->sum / $mask->sum;

	push @iterations, { iter => $iter,
			    nelem => $nelem,
			    clip => sqrt($clip2),
			    sigma => sqrt($variance),
			    center => [@center],
			  };

	$opt->{log}->( $iterations[-1] );

	$done = $opt->{converged}->( $iterations[-2], $iterations[-1] );
    }

    my %results = %{ $iterations[-1] };

    $results{iterations} = \@iterations;
    $results{result} = $done;
    $results{opt} = { map { $_ => $opt->{$_} }
		      grep { ! ref $opt->{$_} }
		      keys %$opt };


    return wantarray()
      ? ($done ? @center : undef )
      : \%results;
}

sub _mask {

    my ( $clip2, $data, $mask, $mask_base, $r2, $is_sparse, $wt, $center ) = @_;

    # determine distance**2 to current center
    _dist2( $r2, $is_sparse, $data, $center );

    if ( defined $clip2 ) {

	$mask .=  defined $mask_base
	       ? ( $r2 <= $clip2 ) & $mask_base
	       : ( $r2 <= $clip2 );
    }

    my $nelem = $mask->sum;

    # generate a masked weight
    $mask *= $wt if defined $wt;

    return $nelem;
}

sub _dist2
{
  my ( $r2, $is_sparse, $data, $center ) = @_;

  if ( $is_sparse )
  {
    $r2 .= 0;
    $r2 += ( $data->[$_] - $center->[$_] ) ** 2 for 0..@$center-1;
  }

  else
  {
    $r2 .= $data->rvals( { Center => $center, Squared => 1 } );
  }
}

sub _centroid
{
  my ( $is_sparse, $data, $wt, $center, $bitmask ) = @_;


  # get centroid for those axes not specified
  if ( $is_sparse )
  {
    my $mdata;
    $mdata = defined $wt ? $bitmask * $wt : $bitmask;

    my $tot_wt = $mdata->sum;

    for my $axis ( 0 .. @$center-1 )
    {
      next if defined $center->[$axis];

      $center->[$axis] = ($data->[$axis] * $mdata )->sum / $tot_wt
    }
  }

  else
  {
    my $mdata = $data * $bitmask;

    my $tot_wt = $mdata->sum;

    for my $axis ( 0 .. @$center-1 )
    {
      next if defined $center->[$axis];
      $center->[$axis] = ( $mdata * $mdata->axisvals($axis) )->sum / $tot_wt;
    }
  }

  return @$center;
}

1;

