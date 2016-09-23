package PDLx::Algorithm::Center::Validate;

use strict;
use warnings;

use Exporter 'import';

our @EXPORT = qw( is_a_positive_integer
		   is_a_positive_number
		   is_a_subroutine_reference
		   is_a
		   is_a_nonempty_PDL
		   is_a_ndim_PDL
);

use feature 'state';
use Memoize;

use Safe::Isa;
use Ref::Util qw<  is_coderef >;
use Scalar::Util 'looks_like_number';

sub is_a_positive_number {

    state $sub = sub {
	my $value = shift // return;
	return if looks_like_number( $value ) && $value > 0;
	return
	  [ 'parameter::value', "<$_[1]> is not a positive number" ];
    };

}

sub is_a_positive_integer {

    state $sub = sub {
	my $value = shift // return;
	return if looks_like_number( $value ) && $value > 0 && int($value) == $value;
	return [ 'parameter::value', "<$_[1]> is not a positive integer" ];
    };
}


sub is_a_subroutine_reference {

    state $sub = sub {
	my $value = shift // return;
	return if is_coderef( $value );
	return [ "parameter::type",  "<$_[1]> must be a subroutine referenece" ];
    };

}

memoize('is_a');
sub is_a {
    my $class = shift;

    return sub {
	my $value = shift // return;
	return if $value->$_isa( $class );
	return [ 'parameter::type', "<$_[1]> must be an object of type <$class>" ] ;
    };

}

memoize('is_a_ndim_PDL');
sub is_a_ndim_PDL {

    my ( $min, $max, $msg ) = @_;

    $msg //=
      "must be a non-empty piddle with dimensions "
      . defined $min && defined $max && $min == $max
      ? "== $min"
      : join( ' and ',
        defined $min ? ">= $min" : (),
        defined $max ? "<= $max" : (),
      );

    return sub {

        my $value = shift // return;
        return
             if $value->$_isa( 'PDL' )
	  && ! $value->isnull
	  && ! $value->isempty
          && ( ! defined $min || $value->ndims >= $min )
          && ( ! defined $max || $value->ndims <= $max );

        return [ "parameter::type", "<$_[1]> $msg" ];
      }
}


sub is_a_nonempty_PDL {

    state $sub = sub {

	my $value = shift // return;
	return if $value->$_isa( 'PDL' ) && ! ( $value->isnull || $value->isempty );
	return [ "parameter::type",  "<$_[1]> may not be empty or null" ];

    }
}


1;
