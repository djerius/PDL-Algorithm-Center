#!perl

use Test::More tests => 1;

use CXC::PDL::Algorithm::Center qw[ maxpix ];

use PDL::GSL::RNG;

my $rng = PDL::GSL::RNG->new('taus');

my ($x, $y ) = $rng->ran_bivariate_gaussian( 1, 2, 0.5, 100000 )->xchg(0,1)->dog;

eval {
    my ( $xc, $yc ) = maxpix( coords => [$x, $y], minmax => 1000 );
};

ok( ! $@ ) or diag $@;


