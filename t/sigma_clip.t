#!perl

use Test::More tests => 1;

use CXC::PDL::Algorithm::Center qw[ sigma_clip ];

use PDL::GSL::RNG;

my $rng = PDL::GSL::RNG->new('taus');

my ($x, $y ) = map { $_ + 0.2 } $rng->ran_bivariate_gaussian( 1, 2, 0.5, 100000 )->xchg(0,1)->dog;

eval {
    my ( $xc, $yc ) = sigma_clip( coords => [$x, $y],
				  center => [ 0.3, 0.5 ],
				  dtol => 0.0001
				);

};

ok( ! $@ );
