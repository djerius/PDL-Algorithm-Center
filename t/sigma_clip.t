#!perl

use strict;
use warnings;


use Test::More;
use Test::Fatal;
use Test::Deep;

use CXC::PDL::Algorithm::Center qw[ sigma_clip ];

use PDL;
use PDL::GSL::RNG;


sub eclass ($) {
    return join( '::', 'CXC::PDL::Algorithm::Center::failure', @_ );
}

# log iterations for debugging poiposes
sub logit {
    my %msg = %{ shift() };
    $msg{center} = [ PDL::Core::topdl( $msg{center} )->list ];
    note explain \%msg;
}

########################################
# interface

my %req = ( nsigma => 1.5 );

# coordinates
subtest "coords" => sub {

    isa_ok( exception { sigma_clip( %req, coords => PDL->null ) },
        eclass( 'parameter::type' ), 'null' );

    isa_ok( exception { sigma_clip( %req, coords => PDL->zeroes( 0 ) ) },
        eclass( 'parameter::type' ), 'empty' );

    isa_ok( exception { sigma_clip( %req, coords => 'scalar' ) },
        eclass( 'parameter::type' ), 'scalar' );

    isa_ok(
        exception { sigma_clip( %req, coords => ['scalar'] ) },
        eclass( 'parameter::type' ),
        'array element not piddle'
    );

    isa_ok(
        exception {
            sigma_clip( %req, coords => [ pdl( 1 ), pdl( 1, 2 ), pdl( 1 ) ] )
        },
        eclass( 'parameter::type' ),
        'unequal number of elements'
    );

    isa_ok(
        exception {
            sigma_clip( %req, coords => [ pdl( [ 1, 2 ], [ 3, 4 ] ) ] )
        },
        eclass( 'parameter::type' ),
        'not 1D'
    );

    is( exception { sigma_clip( %req, coords => pdl( 1, 2 ) ) },
        undef, "pdl(1, 2)" );

    is( exception { sigma_clip( %req, coords => [ pdl( 1, 2 ) ] ) },
        undef, "[ pdl(1, 2) ]" );

    is( exception { sigma_clip( %req, coords => pdl( 1 ) ) }, undef, "pdl(1)" );


};

subtest "center" => sub {

    isa_ok(
        exception {
            sigma_clip(
                %req,
                coords => pdl( 1 ),
                center => 'foo'
              )
        },
        eclass( 'parameter::type' ),
        'center not a 1D piddle'
    );

    is(
        exception {
            sigma_clip(
                %req,
                coords => pdl( [1], [2] ),
                center => pdl( 1, 1 ) )
        },
        undef,
        'center a 1D piddle'
    );


};

subtest "weight" => sub {

    isa_ok(
        exception {
            sigma_clip( %req, weight => [] );
        },
        eclass( 'parameter::type' ),
        'wrong type'
    );

    is(
        exception {
            sigma_clip( %req, weight => pdl( 1 ) );
        },
        undef,
        'pdl(1)'
    );

};



foreach my $field ( 'weight', 'mask' ) {

    subtest "coords + $field" => sub {

        isa_ok(
            exception {
                sigma_clip(
                    %req,
                    coords => [ pdl( 1, 2, 3 ), pdl( 3, 4, 5 ) ],
                    $field => pdl( 1 ),
                  )
            },
            eclass( 'parameter::dimension' ),
            'incorrect dimensions'
        );

        is(
            exception {
                sigma_clip(
                    %req,
                    coords => [ pdl( 1, 2, 3 ), pdl( 3, 4, 5 ) ],
                    $field => pdl( 1,   2, 3 ),
                  )
            },
            undef,
            'matched dimensions'
        );

    };

}

subtest "!( coords || weight)" => sub {

    isa_ok(
        exception { sigma_clip( %req, ) },
        eclass( 'parameter::missing' ),
        'neither coords nor weight'
    );

};


subtest 'log' => sub {

    isa_ok(
        exception {
            sigma_clip(
                %req,
                coords => pdl( 1 ),
                log    => 'string',
              )
        },
        eclass( 'parameter::type' ),
        'log = string'
    );

};


for my $field ( qw( clip nsigma dtol ) ) {

    subtest $field => sub {

        isa_ok(
            exception {
                sigma_clip( %req, coords => pdl( 1 ), $field => 'foo' )
            },
            eclass( 'parameter::value' ),
            qq/$field = string/
        );

        isa_ok(
            exception {
                sigma_clip( %req, coords => pdl( 1 ), $field => -1 )
            },
            eclass( 'parameter::value' ),
            qq/$field = -1/
        );

        isa_ok(
            exception {
                sigma_clip( %req, coords => pdl( 1 ), $field => 0 )
            },
            eclass( 'parameter::value' ),
            qq/$field = 0/
        );

        is(
            exception {
                sigma_clip( %req, coords => pdl( 1, 1 ), $field => 1 )
            },
            undef,
            qq/$field = 1/
        );

    };


}

########################################
# operations

# let's try clipping!
subtest 'coords + clip results' => sub {

    my $rng = PDL::GSL::RNG->new( 'taus' );

    # so tests are reproducible
    my $nelem = 100000;
    $rng->set_seed( 1 );
    srand( 1 );

    # generate a bunch of coordinates
    my $coords = $rng->ran_bivariate_gaussian( 10, 8, 0.5, $nelem );
    my $center = pdl( 0.5, 0.5 );

    # calculate sigma for those inside of a radius of 10
    my $inside = $coords->xchg( 0, 1 )
      ->whereND( dsumover( ( $coords - $center )**2 ) < 100 )->xchg( 0, 1 );

    my $ninside = $inside->dim( 1 );
    my $sigma = sqrt( dsum( ( $inside - $center )**2 ) / $ninside );

    my $results = sigma_clip(
        coords  => $coords,
        center  => $center,
        clip    => 10,
        dtol    => 0.00001,
        iterlim => 100,
        nsigma  => 1.5,
    );

    ok( $results->{success}, "successful centering" );

    # make sure iteration 0 agrees with the above calculations
    cmp_deeply(
        $results->{iterations}[0],
        methods(
            sigma  => num( $sigma, .0000001 ),
            nelem  => $ninside,
            weight => $ninside,
        ),
        "iteration 0",
    );

    # and that the last one agrees with previous fiducial runs, to see
    # if something has broken
    cmp_deeply(
        $results->{iterations}[-1],
        methods(
            iter   => 56,
            dist   => 0,
            nelem  => 43597,
            weight => 43597,
            center => code(
                sub {
                    all(
                        $_[0]->approx(
                            pdl( 0.0126755884280886, 0.0337090322699186 ) ) );
                } )
        ),
        "iteration -1",
    );


};

# let's try masking!
subtest 'coords + mask results' => sub {

    my $rng = PDL::GSL::RNG->new( 'taus' );

    # so tests are reproducible
    my $nelem = 100000;
    $rng->set_seed( 1 );
    srand( 1 );

    # generate a bunch of coordinates
    my $coords = $rng->ran_bivariate_gaussian( 10, 8, 0.5, $nelem );

    # calculate sigma for those inside of a radius of 10
    my $center = pdl( 0.5, 0.5 );
    my $mask = dsumover( ( $coords - $center )**2 ) < 100;
    my $inside = $coords->xchg( 0, 1 )->whereND( $mask )->xchg( 0, 1 );

    my $ninside = $inside->dim( 1 );
    my $sigma = sqrt( dsum( ( $inside - $center )**2 ) / $ninside );

    my $results = sigma_clip(
        coords  => $coords,
        center  => $center,
        mask    => $mask,
        iterlim => 100,
        dtol    => 0.00001,
        nsigma  => 1.5,
    );

    # make sure iteration 0 agrees with the above calculations
    cmp_deeply(
        $results->{iterations}[0],
        methods(
            sigma  => num( $sigma, .0000001 ),
            nelem  => $ninside,
            weight => $ninside,
        ),
        "iteration 0",
    );

    # and that the last one agrees with previous fiducial runs, to see
    # if something has broken
    cmp_deeply(
        $results->{iterations}[-1],
        methods(
            iter   => 56,
            dist   => 0,
            nelem  => 43597,
            weight => 43597,
            center => code(
                sub {
                    all(
                        $_[0]->approx(
                            pdl( 0.0126755884280886, 0.0337090322699186 ) ) );
                } )
        ),
        "iteration -1",
    );

};

# Let's try weighting!

subtest 'coords + clip + weight results' => sub {

    my $rng = PDL::GSL::RNG->new( 'taus' );

    # so tests are reproducible
    my $nelem = 100000;
    $rng->set_seed( 1 );
    srand( 1 );

    # generate a bunch of coordinates
    my $center = pdl( 0.5, 0.5 );
    my $coords = $rng->ran_bivariate_gaussian( 10, 8, 0.5, $nelem );

    # calculate sigma for those inside of a radius of 10
    my $weight       = zeroes( $nelem ) + 2;
    my $mask         = dsumover( ( $coords - $center )**2 ) < 100;
    my $inside       = $coords->xchg( 0, 1 )->whereND( $mask )->xchg( 0, 1 );
    my $inside_nelem = $inside->dim( -1 );

    my $inside_weight     = $weight->where( $mask );
    my $inside_weight_sum = $inside_weight->dsum;

    my $sigma
      = sqrt( dsum( $inside_weight * dsumover( ( $inside - $center )**2 ) )
          / $inside_weight_sum );

    my $results = sigma_clip(
        coords  => $coords,
        center  => $center,
        nsigma  => 1.5,
        clip    => 10,
        iterlim => 100,
        dtol    => 0.00001,
        weight  => $weight,
    );

    # make sure iteration 0 agrees with the above calculations
    cmp_deeply(
        $results->{iterations}[0],
        methods(
            sigma  => num( $sigma, .0000001 ),
            nelem  => $inside_nelem,
            weight => $inside_weight_sum,
        ),
        "iteration 0",
    );

    # and that the last one agrees with previous fiducial runs, to see
    # if something has broken
    cmp_deeply(
        $results->{iterations}[-1],
        methods(
            iter   => 56,
            dist   => 0,
            nelem  => 43597,
            weight => 43597 * 2,
            center => code(
                sub {
                    all(
                        $_[0]->approx(
                            pdl( 0.0126755884280886, 0.0337090322699186 ) ) );
                } )
        ),
        "iteration 0",
    );

};


done_testing;

