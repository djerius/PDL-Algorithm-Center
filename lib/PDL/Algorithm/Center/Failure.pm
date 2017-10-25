package PDL::Algorithm::Center::Failure;

use strict;
use warnings;

our $VERSION = '0.05';

use custom::failures();

use Package::Stash;

use Exporter 'import';
our @EXPORT;

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

        $name = "${name}_failure";

        $stash->add_symbol( "&$name", sub () { __PACKAGE__ . "::$failure" } );

        push @EXPORT, $name;
    }

}

1;
