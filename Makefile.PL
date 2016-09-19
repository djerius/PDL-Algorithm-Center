use strict;
use warnings;
use inc::Module::Install;

name 'CXC-PDL-Algorithm-Center';

license 'gpl3';

auto_license( holder => 'Smithsonian Astrophysical Observatory' );

all_from 'lib/CXC/PDL/Algorithm/Center.pm';

readme_from;
readme_from '', { format => 'md' };;

author_tests( 'xt' );
cpanfile;

WriteAll;
