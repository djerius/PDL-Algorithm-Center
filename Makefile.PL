use strict;
use warnings;
use inc::Module::Install;

name 'PDLx-Algorithm-Center';

license 'gpl3';

auto_license( holder => 'Smithsonian Astrophysical Observatory' );

all_from 'lib/PDLx/Algorithm/Center.pm';

readme_from;
readme_from '', { format => 'md' };;

resources (
   license    => 'http://www.gnu.org/licenses/gpl-3.0.html',
   repository => 'https://github.com/djerius/PDLx-Algorithm-Center',
   bugtracker => 'https://rt.cpan.org/Public/Dist/Display.html?Name=PDLx-Algorithm-Center',
);

author_tests( 'xt' );
cpanfile;

WriteAll;
