use strict;
use warnings;
use inc::Module::Install;

name 'CXC-PDL-Algorithm-Center';

all_from 'lib/CXC/PDL/Algorithm/Center.pm';

requires 'PDL';
requires 'Img2D::Events';
requires 'PDL::Image2D';
requires 'PDL::Transform';

build_requires 'Test::More';

auto_install;

WriteAll;
