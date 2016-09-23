#! perl

requires 'PDL';
requires 'PDL::Image2D';
requires 'PDL::Transform';
requires 'Ref::Util';
requires 'custom::failures';
requires 'Memoize';
requires 'Package::Stash';
requires 'Safe::Isa';
requires 'Scalar::Util';
requires 'Validate::Tiny';

on test => sub {

   requires 'Test::More';
   requires 'Test::Deep';
   requires 'Test::Fatal';
   requires 'PDL::GSL::RNG';

};

on develop => sub {

    requires 'Module::Install';
    requires 'Module::Install::AuthorRequires';
    requires 'Module::Install::AuthorTests';
    requires 'Module::Install::AutoLicense';
    requires 'Module::Install::CPANfile';
    requires 'Module::Install::ReadmeFromPod';

    requires 'Test::Fixme';
    requires 'Test::More';
    requires 'Test::NoBreakpoints';
    requires 'Test::Pod';
    requires 'Test::Pod::Coverage';
    requires 'Test::Perl::Critic';
    requires 'Test::CPAN::Changes';
    requires 'Test::CPAN::Meta';
    requires 'Test::CPAN::Meta::JSON';
};
