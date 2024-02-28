#!/usr/bin/perl
#
use strict;
use warnings;
use Cwd qw( abs_path cwd );
#
my $current_dir = cwd();
print "1: $current_dir\n";
#
my $script_path = abs_path($0);
print "1: $script_path\n";

