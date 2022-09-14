use strict;
use warnings;
use utf8;

unless (@ARGV == 1){
  print "1 argument is needed.\n";
  exit(1);
}

my $rev = reverse $ARGV[0];
$rev =~ tr/AaTtCcGg/TtAaGgCc/;
print $rev;
