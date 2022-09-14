use strict;
use warnings;
use utf8;

my $line;
my $input;

if (@ARGV == 1){
  $input = $ARGV[0];
}else{
  print "1 arguments are needed.\n";
  exit(1);
}

my $idnum = 0; #id number of primer candidate, from 0
my $chr;
my $start_pos;
my @fwd = ();
my @rev = ();
my @fwd_pos = ();
my @rev_pos = ();

open my $fh_in, '<', $input
  or die "Can not open input file";
while($line = <$fh_in>){
  chomp($line);
  my @tmp = split(/=/, $line);
  unless(defined $tmp[0]){
    last;
  }
  if($tmp[0] eq "SEQUENCE_ID"){
    my @tmp2 = split(/:/, $tmp[1]);
    $chr = $tmp2[0];
    my @tmp3 = split(/-/, $tmp2[1]);
    $start_pos = $tmp3[0];
  }
  if($tmp[0] eq "PRIMER_LEFT_".${idnum}."_SEQUENCE") {
    push(@fwd, $tmp[1]);
  }
  if($tmp[0] eq "PRIMER_RIGHT_".${idnum}."_SEQUENCE") {
    push(@rev, $tmp[1]);
  }
  if($tmp[0] eq "PRIMER_LEFT_".${idnum}) {
    my @tmp2 = split(/,/, $tmp[1]);
    push(@fwd_pos, $tmp2[0]);
  }
  if($tmp[0] eq "PRIMER_RIGHT_".${idnum}) {
    my @tmp2 = split(/,/, $tmp[1]);
    push(@rev_pos, $tmp2[0]);
    ${idnum}++;
  }
}
close $fh_in;

my $check = @fwd;
if($check == 0){
  exit(0);
}

(my $fn = $input) =~ s/^(.*)\..*$/$1/; #remove extension
my $output = $fn."_query.fasta";
open my $fh_out, '>', $output;

for (my $i = 0; $i < $check; $i++){
  my $fwd_pos_all = $start_pos + $fwd_pos[$i];
  my $rev_pos_all = $start_pos + $rev_pos[$i];
  print $fh_out $chr.":".$fwd_pos_all."-".$rev_pos_all."_".$i."\n";
  print $fh_out $fwd[$i]."NNNNNNNNNNNNNNNNNNNN".$rev[$i]."\n";
}

close $fh_out;
