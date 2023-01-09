use strict;
use warnings;
use utf8;

my ($INVCF, $FAI, $TARGET_NUM, $INTERVAL,$col_x,$col_y);
my $line;
my $before;
my $column;

if (@ARGV == 6){
  $INVCF = $ARGV[0];
  $FAI = $ARGV[1];
  $TARGET_NUM = $ARGV[2];
  $INTERVAL = $ARGV[3];
  $col_x = $ARGV[4];
  $col_y = $ARGV[5];
}else{
  print "6 arguments are needed.\n";
  exit(1);
}
(my $fn = $INVCF) =~ s/^(.*)\..*$/$1/; #remove extension
my $OUTVCF = "${fn}_select_${TARGET_NUM}_mk_or_${INTERVAL}_bp.vcf";

open my $fh_invcf, '<', $INVCF
  or die "Can not open file ${INVCF}\n";
open my $fh_outvcf, '>', $OUTVCF;

my @filtered = ();
while($line = <$fh_invcf>){
  if(substr($line,0,1) eq "#"){
    #add information to end of header.
    unless(substr($line,0,2) eq "##"){
      $column = $line;
      chomp($column);
      print $fh_outvcf "##select_primer_by=select_primers_to_use.pl target number ${TARGET_NUM}; position interval ${INTERVAL}\n";
    }
    print $fh_outvcf $line;
    next;
  }

  chomp($line);
  my @tmp = split(/\t/, $line);
  #memory only varinats with primer.
  if($tmp[6] eq "PASS"){
    push(@filtered, \@tmp);
  }
}
close $fh_invcf;
$before = @filtered;

#read chromosome information
open my $fh_fai, '<', $FAI
  or die "Can not open file ${INVCF}\n";
my @chroms = ();
while($line = <$fh_fai>){
  chomp($line);
  my @tmp = split(/\t/, $line);
  push(@chroms, $tmp[0]);
}
close $fh_fai;

#Filter out closely spaced variants.
my $chr_num = @chroms;
my $num = @filtered;
my @spaces;
my $min_space = 9**9**9;
while(${num} > ${TARGET_NUM}){
  @spaces = (9**9**9);
  #9**9**9 must be over flow, so this means "Infinity"
  #length of @spaces will be $num + 1. these numbers indicates positional interval of each variants.
  $min_space = 9**9**9;
  my $min_index = 0;
  for(my $i = 0; $i < ($num - 1); $i++){
    if($filtered[$i][0] ne $filtered[$i+1][0]){
      push(@spaces, 9**9**9);
    }else{
      push(@spaces, $filtered[$i+1][1] - $filtered[$i][1]);
    }

    #check the interval is minimum or not.
    if($spaces[$i+1] < $min_space){
      $min_space = $spaces[$i+1];
      $min_index = $i + 1;
    }
  }
  push(@spaces, 9**9**9);

  if($min_space > $INTERVAL){
    last;
  }

  #Which of the two variants that make up the shortest interval?
  #the one with shorter interval to their opposite side variants will be deleted.
  if($spaces[$min_index - 1] < $spaces[$min_index + 1]){
    splice(@filtered, $min_index - 1, 1);
  }else{
    splice(@filtered, $min_index, 1);
  }
  $num = @filtered;
}

foreach my $n (@filtered){
  print $fh_outvcf join("\t", @$n)."\n";
}
close $fh_outvcf;


my $OUTTXT = "${fn}_select_${TARGET_NUM}_mk_or_${INTERVAL}_bp.txt";
open my $fh_outtxt, '>', $OUTTXT;
my @tmp2 = split(/\t/, $column);
print $fh_outtxt "Chr\tPos\tID\tLeft\tRight\tL_TM\tR_TM\tProduct_Size\t".$tmp2[$col_x]."\t".$tmp2[$col_y]."\n";
for(my $i = 0; $i < $num; $i++){
  my @tmp3 = split(/;/, $filtered[$i][7]);
  my @o_line = ($filtered[$i][0], $filtered[$i][1]);
  foreach my $str (@tmp3){
    if(substr($str,0,6) eq "PRIMER"){
      push(@o_line, split(/\|/, substr($str,7)));
    }
  }
  if(substr($filtered[$i][$col_x], 0, 3) eq "0/0"){
    push(@o_line, "Ref");
  }else{
    push(@o_line, "Alt");
  }
  if(substr($filtered[$i][$col_y], 0, 3) eq "0/0"){
    push(@o_line, "Ref");
  }else{
    push(@o_line, "Alt");
  }
  print $fh_outtxt join("\t", @o_line)."\n";
}
close $fh_outtxt;


print "Selecting variants in VCF with primer set has finished.\n";
print "${num} variants selected from ${before} variants.\n";
print "The shortest marker interval is ${min_space} bp.\n";
