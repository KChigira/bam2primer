use strict;
use warnings;
use utf8;

my $INVCF;
my $line;
my $column;

if (@ARGV == 1){
  $INVCF = $ARGV[0];
}else{
  print "1 arguments are needed.\n";
  exit(1);
}
(my $fn = $INVCF) =~ s/^(.*)\..*$/$1/; #remove extension
my $OUTVCF = "${fn}_available.vcf";

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
      print $fh_outvcf "##select_indels_by=select_primers_indel.pl\n";
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
my $num = @filtered;

foreach my $n (@filtered){
  print $fh_outvcf join("\t", @$n)."\n";
}
close $fh_outvcf;


my $OUTTXT = "${fn}_available.txt";
open my $fh_outtxt, '>', $OUTTXT;
my @tmp2 = split(/\t/, $column);
print $fh_outtxt "Chr\tPos\tID\tLeft\tRight\tL_TM\tR_TM\tProduct_Size_Ref\t".$tmp2[9]."\t".$tmp2[10]."\tInDel_Size\n";
for(my $i = 0; $i < $num; $i++){
  my @tmp3 = split(/;/, $filtered[$i][7]);
  my @o_line = ($filtered[$i][0], $filtered[$i][1]);
  foreach my $str (@tmp3){
    if(substr($str,0,6) eq "PRIMER"){
      push(@o_line, split(/\|/, substr($str,7)));
    }
  }
  my $prd_len_ref = $o_line[-1];
  my $indel_size = length($filtered[$i][3]) - length($filtered[$i][4]);
  #positive value means Deletion, negative value means insertion.
  if(substr($filtered[$i][9], 0, 3) eq "0/0"){
    push(@o_line, $prd_len_ref);
  }else{
    push(@o_line, $prd_len_ref - $indel_size);
  }
  if(substr($filtered[$i][10], 0, 3) eq "0/0"){
    push(@o_line, $prd_len_ref);
  }else{
    push(@o_line, $prd_len_ref - $indel_size);
  }
  push(@o_line, abs($indel_size));
  print $fh_outtxt join("\t", @o_line)."\n";
}
close $fh_outtxt;

print "${num} variants.\n";

