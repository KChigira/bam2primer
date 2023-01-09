use strict;
#use warnings;
use utf8;

my $line;
my $input;
my $query;
my $max_mm; #max missmatch bases number. *default 5
my $max_mm_3t; #max missmatch bases number in 3' terminal end (5bp). *default 1
my $product; #Upper size limit for unintended PCR products. *default 4000
my $output;
my @q_data = ();
my @data = ();
my @candidate = ();

if (@ARGV == 5){
  $input = $ARGV[0];
  $query = $ARGV[1];
  $max_mm = $ARGV[2];
  $max_mm_3t = $ARGV[3];
  $product = $ARGV[4];
}else{
  print "5 arguments are needed.\n";
  exit(1);
}

open my $fh_q, '<', $query
  or die "Can not open query file";
while($line = <$fh_q>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    my @str = (substr($line,1), "");
    push(@q_data, \@str);
  }else{
    my $last = @q_data;
    if($last > 0){
      $q_data[-1][1] = $q_data[-1][1].$line;
    }
  }
}
close $fh_q;


my $primer_num = 0; #This is index number for @q_data.
my @q_seq = split(/N/, $q_data[0][1]);
my @p_pos_L = (1, length($q_seq[0]));
my @p_pos_R = ($p_pos_L[1]+20+1, $p_pos_L[1]+20+length($q_seq[-1]));

open my $fh_in, '<', $input
  or die "Can not open input file";
while($line = <$fh_in>){
  chomp($line);
  my @tmp = split(/\t/, $line);

  while($tmp[0] ne $q_data[$primer_num][0]){
    $primer_num++;
    @q_seq = split(/N/, $q_data[0][1]);
    @p_pos_L = (1, length($q_seq[0]));
    @p_pos_R = ($p_pos_L[1]+20+1, $p_pos_L[1]+20+length($q_seq[-1]));
  }

  my $mm;
  if($tmp[7] < ($p_pos_R[1] / 2)) {
    #$tmp[6] means the position of 5' terminal in query sequence.
    #$tmp[7] means the position of 3' terminal in query sequence.
    $mm = ($tmp[6] - $p_pos_L[0]) + ($p_pos_L[1] - $tmp[7]);
  }else{
    $mm = ($tmp[6] - $p_pos_R[0]) + ($p_pos_R[1] - $tmp[7]);
  }
  $mm = $mm + $tmp[4];
  #$tmp[4] means the number of missmach.
  if($mm > $max_mm){
    next;
  }

  my $mm_3t; #missmatch number of 3' terminal 5bp
  if($tmp[7] < ($p_pos_R[1] / 2)) {
    #$tmp[7] means the position of 3' terminal in query sequence.
    $mm_3t = $p_pos_L[1] - $tmp[7];
  }else{
    $mm_3t = $p_pos_R[1] - $tmp[7];
  }

  #compare characters of 5bp in 3' teminal.
  my @qseq = split(//, substr($tmp[12], -5, 5));
  my @sseq = split(//, substr($tmp[13], -5, 5));
  for (my $i = 0; $i < (5-$mm_3t); $i++){
    if($qseq[$i] ne $sseq[$i]){
      $mm_3t++;
    }
  }
  if($mm_3t > $max_mm_3t){
    next;
  }

  push(@tmp, $primer_num); #For sort by primer sets.
  push(@data, \@tmp);
}
close $fh_in;

@data = sort { $a->[8] <=> $b->[8] } @data; #sort by position
@data = sort { $a->[1] cmp $b->[1] } @data; #sort by chromosome
@data = sort { $a->[15] <=> $b->[15] } @data; #sort by primers
my $nrow = @data;

for (my $i = 0; $i < $nrow; $i++){
  if($data[$i][14] eq "plus"){
    my $pri = $data[$i][0];
    my $chr = $data[$i][1];
    my $pos = $data[$i][8];
    my $cnt = 0;
    my $flag = 1;
    while($flag){
      $cnt++;
      if($i + $cnt >= $nrow){
        $flag = 0;
      }elsif($data[$i + $cnt][0] ne $pri){
        $flag = 0;
      }elsif($data[$i + $cnt][1] ne $chr){
        $flag = 0;
      }elsif($data[$i + $cnt][14] eq "minus"){
        if($data[$i + $cnt][8] - $pos < $product){
          push(@candidate, $data[$i]);
          push(@candidate, $data[$i + $cnt]);
        }else{
          $flag = 0;
        }
      }
    }
  }
}

(my $fn = $input) =~ s/^(.*)\..*$/$1/; #remove extension
$output = $fn."_filtered.txt";
open my $fh_out, '>', $output;
foreach my $item(@candidate){
  print $fh_out join("\t", @$item)."\n";
}
close $fh_out;
