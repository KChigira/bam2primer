use strict;
use warnings;
use utf8;

my $line;
my $fasta;
my $template;
my $scope;
my $margin;
my $outdir;

my @seq_name = ();
my @sequence = ();
my @fn =();
my $seq;
my $nrow = 0;
my $target_start;
my $target_length;

if (@ARGV == 5){
  $fasta = $ARGV[0];
  $template = $ARGV[1];
  $scope = $ARGV[2];
  $margin = $ARGV[3];
  $outdir = $ARGV[4];
}else{
  print "5 arguments are needed.\n";
  exit(1);
}

open my $fh_fasta, '<', ${fasta}
  or die "Can not open fasta file";
#read sequence from fasta input
while($line = <$fh_fasta>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    push(@seq_name, $line);
    if(defined $seq){
      push(@sequence, $seq);
    }
    $seq = "";
  }else{
    $seq = ${seq}.${line};
  }
}
if(defined $seq){
  push(@sequence, $seq);
}
close $fh_fasta;


my $num_target = @seq_name;
open my $fh_idx, '>', ${outdir}."/index.txt";

for (my $i = 0; $i < $num_target; $i++){
  open my $fh_temp, '<', ${template}
    or die "Can not open template file";
  my $fn = "format".sprintf("%08d", $i).".txt";
  open my $fh_out, '>', ${outdir}."/".$fn;
  print $fh_idx sprintf("%08d", $i)."\n";

  #caliculate target sequence start position and length
  $target_start = $scope - $margin + 1;
  $target_length = (length $sequence[$i]) - ($scope * 2) + ($margin * 2);

  $nrow = 0;
  while($line = <$fh_temp>){
    $nrow = $nrow + 1;
    chomp($line);

    if($nrow == 1){
      print $fh_out $line.$seq_name[$i]."\n";
    }elsif($nrow == 2){
      print $fh_out $line.$sequence[$i]."\n";
    }elsif($nrow == 3){
      print $fh_out "${line}${target_start},${target_length}\n";
  #print substr($seq,$target_start-1,$target_length)."\n";
    }else{
      print $fh_out "${line}\n";
    }
  }
  close $fh_temp;
  close $fh_out;
}

close $fh_idx;
