use strict;
use warnings;
use utf8;

my $line;
my $vcf;
my $fasta;
my $template;
my $outdir;

my @data_vcf = ();
my @seq_name = ();
my @sequence = ();
my @fn =();
my $seq;

my $scope;
my $margin;

if (@ARGV == 6){
  $vcf = $ARGV[0];
  $fasta = $ARGV[1];
  $template = $ARGV[2];
  $outdir = $ARGV[3];
  $scope = $ARGV[4];
  $margin = $ARGV[5];
}else{
  print "6 arguments are needed.\n";
  exit(1);
}


open my $fh_vcf, '<', ${vcf}
  or die "Can not open vcf file";
while($line = <$fh_vcf>){
  chomp($line);
  if(substr($line,0,1) eq "#"){
    next;
  }
  
  #extract information of (Ref. length, Alt. length, Former Obstacle, Later obstacle)
  #-> $data_vcf
  my @inf = ();
  
  my @tmp = split(/\t/, $line);
  push(@inf, length $tmp[3]);
  push(@inf, length $tmp[4]);
  
  my @tmp2 = split(/;/, $tmp[7]); #INFO column
  foreach my $obj (@tmp2){
    my @tmp3 = split(/=/, $obj);
    if($tmp3[0] eq "OBST"){
      my @tmp4 = split(/,/, $tmp3[1]);
      @inf = (@inf, @tmp4);
    }
  }
  
  push(@data_vcf, \@inf);
}
close $fh_vcf;


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
  my $target_start = $scope - $margin + 1;
  my $target_length = $data_vcf[$i][0];
  #caliculate excluded region
  my @exclude = ();
  if($data_vcf[$i][2] <= $scope){
    @exclude = (@exclude, 1, ($scope - $data_vcf[$i][2] + 1));
  }
  if($data_vcf[$i][3] <= $scope+$target_length-1){
    @exclude = (@exclude, ($scope + $target_length + $data_vcf[$i][3]), 
                          ($scope - $data_vcf[$i][3]) + 1);
  }
  #caliculate max, min, opt length of PCR products
  my ($max,$opt,$min);
  my $x = abs($data_vcf[$i][0] - $data_vcf[$i][1]);
  if($data_vcf[$i][0] > $data_vcf[$i][1]){
    #Deletion
    $max = int(300 * $x / (10 + $x));
    $opt = int(0.8 * 300 * $x / (10 + $x));
    if($x < 10){
      $min = 5 * $x + 25 + $x;
    }else{
      $min = 75 + $x;
    }
  }else{
    #Insertion
    $max = int(300 * $x / (10 + $x)) - $x;
    $opt = int(0.8 * 300 * $x / (10 + $x)) - $x;
    if($x < 10){
      $min = 5 * $x + 25;
    }else{
      $min = 75;
    }
  }

  my $nrow = 0;
  while($line = <$fh_temp>){
    $nrow = $nrow + 1;
    chomp($line);
    if($nrow == 1){
      print $fh_out $line.$seq_name[$i]."\n";
    }elsif($nrow == 2){
      print $fh_out $line.$sequence[$i]."\n";
    }elsif($nrow == 3){
      print $fh_out "${line}${target_start},${target_length}\n";
    }elsif($nrow == 4){
#      my $len_exclude = @exclude;
#      if($len_exclude == 2){
#        print $fh_out $line.$exclude[0].",".$exclude[1]."\n";
#      }elsif($len_exclude == 4){
#        print $fh_out $line.$exclude[0].",".$exclude[1]." ".$exclude[2].",".$exclude[3]."\n";
#      }
    }elsif($nrow == 13){
      print $fh_out $line.$opt."\n";
    }elsif($nrow == 14){
      print $fh_out $line.$min."-".$max."\n";
    }else{
      print $fh_out "${line}\n";
    }
  }
  close $fh_temp;
  close $fh_out;
}

close $fh_idx;
