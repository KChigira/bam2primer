use strict;
use warnings;
use utf8;

my ($IN, $OUT);
my ($min_dp, $max_dp);
my $between;
my ($col_x, $col_y);
my @pos_chr;
my (@pos_s, @pos_e);
my $pos_num;
my $line;
my @tmp_1 = ();
my @tmp_2 = ();
my @tmp_3 = ();
my (@genotype, @depth);
my $current_chr;
my $cnt = 0;

if (@ARGV == 10){
  $IN = $ARGV[0];
  $OUT = $ARGV[1];
  $min_dp = $ARGV[2];
  $max_dp = $ARGV[3];
  $between = $ARGV[4];
  $col_x = $ARGV[5];
  $col_y = $ARGV[6];
  @pos_chr = split(/,/, $ARGV[7]);
  @pos_s = split(/,/, $ARGV[8]);
  @pos_e = split(/,/, $ARGV[9]);
  $pos_num = @pos_chr;
}else{
  print "10 arguments are needed.\n";
  exit(1);
}

open my $fh_in, '<', $IN
  or die "Can not open file ${IN}.";
open my $fh_out, '>', $OUT;

while($line = <$fh_in>){
  if(substr($line,0,1) eq "#"){
    print $fh_out $line;
    next;
  }

  chomp($line);
  @tmp_3 = @tmp_2;  # 1 row forward line
  @tmp_2 = @tmp_1;  # a line checked in this turn
  @tmp_1 = split(/\t/, $line);  # next line

  #check the vcf file availavility
  my $tmp_1_len = @tmp_1;
  if($tmp_1_len < 11) {
    print "VCF file must have geno data of 2 or more samples.\n";
    exit(1);
  }

  #This is for 1st turn
  if(!(defined $current_chr)) {
    $current_chr = $tmp_1[0];
  }
  #Empty array means FALSE.
  #In 1st or 2nd turn, we cannot evaluate, so skip
  if (!(@tmp_2 && @tmp_3)) {next;}
  #when chromosome was changed, reset temporary values.
  if ($tmp_1[0] ne $current_chr) {
    @tmp_2 = ();
    @tmp_3 = ();
    $current_chr = $tmp_1[0];
    next;
  }

  #if it's multi allelic, next.
  my @allele = split(/,/, $tmp_2[4]);
  my $len = @allele;
  if($len != 1) {next;}

  #if it's indel, skip
  if((length $tmp_2[3] != 1) || (length $tmp_2[4] != 1)) {next;}

  #if filter is not "PASS", next.
  if($tmp_2[6] ne "PASS") {next;}

  #if variants in [$BETWEEN]bp front or back of the SNP exists, continued.
  if(($tmp_1[1] - $tmp_2[1]) <= $between ||
     ($tmp_2[1] - $tmp_3[1]) <= $between) {next;}

  #Select only designated chromosome and position
  my $flag = 0;
  for (my $i = 0; $i < $pos_num; $i++){
    if($tmp_2[0] eq $pos_chr[$i] && $tmp_2[1] >= $pos_s[$i] && $tmp_2[1] <= $pos_e[$i]) {
      $flag = 1;
    }
  }
  if($flag == 0) {next;}

  #extract genotype and depth
  my @geno_format = split(/:/, $tmp_2[8]);
  my $formatlen = @geno_format;
  my @geno_a = split(/:/, $tmp_2[${col_x}]);
  my @geno_b = split(/:/, $tmp_2[${col_y}]);
  for (my $i = 0; $i < $formatlen; $i++){
    if($geno_format[$i] eq "GT"){
      @genotype = ($geno_a[$i], $geno_b[$i]);
    }elsif($geno_format[$i] eq "DP"){
      @depth = ($geno_a[$i], $geno_b[$i]);
    }
  }

  #when genotype is different between samples
  #and depth is within the range, write output vcf
  if(($genotype[0] eq "0/0" && $genotype[1] eq "1/1") ||
     ($genotype[0] eq "1/1" && $genotype[1] eq "0/0")){
    if($depth[0] >= $min_dp && $depth[0] <= $max_dp &&
       $depth[1] >= $min_dp && $depth[1] <= $max_dp){
         print $fh_out join("\t", @tmp_2)."\n";
         $cnt = $cnt + 1;
    }
  }
}

close $fh_in;
close $fh_out;

print "${cnt} variants are selected.\n";
