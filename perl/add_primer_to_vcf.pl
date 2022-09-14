use strict;
use warnings;
use utf8;

my ($QUERY_DIR, $BLAST_DIR, $INVCF, $OUTVCF);
my $line;
my @data = ();

if (@ARGV == 4){
  $QUERY_DIR = $ARGV[0];
  $BLAST_DIR = $ARGV[1];
  $INVCF = $ARGV[2];
  $OUTVCF = $ARGV[3];
}else{
  print "4 arguments are needed.\n";
  exit(1);
}

open my $fh, '<', $QUERY_DIR."/index.txt"
  or die "Can not open file ${QUERY_DIR}/index.txt.\n";
while($line = <$fh>){
  chomp($line);
  my @tmp =();
  my $fh2;
  unless(open $fh2, '<', $QUERY_DIR."/result".$line."_query.fasta") {
    push(@data, 0); #"0" will be flag for "no primer was found"
    next;
  }
  #When the query file exists (= some primers were found by primer3)
  #collect data of query.
  my $line2;
  my @query = ();
  while($line2 = <$fh2>){
    chomp($line2);
    if(substr($line2,0,1) eq ">"){
      push(@query, substr($line2,1));
    }
  }
  close $fh2;

  #collect data of blastn.
  open my $fh3, '<', $BLAST_DIR."/result".$line."_specif_filtered.txt"
    or die "Can not open file ${BLAST_DIR}/result${line}_specif_filtered.txt.\n";
  my $line3;
  my @specif = ();
  while($line3 = <$fh3>){
    chomp($line3);
    my @tmp3 = split(/\t/, $line3);
    push(@specif, \@tmp3);
  }
  close $fh3;

  my $nset = @query; #Number of primer set
  my $nrow = @specif;
  my $use_prmr = -1; #"-1" means "No adequete primer set"
  #If Useful primer set was found, this parameter will be natural number.
  for (my $i = 0; $i < $nset; $i++){
    my $hit_count = 0;
    my $des_count = 0;
    my $num = -1;
    for (my $j = 0; $j < $nrow; $j++){
      if($query[$i] eq $specif[$j][0]){
        $hit_count++;
        #Obtain chr name and position by parsing the name string.
        my @str1 = split(/:/, $specif[$j][0]);
        my @str2 = split(/_/, $str1[1]);
        my @pos = split(/-/, $str2[0]); #[0] start position, [1]end position.
        my $chr = $str1[0]; #chromosome
        $num = $str2[1]; #primer candidate number
        if($specif[$j][14] eq "plus"){
          if($specif[$j][8] == $pos[0]){
            $des_count++;
          }
        }else{ #if($specif[$j][14] eq "minus")
          if($specif[$j][8] == $pos[1]){
            $des_count++;
          }
        }
      }
    }

    #If the sequences were hit in only desire position, use the primer set.
    if($hit_count == 2 && $des_count == 2){
      $use_prmr = $num;
      last;
    }
  }

  if($use_prmr == -1){
    push(@data, 1); #"1" means "primer set was found but they are all non-specific"
    next;
  }

  open my $fh4, '<', $QUERY_DIR."/result".$line.".txt"
    or die "Can not open file ${QUERY_DIR}/result${line}.txt.\n";
  my $line4;
  my @primer_data = ();

  while($line4 = <$fh4>){
    chomp($line4);
    my @tmp4 = split(/=/, $line4);
    unless(defined $tmp4[0]){
      last;
    }
    if($tmp4[0] eq "SEQUENCE_ID"){
      push(@primer_data, substr($tmp4[1], 1));
    }
    if($tmp4[0] eq "PRIMER_LEFT_${use_prmr}_SEQUENCE") {
      push(@primer_data, $tmp4[1]);
    }
    if($tmp4[0] eq "PRIMER_RIGHT_${use_prmr}_SEQUENCE") {
      push(@primer_data, $tmp4[1]);
    }
    if($tmp4[0] eq "PRIMER_LEFT_${use_prmr}_TM") {
      push(@primer_data, $tmp4[1]);
    }
    if($tmp4[0] eq "PRIMER_RIGHT_${use_prmr}_TM") {
      push(@primer_data, $tmp4[1]);
    }
    if($tmp4[0] eq "PRIMER_PAIR_${use_prmr}_PRODUCT_SIZE") {
      push(@primer_data, $tmp4[1]);
    }
  }
  close $fh4;
  push(@data, \@primer_data);

}
close $fh;

#write primer data to input vcf.
open my $fh_in, '<', $INVCF
  or die "Can not open file ${INVCF}.\n";
open my $fh_out, '>', $OUTVCF;

my $header_flag = 0;
my $count = 0;
my $c_no_prmr = 0;
my $c_non_spcif = 0;
my $success = 0;
while($line = <$fh_in>){
  if(substr($line,0,1) eq "#"){
    #add information to end of ##INFO row in header.
    if($header_flag){
      if(substr($line,0,6) ne "##INFO"){
        print $fh_out "##INFO=<ID=PRIMER,Number=,Type=String,Description=\"Specific primer: \'ID | Left primer | Right primer | Left TM | Right TM | Product length\'\">\n";
        $header_flag = 0;
      }
    }else{
      if(substr($line,0,6) eq "##INFO"){
        $header_flag = 1;
      }
    }
    print $fh_out $line;
    next;
  }

  chomp($line);
  my @tmp = split(/\t/, $line);

  if($data[$count] == 0){
    $tmp[6] = "NO_PRIMER";
    $c_no_prmr++;
  }elsif($data[$count] == 1){
    $tmp[6] = "NON_SPECIFIC";
    $c_non_spcif++;
  }else{
    my $data_ref = $data[$count];
    $tmp[7] = $tmp[7].";PRIMER=".join("|", @$data_ref);
    $success++;
  }
  print $fh_out join("\t", @tmp)."\n";
  $count++;
}
close $fh_in;
close $fh_out;

print "Making VCF with primer set has finished.\n";
print "${success} variants have primer, ${c_no_prmr} variants have no primer candidate, ${c_non_spcif} variants have only non-specific primer in ${count} variants.\n";
