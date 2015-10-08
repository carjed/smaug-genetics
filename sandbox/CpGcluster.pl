#!/usr/bin/perl -w 

###########################################################################
###########################################################################

#              ********************************
#              **  CpGcluster (version 1.0)  **
#              ******************************** 

#   Laboratorio de Genómica Evolutiva y BioInformática
#    Universidad de Granada, Departamento de Genética

#              Web: http://bioinfo2.ugr.es
#              CGI: http://bioinfo2.ugr.es/CpGcluster


#  For questions, feedback, etc please contact to: José L. Oliver (oliver@ugr.es)
#                                                  Michael Hackenberg (mlhack@gmail.com)

# To see the options of the algorithm, please launch CpGcluster without any command line arguments
# If you use CpGcluster, please cite...

############################################################################
############################################################################
use strict;

my ($output, $sequence, $d, $plimit, $getd ,@dist_n);

&GetDefault();

################################# 
### Parameters ##################
#################################
my $maxN = 10; ## maximal number of Ns
#################################

print "\n      ***            Reading Sequence          ***\n";
(my $seqst, my $ID) = &GetFA($sequence);

############ Get Sequence features ##################
(my $seqlen, my $cpgnr) = &GetSeqFeat($seqst);
my $seqlenbruto = length($seqst);
my $Ndach = $seqlen-($cpgnr+1);
my $prob = $cpgnr/$Ndach;
my $CPGmax = int($cpgnr/2);

######################################################
print "      ***        Getting CpG Coordinates       ***\n";
# Get Coordinates of CpGs
my @cod = &GetCoords($seqst,"CG");

# sort numerically ascending
my @dd = sort {$a <=> $b} @dist_n;
$d = &GetPerc(\@dd,$getd);
print "\n      ***   Setting threshold distance to $d   ***\n\n";

### get protoislands
print "      ***         Detecting CpG clusters       ***\n";
my @protoislas = &GetProtoIslas(\@cod);

print "      ***          Calculating P-values        ***\n";
@protoislas = &CalcPvalNB(\@protoislas,$prob);

## Get Features like the obs/esp, clustering etc....
@protoislas = &GetCGI_features(\@protoislas);

&OUT(\@protoislas);



####################################################################
#######   SUBFUNCTIONS   ###########################################
###################################################################

sub GetDefault{
  print "\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------                                                         ---------\n";
  print "---------   Laboratorio de Genomica Evolutiva y BioInformatica    ---------\n";
  print "---------    Universidad de Granada, Departamento de Genetica     ---------\n";
print "---------                                                         ---------\n";
  print "---------               Web: http://bioinfo2.ugr.es               ---------\n";
  print "---------        CGI: http://bioinfo2.ugr.es/CpGcluster           ---------\n";
  print "---------                                                         ---------\n";
  print "---------              CpGcluster (1.0) 3/29/06                   ---------\n";
  print "---------                                                         ---------\n";
  print "---------------------------------------------------------------------------\n";
  print "---------------------------------------------------------------------------\n";
  print "\n";

  if($ARGV[0]){
    $sequence = $ARGV[0];
    if(-e $sequence){
    }
    else{
      die "Cannot find the input sequence: $sequence\n";
    }
  
    if($ARGV[1]){
      $getd = $ARGV[1];
      if($getd < 0 or $getd > 100){
	die "The Percentile must be between 0 and 100\n";
      }
    }
    else{
      print "Missing parameters!\n\n";
      print "Example for the usage of CpGcluster:\n";
      die "perl CpGcluster.pl [sequence]  [d]  [P-value]  [output]\n\n";
    }
    if($ARGV[2]){
      $plimit = $ARGV[2];
      if($plimit > 1){
	die "The maximal P-value you have choosen is higher than 1!\nPlease revise the order of the input parameters\n";
      }
    }
    else{
      print "Missing parameters!\n\n";
      print "Example for the usage of CpGcluster:\n";
      die "perl CpGcluster.pl [sequence]  [d]  [P-value]  [output]\n\n";
      
    }
    if($ARGV[3]){
      $output = $ARGV[3];
    }
    elsif($sequence){
      my @f = split(/\//,$sequence);
      my @f1 = split(/\./,$f[$#f]);
      $f[$#f] = "$f1[0].cpg";
      $output = join('/',@f);
    }
  }
  if(!($sequence)){

    print "Example for the usage of CpGcluster:\n\n";
    print "perl CpGcluster.pl [sequence]  [d]  [P-value]  [output]\n\n";
    
    print "sequence:   The input sequence in FASTA format\n\n";

    print "d:          The threshold distance on basis of a given percentile\n";
    print "            For example: d=25 calculates the percentile 25 of the genomic\n";
    print "            CpG distance distribution and takes this value as the threshold\n";
    print "            distance\n";
    print "            The recommended value is 50 (median distance)\n\n";
    print "P-value:    The maximal P-value under which a CpG cluster is considered as a\n";
    print "            CpG island\n";
    print "            The recommended limit is 1E-5\n\n";

    print "output:     The output file name\n";
    print "            If ommited; CpGcluster creates the output file in the\n"; 
    print "            input file directory with *.cpg extension\n\n";
    die "\n";
  }
  else{
    print "\n---------------------------------------------------------------------------\n";
    print "CpCcluster runs with the following parameters:\n\n";
    print "[sequence]    $sequence\n";
    print "[d]           $getd\n";
    print "[P-value]     $plimit\n";
    print "[output]      $output\n";
    print "---------------------------------------------------------------------------\n\n";
  }
}

# Read Sequence
sub GetFA{

  my $seqst_temp = "";
  open (I,$_[0]) or die "Can't open $_[0]";
  my $z = <I>;
  my $tes = substr($z,0,1);
  if($tes ne ">"){
    die "Sequence seems not to be in fasta format !!!!";
  }
  my @z = split(/\s+/," $z");
  $z = $z[1];
  $z=~ s/>//;
  $z =~ s/[\n\t\f\r\s]//g;
  my $ID_temp = $z;
  while($z = <I>){
    $z =~ s/[\n\t\f\r_0-9\s]//g;
    $seqst_temp .= $z;
  }
  return ($seqst_temp,$ID_temp);
}

sub OUT {
  my $c=0;

  open (OO,">$output") or die "could not open $output";
  
  while($_[0]->[$c]){
    printf OO "%d,%d\n",$_[0]->[$c]->[0],$_[0]->[$c]->[1]; 
    $c++;
  }
  close(OO);
  open(O,">$output-log.txt") or die "can't open $output-log.txt";
  print O "Basic statistics of the input sequence: $ID\n";
  printf O "Length: %d\n",$seqlen;
  printf O "Length without Ns: %d\n",$seqlenbruto;
  my $fg = $seqst =~ s/g/g/ig;
  my $fc = $seqst =~ s/c/c/ig;
  my $fa = $seqst =~ s/a/a/ig;
  my $ft = $seqst =~ s/t/t/ig;
  printf O "GC content: %0.3f\n",100*($fg+$fc)/$seqlen;
  printf O "Number of CpGs in sequence: %d\n",$cpgnr;
  printf O "Probability to find a CpG: %.4f\n\n",$prob;
  print O "Parameters used:\n";

  printf O "p-value threshold: $plimit\n";
  if($getd){
    print O "Distance threshold (percentile): $getd\n";
  }

  close(O);
}

## Get CpG cluster
sub GetProtoIslas{

  my @coord = @{$_[0]};
  my @t;
  my ($start, $end);
  my $des = "no";
  for(my $i = 0; $i <= $#coord - 1; $i++){
    
    my $dist = $coord[$i+1]  - ($coord[$i] + 1);

    if($dist <= $d){
      if($des eq "no"){
	$start = $coord[$i];
      }
      $end = $coord[$i+1]+1;
      $des = "yes";
    }
    elsif($dist > $d and $des eq "yes"){
      $des = "no";
      my @f = ($start, $end);
      push @t,\@f;
    }
  }
  if($des eq "yes"){
    my @f = ($start, $end);
    push @t,\@f;
  }
  return @t;
}

sub GetCGI_features{

  my @temp;
  my $c=0;
  while(defined($_[0]->[$c])){
    if($_[0]->[$c]->[4] < $plimit){

      my $len = $_[0]->[$c]->[1] - $_[0]->[$c]->[0] +1;
      (my $oe, my $cpgseq, my $gccont)= &CalcObsEsp($_[0]->[$c]->[0] -1,$len);
      my $coord1 = &GetCoord($cpgseq);
      (my $clust, my $meandist) = &GetClust($coord1);
      
      my $pval = $_[0]->[$c]->[4];
      $_[0]->[$c]->[4] = $oe;
      $_[0]->[$c]->[5] = $meandist;
      $_[0]->[$c]->[6] = $clust;
      $_[0]->[$c]->[7] = $gccont;
      $_[0]->[$c]->[8] = $pval;
      my @t = @{$_[0]->[$c]};
      push @temp,\@t;
    }
    $c++;
  }
  return @temp;
}

sub CalcObsEsp{

  my $cpgseq = substr ($seqst,$_[0],$_[1]);
  my $fc = $cpgseq =~ s/c/c/ig;
  my $fg = $cpgseq =~ s/g/g/ig;
  my $CGICpG = $cpgseq =~ s/CG/CG/ig;
  my $e = $fc*$fg;
  my $oet = $CGICpG*length($cpgseq)/$e;
  return ($oet,$cpgseq,$fc+$fg);
}

sub GetSeqFeat{

  my $n = $_[0];
  my $CpGnr = $n =~ s/CG/Cg/ig;
  my $NN = $n =~ s/N/N/ig;  
  my $seqlen = length($n)-$NN;
  return ($seqlen,$CpGnr);
}


sub GetCoords{

  my $n = $_[0];
  
  $n.="j";
  my @f =split(/$_[1]/i,$n);
  my @t;

  my $lencount = 0;
  for(my $i = 0; $i < $#f; $i++){
    $lencount += length($f[$i]);
    $t[$i] = $lencount + 1;
    $lencount+=2;
    my $nnr = $f[$i] =~ s/n/n/ig;
    if($nnr <= $maxN){
      push @dist_n,length($f[$i])+1;
    }
  }
  return @t;
}

sub GetPerc{

  my @t = @{$_[0]};
  my $totnr = @t;
  for(my $i = 0; $i <= $#t; $i++){
    if(100*$i/$totnr >= $_[1]){
      return $t[$i];
    }
  }
  return $t[$#t];
}


##########################################################################
############## SUBFUNCTIONS for P-value calculations  ####################
##########################################################################

######################################################################
## *** Calculates the negative binomial distribution

sub CalcPvalNB{

  my ($pval, @temp);
  my $c=0;
  my @islas = @{$_[0]};
  for(my $i = 0; $i <= $#islas; $i++){

    my $l = $islas[$i]->[1] - $islas[$i]->[0] + 1;
    my $str = substr ($seqst,$islas[$i]->[0]-1,$l);
    my $cpg = $str =~ s/cg/cg/ig;
    my $pval = &GetNB($l-(2*$cpg),$cpg-1,$_[1]);
#    print $cpg,"\n$str";<STDIN>;
    $pval = sprintf("%.5e",$pval);
    my @t = ($islas[$i]->[0],$islas[$i]->[1],$l,$cpg,$pval);
    push @temp, \@t;;
    $c++;
  }
  return @temp;
}
sub GetNB{
  
  my $pval = 0;
  for(my $j = 0; $j <= $_[0]; $j++){
    my $ptemp = &FactorialNB($j,$_[1]) + $_[1]*log($_[2]) + $j*log(1.0-$_[2]);
    $ptemp = exp($ptemp);
    $pval += $ptemp;
    }
  return $pval;
}


sub FactorialNB{
  my $stop = $_[0]+$_[1]-1;
  my $l1 = 0;
  my $l2 = 0;
  for(my $i = $_[0]+1;$i <= $stop; $i++){
    $l1 +=log($i);
  }
  for(my $i = 1;$i <= $_[1]-1; $i++){
    $l2 +=log($i);
  }
#  printf "$_[0] $_[1] %f\n",$l1-$l2;
  return $l1-$l2;
}
##########################################################################
############## SUBFUNCTIONS for Clustering Calculation  ####################
##########################################################################

sub GetClust{

  my @t = @{$_[0]};
  my $dist = &GetDist($_[0]);
  my $mean = &Normalize($dist);
  return (1, $mean);
}

sub Normalize{

  my @d = @{$_[0]};
  my $tot;
  for(my $i = 0; $i <= $#d; $i++){
    $tot+=$d[$i];
  }
  my $mean = $tot/@d;
  return $mean;
}

sub GetDist{

  my @dist;
  my @d = @{$_[0]};
  for(my $i = 0; $i < $#d; $i++){
    my @f = split (/\-/,$d[$i]);
    my @s = split (/\-/,$d[$i+1]);
    my $dist = $s[0]-$f[1];
    push @dist,$dist;
  }
  return \@dist;
}

sub GetCoord{

  my @coord;
  my @c = split (//,$_[0]);
  for(my $i = 0; $i < $#c; $i++){
    if($c[$i] =~ /[cC]/ and $c[$i+1] =~ /[gG]/){
      my $str = $i.'-'.eval($i+1);
      push @coord,$str;
    }
  }
  return \@coord;
}

