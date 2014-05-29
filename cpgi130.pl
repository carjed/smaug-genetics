## LICENCE AGREEMENT
##
## The intellectual property right of CpG island search script belongs 
## to Daiya Takai and Peter A. Jones.  Use of the this script is free for 
## academic users. Daiya Takai is responsible for the scientific content 
## of the script.
##
## Integration of this script into another script/program  will require
## explicit permission of authors. Such permission will only be granted 
## if there is a valid scientific or technical reason to encapsulate the 
## entire content of the script into a new resource. 
##
## Integration of any part of the script into any commercial product is only 
## permitted by the agreement with the authors.
##
## Publication of results obtained with the aid of the script should cite:
##
## Takai D and Jones PA. Comprehensive analysis of CpG islands in human 
## chromoseme 21 and 22. PNAS 2002 99(6):3740-5
##
## If you have questions, please send an email to the following address: 
## Daiya Takai : takai_d@ccnt.hsc.usc.edu
##
## This script is reviced on May 21st, 2002.

$usage = <<EOT;
usage cpgi130 [sequencefile] [GCC=X] [OE=X] [LENGTH=X] [HELP]
example cpgi130 NT_000000.seq LENGTH=200
EOT

$helpmenu = <<EOT;
 GCC          set \%GC of lower limit (50-70,default value:55\%)
 OE           set ObsCpG/ExpCpG of lower limit (0.60-1.00,default value:0.65)
 LENGTH       set length of lower limit (200-1500,default value:500bp)
EOT

$credit = <<EOT;
CpG island searcher command line version
Ver 1.3 released 05/21/03
by Takai D. & Jones PA.

EOT

$GCC=55;
$OE=0.65;
$LENGTH=500;

if (@ARGV == 0){
	print $usage;
	exit;
}
$filename= shift @ARGV;
	if ($filename =~ /HELP/i){
		print $helpmenu;
		exit;
	}
while(@ARGV){
	$commandlineoption=shift @ARGV;
		if ($commandlineoption =~ /HELP/i){
			print $helpmenu;
			exit;
		}
	@commandlineoption=split ("=",$commandlineoption);
		if ($commandlineoption[0] =~ /GCC/i){
			if (($commandlineoption[1] >=50) and ($commandlineoption[1] <=70)){
				$GCC=$commandlineoption[1];
			}else{
				$GCC=55;
			}
		}elsif($commandlineoption[0] =~ /OE/i){
			if (($commandlineoption[1] >=0.60) and ($commandlineoption[1] <=1.0)){
				$OE=$commandlineoption[1];
			}else{
				$OE=0.65;
			}
		}elsif($commandlineoption[0] =~ /LENGTH/i){
			if (($commandlineoption[1] >=200) and ($commandlineoption[1] <=1500)){
				$LENGTH=$commandlineoption[1];
			}else{
				$LENGTH=500;
			}
		}
}
print $credit;
&extract($filename,$GCC,$OE,$LENGTH);

sub extract($$$$) {
	$sequence=undef;
	open (IN,$_[0]) || die "$_[0]: $!";
		while(<IN>){
			chomp;
			s/\r\n//g; # remove return
			s/\r//g;
			s/\n//g;
			s/\s//g;
			s/\d//g;
			$sequence=$sequence.$_;
		}
	close(IN);

	@filename=split ('\.',$_[0]);
	$exportfile=$filename[0].'.csv';

open(OUT, ">$exportfile") || die "$exportfile: $!";


$resolution=$_[3];													#definition of resolution
$window_cg=7;
$window_cgr=$_[1];
$window_cgs=$_[2];
my($criteria_cgr)=$_[1];
my($criteria_cgs)=$_[2];
my($criteria_length)=$_[3];


$seqlength=length $sequence;
$CpGno=1;
$CpGstart=0;														#CpGstart+1 is equal to the 
																	#start position of sequences.

	while($CpGstart<=$seqlength-$resolution){						#determine 5' border of CpG island
		$sample1=substr $sequence, $CpGstart, $resolution;			#get a sample
		($cgr,$cgs,$cg)=&get_parameter($sample1);
#		print "$CpGstart\n";

			if (($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg)){
				
				$start=$CpGstart+1;
				$faultstart=$start;									#temporal 5' border
				$faultend=$start+$resolution-1;						#temporal 3' border
					if ($CpGstart<=$seqlength-$resolution*2){
						$CpGstart=$CpGstart+$resolution;			#shift the window to next portion
					}else{
						$CpGstart=$seqlength-$resolution;
					}
enddefinition:
					while($CpGstart<=$seqlength-$resolution){
						$sample2=substr $sequence, $CpGstart, $resolution;
						($cgr,$cgs,$cg)=&get_parameter($sample2);
							if( (($cgr < $window_cgr) or ($cgs < $window_cgs) or ($cg < $window_cg)) or (($CpGstart == $seqlength-$resolution) and ($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg))){
								$tempCpGstart=$CpGstart;			#when the window does not meet the criteria
																	#then shift the window back to 5' side
																	#until the window meets the cirteria
																	
									while ($tempCpGstart >= $CpGstart-$resolution){
										$sample3=substr $sequence, $tempCpGstart, $resolution;	
										($cgr,$cgs,$cg)=&get_parameter($sample3);
											if (($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg)){
												$end=$tempCpGstart+$resolution;
												$cpglength=$end-$start+1;
												$flag=0;
												
													while($cpglength>=$resolution){
														$sample4=substr $sequence, $start-1, $cpglength;
														($cgr,$cgs,$cg)=&get_parameter($sample4);

															if (($cgr>=$criteria_cgr) and ($cgs>=$criteria_cgs) and ($cpglength>=$criteria_length)){
																$cgr=(int($cgr*10))/10;
																$cgs=(int($cgs*1000))/1000;
															#	print "<FONT SIZE=-1>CpG island $CpGno start=$start,end=$end,\%GC=$cgr,ObsCpG/ExpCpG=$cgs,Length=$cpglength</FONT><BR>\n";
																print OUT "$start,$end,$cgr,$cgs,$cpglength\n";
																$CpGstart=$end-1;
																$CpGno++;
																last enddefinition;
															}elsif (($end<=$faultend) or ($cpglength<=$resolution)){
																$cpglength=$faultend-$faultstart+1;
																$sample5=substr $sequence, $faultstart-1, $cpglength;
																($cgr,$cgs,$cg)=&get_parameter($sample5);
																$CpGstart=$faultend-1;
																$cgr=(int($cgr*10))/10;
																$cgs=(int($cgs*1000))/1000;
												#					if ($cpglength<$criteria_length){
												#							last enddefinition;
												#						}elsif(($cgr<$criteria_cgr) or ($cgs<$criteria_cgs)){
												#							$faultstart=$faultstart+$resolution/2;
												#							$faultend=$faultend-$resolution/2;
												#						}else{
																			print OUT "$faultstart,$faultend,$cgr,$cgs,$cpglength\n";
																			$CpGno++;
																			last enddefinition;
												#						}
															}elsif ($flag==0) {
																$end--;		#decrement 3' border
																$cpglength=$end-$start+1;
																$flag=1;
															}elsif ($flag==1) {
																$start++;	#increment 5' border
																$cpglength=$end-$start+1;
																$flag=0;
															}
													}
											}
									}continue{
										$tempCpGstart--;
									}
							}
							
							
							
							
					}continue{
						if ($CpGstart <= $seqlength-$resolution*2){
							$CpGstart=$CpGstart+$resolution;
							$faultend=$faultend+$resolution;
						}else{
							$CpGstart=$seqlength-$resolution;
							$faultend=$seqlength;
						}
					}
			}
	}continue{
		$CpGstart=$CpGstart+1;
	}
			

close(OUT);

open(IN,$exportfile)|| die "$exportfile: $!";
$j=0;
	while(<IN>){
		chomp;
		@temp=split(",",$_);
			for ($k=0;$k<=4;$k++){
				$cpginf[$j]->[$k]=$temp[$k];
			}
		$j++;
	}
close(IN);
	for ($l=1;$l<=$j-1;$l++){
		if ($cpginf[$l]->[0]-$cpginf[$l-1]->[1]<=100){
			$length=$cpginf[$l]->[1]-$cpginf[$l-1]->[0]+1;
			$samplex=substr $sequence,$cpginf[$l-1]->[0]-1,$length;
			($cgr,$cgs,$cg)=&get_parameter($samplex);
				if (($cgr>=$criteria_cgr) and ($cgs>=$criteria_cgs) and ($length>=$criteria_length)){
	
					$cpginf[$l-1]->[1]=$cpginf[$l]->[1];
					$cpginf[$l-1]->[2]=(int($cgr*10))/10;
					$cpginf[$l-1]->[3]=(int($cgs*1000))/1000;
					$cpginf[$l-1]->[4]=$length;
					
						for ($p=$l;$p<=$j-1;$p++){
							for ($q=0;$q<=4;$q++){
								$cpginf[$p]->[$q]=$cpginf[$p+1]->[$q];
							}
						}
						
						for ($q=0;$q<=4;$q++){
							$cpginf[$j-1]->[$q]="";
						}
					$j--;
					$l--;
				}
			
		}
	}

print "Selected lower limits: \%GC=$_[1], ObsCpG/ExpCpG=$_[2], Length=$_[3]\n";
$contigname=substr $exportfile, 0, 8;
open(OUT, '>', $exportfile) or die "$exportfile: $!";
$q=1;
	for ($m=1;$m<=$j;$m++){
			if ($cpginf[$m-1]->[4]>=$criteria_length){
				print OUT "$cpginf[$m-1]->[0],$cpginf[$m-1]->[1]\n";
#				print OUT "CpG island "."$q,$cpginf[$m-1]->[0],$cpginf[$m-1]->[1],$cpginf[$m-1]->[2],$cpginf[$m-1]->[3],$cpginf[$m-1]->[4]\n";
#				print "$contigname, CpG island $q, start=$cpginf[$m-1]->[0], end=$cpginf[$m-1]->[1], \%GC=$cpginf[$m-1]->[2], ObsCpG/ExpCpG=$cpginf[$m-1]->[3], Length=$cpginf[$m-1]->[4]\n";
#				print "$contigname, CpG island $q, $cpginf[$m-1]->[0], $cpginf[$m-1]->[1], $cpginf[$m-1]->[2], $cpginf[$m-1]->[3], $cpginf[$m-1]->[4]\n";
				$q++;
			}
	}
if ($q==1){
	print "No CpG islands was found.\n";
}

#close(OUT);
#unlink "$exportfile";

sub get_parameter ($){
	my $cgs;
	my($sequence)=@_[0];
	my($seqlength)=length $sequence;
	my($c)=($sequence =~s/c/c/ig);
	my($g)=($sequence =~s/g/g/ig);
	my($cg)=($sequence =~s/cg/cg/ig);
	my($cgr)=($c+$g)/$seqlength*100;
		if (($c==0) or ($g==0)) {
			$cgs=0;
		}else{
			$cgs=$cg*$seqlength/$c/$g;
		}
	($cgr,$cgs,$cg)	
}

}

