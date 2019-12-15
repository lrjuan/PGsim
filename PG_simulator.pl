#!/usr/bin/perl -w

use List::Util qw(sum);

if(@ARGV<1 || $ARGV[0] eq '--help' || $ARGV[0] eq '-h'){
	print STDERR "Usage: perl PG_simulator.pl Personal_genome_ID (consistent with PG_planner)\n";
}else{
#Initializing Parameter
$OVR = 0;
$KVR = 0;
$CVR = 0;
$pnvNum = 0;
$SVN = 0;
$PVN = 0;
$PLalpha = 1.8;

#Reading Parameter

open(CONF,"<$ARGV[0].pgsim.allparams.conf") || die "No such file: $ARGV[0].pgsim.allparams.conf";
while(<CONF>){
	chomp;
	@temp=split(/=/);
	${$temp[0]}=$temp[1];
}
close(CONF);
print STDERR "Loading reference genome...\n";

if($Genome=~/gz$/){
	open(REF,"gzip -dc $Genome|") || die "No such file: $Genome\n";
}else{
	open(REF,"<$Genome") || die "No such file: $Genome\n";
}
my %ref;
our @chrs;
$chr='chr';
while(<REF>){
	chomp;
	if (/^>/){
		$chrs[$#chrs][1]=length($ref{$chr}) if($chr ne 'chr');
		$chr=XYM_trans(substr $_,1);
		push(@chrs,[$chr,0]);
	}
	else{
		$ref{$chr}.=$_;
	}
}
$chrs[$#chrs][1]=length($ref{$chr});
close(REF);

print STDERR "Simulating personal genome...\n";

our %NV;
our %SV;
our %PV;

if($pnvNum+$CDN>0){
	open(NVAR,"<$ARGV[0].nv.loc") || die "No such file: $ARGV[0].nv.loc";
	while(<NVAR>){
		chomp;
		@temp = split(/\t/);
		$temp[0] = XYM_trans($temp[0]);
		$NV{$temp[0]} = [] if(!exists($NV{$temp[0]}));
		push(@{$NV{$temp[0]}},[$temp[1],$temp[2]]);
	}
	close(NVAR);
	foreach $chr (@chrs){
		@{$NV{$chr->[0]}}=sort {$a->[0] <=> $b->[0]} @{$NV{$chr->[0]}} if(exists($NV{$chr->[0]}));
	}
	delete $NV{24} if ($Gender eq 'Female' && exists($NV{24}));
}
if($SVN>0){
	open(SVAR,"<$ARGV[0].$svDB.pgsim.vcf") || die "No such file: $ARGV[0].$svDB.pgsim.vcf";
	while(<SVAR>){
		next if(/^#/);
		chomp;
		@temp = split(/\t/);
		$temp[0] = XYM_trans($temp[0]);
		@temp7 = split(/;/,$temp[7]);
		my $end = -1;
		foreach $info (@temp7) {
			if($info=~/^END=/){
				$end = substr($info,4);
				last;
			}
		}
		$end = $end>$temp[1]?$end:$temp[1];
		$SV{$temp[0]} = [] if(!exists($SV{$temp[0]}));
		push(@{$SV{$temp[0]}},[$temp[1],$end,$_]);
	}
	close(SVAR);
	foreach $chr (@chrs){
		@{$SV{$chr->[0]}}=sort {$a->[0] <=> $b->[0]} @{$SV{$chr->[0]}} if(exists($SV{$chr->[0]}));
	}
	delete $SV{24} if ($Gender eq 'Female' && exists($SV{24}));
}
if($PVN>0){
	open(PVAR,"<$ARGV[0].$pvDB.pgsim.vcf") || die "No such file: $ARGV[0].$pvDB.pgsim.vcf";
	while(<PVAR>){
		next if(/^#/);
		chomp;
		@temp = split(/\t/);
		$temp[0] = XYM_trans($temp[0]);
		$PV{$temp[0]} = [] if(!exists($PV{$temp[0]}));
		push(@{$PV{$temp[0]}},[$temp[1],$temp[1]+length($temp[3])-1,$_]);
	}
	close(PVAR);
	foreach $chr (@chrs){
		@{$PV{$chr->[0]}}=sort {$a->[0] <=> $b->[0]} @{$PV{$chr->[0]}} if(exists($PV{$chr->[0]}));
	}
	delete $PV{24} if ($Gender eq 'Female' && exists($PV{24}));
}

open(VCFOUT,">$ARGV[0].vcf") || die;

print VCFOUT "##fileformat=VCFv4.1\n###source=Personal Genome simulator (PGsim)\n###reference=$Genome\n";
print VCFOUT "###Overall variation rate=$OVR\n";
print VCFOUT "###Known variation rate=$KVR\n";
print VCFOUT "###Common variation rate=$CVR\n";
print VCFOUT "###Common variation AF threshold=$CVT\n";
print VCFOUT "###Novel indel rate=$NIDR\n";
print VCFOUT "###Novel indel length limit=$NIDL\n";
print VCFOUT "###Number of novel variants in CDS region=$CDN\n";
print VCFOUT "###Frame shifting ratio=$FSR\n";
print VCFOUT "###Transition/Transversion Ratio=$TiTv\n";
print VCFOUT "###Transition/Transversion Ratio in Coding Region=$TiTvC\n";
print VCFOUT "##ALT=<ID=INS,Description=\"Insertion\">\n";
print VCFOUT "##ALT=<ID=INV,Description=\"Invertion\">\n";
print VCFOUT "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print VCFOUT "##ALT=<ID=DUP,Description=\"Duplication\">\n";
print VCFOUT "##ALT=<ID=CNV,Description=\"Copy number variation\">\n";
print VCFOUT "##INFO=<ID=COMVAR,Number=0,Type=Flag,Description=\"Common Variant\">\n";
print VCFOUT "##INFO=<ID=NCMVAR,Number=0,Type=Flag,Description=\"Non-common Variant\">\n";
print VCFOUT "##INFO=<ID=NCNV,Number=0,Type=Flag,Description=\"Novel Non-coding Variant\">\n";
print VCFOUT "##INFO=<ID=CDSNV,Number=0,Type=Flag,Description=\"Novel Coding Variant\">\n";
print VCFOUT "##INFO=<ID=ALLELES,Type=String,Description=\"Original alleles in the databases, corresponding to AF\">\n";
print VCFOUT "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise SV\">\n";
print VCFOUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"SV end\">\n";
print VCFOUT "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n";
print VCFOUT "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">\n";
print VCFOUT "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval of POS\">\n";
print VCFOUT "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval of END\">\n";
print VCFOUT "##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval of SVLEN\">\n";
print VCFOUT "##FORMAT=<ID=GT,Number=1,Type=GT,Description=\"Genotype\">\n";
print VCFOUT "##SAMPLE=<ID=$ARGV[0],Description=\"Individual, $Gender, $Population\">\n";
print VCFOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ARGV[0]\n";

if($commonDB=~/gz$/){
	open(COM,"gzip -dc $commonDB|") || die "No such file: $commonDB\n";
}else{
	open(COM,"<$commonDB") || die "No such file: $commonDB\n";
}
if($allvarDB=~/gz$/){
	open(VAR,"gzip -dc $allvarDB|") || die "No such file: $allvarDB\n";
}else{
	open(VAR,"<$allvarDB") || die "No such file: $allvarDB\n";
}

our @pointers;
$pointers[0]=<COM>;		#cmpointer
$pointers[1]=<VAR>;		#dbpointer
$pointers[2]=0;			#nvpointer
$pointers[3]=0;			#pvpointer
$pointers[4]=0;			#svpointer
while($pointers[0]=~/^##/){$pointers[0]=<COM>;}		
while($pointers[1]=~/^##/){$pointers[1]=<VAR>;}	

our $currentCHR = -1;
our $currentPOS = 0;

our @common_cache;
$common_cache[0]=-1;
$common_cache[1]=0;
our @allvar_cache;
$allvar_cache[0]=-1;
$allvar_cache[1]=0;

my @vcfline;
while($pointers[0]||$pointers[1]||scalar(keys %NV)||scalar(keys %SV)||scalar(keys %PV)){
	$vcfline[0]='chr0';
	$vcfline[1]='0';
	$vcfline[2]='.';
	$vcfline[3]='N';
	$vcfline[4]='.';
	$vcfline[5]='.';
	$vcfline[6]='.';
	$vcfline[7]='.';
	$vcfline[8]='GT';
	$vcfline[9]='./.';
	$p = coord_regist();

	if($p==4){
		@temp=split(/\t/,$SV{$currentCHR}[$pointers[4]][2]);

		$currentPOS = $SV{$currentCHR}[$pointers[4]][1];

		my $temp_vcf9="";
		if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male'){
			$temp_vcf9='1';
		}elsif(rand()<0.5){
			$temp_vcf9='0|1';
		}else{
			$temp_vcf9='1|0';
		}

		$temp[0] = XYM_trans_rev($temp[0]);
		print VCFOUT join("\t",@temp)."\tGT\t$temp_vcf9\n";

		$vcfline[0]='chr0';
	}
	elsif($p==3){
		@temp=split(/\t/,$PV{$currentCHR}[$pointers[3]][2]);
		
		$currentPOS = $PV{$currentCHR}[$pointers[3]][1];

		my $temp_vcf9="";
		if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male'){
			$temp_vcf9='1';
		}elsif(rand()<0.5){
			$temp_vcf9='0|1';
		}else{
			$temp_vcf9='1|0';
		}

		$temp[0] = XYM_trans_rev($temp[0]);
		print VCFOUT join("\t",@temp)."\tGT\t$temp_vcf9\n";

		$vcfline[0]='chr0';
	}
	elsif($p==0){
		chomp($pointers[0]);
		@temp=split(/\t/,$pointers[0]);

		$currentPOS = $common_cache[1]+length($temp[3])-1;

		$allele[0]='';
		$allele[1]='';
		my @afs=@{$common_cache[2]};
		my @alt_temp=@{$common_cache[3]};

		$sampling1 = rand();
		$sampling2 = rand();
		for($i=@afs-@alt_temp;$i<@afs;$i++){
			if($sampling1>=0 && $afs[$i] ne '.'){
				$sampling1-=$afs[$i]*$paC; 
				$allele[0]=$alt_temp[$i-@afs+@alt_temp] if($sampling1<0);
			}
			if($sampling2>=0 && $afs[$i] ne '.'){
				$sampling2-=$afs[$i]*$paC;
				$allele[1]=$alt_temp[$i-@afs+@alt_temp] if($sampling2<0);
			}
		}

		$allele[1] = '' if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male');

		if ("$allele[0]$allele[1]" ne ''){
			$vcfline[0]=XYM_trans_rev($temp[0]);
			$vcfline[1]=$temp[1];
			$vcfline[2]=$temp[2];
			$vcfline[3]=$temp[3];
			$vcfline[7]='COMVAR;ALLELES='.join(',',@alt_temp).";".$temp[7];
			if ($allele[0] eq ''){
				$vcfline[4]=$allele[1];
				$vcfline[9]='0|1';
			}
			elsif ($allele[1] eq ''){
				$vcfline[4]=$allele[0];
				$vcfline[9]='1|0';
				$vcfline[9]='1' if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male');
			}
			elsif ($allele[0] eq $allele[1]){
				$vcfline[4]=$allele[0];
				$vcfline[9]='1|1';
			}
			else {
				$vcfline[4]=$allele[0].','.$allele[1];
				$vcfline[9]='1|2';
			}
		}
		else{
			$vcfline[0]='chr0';
		}
	}
	elsif($p==1){
		chomp($pointers[1]);
		@temp=split(/\t/,$pointers[1]);
		my @alt_temp=split(/,/,$temp[4]);

		$currentPOS = $allvar_cache[1]+length($temp[3])-1;

		$allele[0]='';
		$allele[1]='';

		$allele[0]=$alt_temp[int(rand(@alt_temp))] if(rand()<$paAF);
		$allele[1]=$alt_temp[int(rand(@alt_temp))] if(rand()<$paAF);

		$allele[1] = '' if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male');

		if ("$allele[0]$allele[1]" ne ''){
			$vcfline[0]=XYM_trans_rev($temp[0]);
			$vcfline[1]=$temp[1];
			$vcfline[2]=$temp[2];
			$vcfline[3]=$temp[3];
			$vcfline[7]='NCMVAR;ALLELES='.join(',',@alt_temp).";".$temp[7];
			if ($allele[0] eq ''){
				$vcfline[4]=$allele[1];
				$vcfline[9]='0|1';
			}
			elsif ($allele[1] eq ''){
				$vcfline[4]=$allele[0];
				$vcfline[9]='1|0';
				$vcfline[9]='1' if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male');
			}
			elsif ($allele[0] eq $allele[1]){
				$vcfline[4]=$allele[0];
				$vcfline[9]='1|1';
			}
			else {
				$vcfline[4]=$allele[0].','.$allele[1];
				$vcfline[9]='1|2';
			}
		}
		else{
			$vcfline[0]='chr0';
		}
	}
	elsif($p==2){
		$vcfline[0]=XYM_trans_rev($currentCHR);
		$vcfline[1]=$NV{$currentCHR}[$pointers[2]][0];
		$vcfline[2]='.';
		$vcfline[7]=$NV{$currentCHR}[$pointers[2]][1];
		if(rand()<0.5){
			$vcfline[9]='0|1';
		}else{
			$vcfline[9]='1|0';
		}
		$vcfline[9]='1' if($currentCHR == 25 || ($currentCHR == 23 || $currentCHR == 24) && $Gender eq 'Male');

		if(rand()<$NIDR){#indel
			$indel_len=int(rand()**(1/(1-$PLalpha)));
			$indel_len=$indel_len%$NIDL+1 if($indel_len>$NIDL);
			if($indel_len%3!=0 && rand()>$FSR && $vcfline[7] eq 'CDSNV'){
				$vcfline[0]='chr0';
			}else{
				if(rand()<0.5){#deletion
					$vcfline[3]=substr($ref{$currentCHR},$vcfline[1]-1,$indel_len+1);
					$vcfline[4]=substr($vcfline[3],0,1);
				}else{#insertion
					$vcfline[3]=substr($ref{$currentCHR},$vcfline[1]-1,1);
					$vcfline[4]='';
					for($i=0;$i<$indel_len;$i++){
						$vcfline[4].=int(rand(10));
					}
					$vcfline[4]=~s/0/A/g;
					$vcfline[4]=~s/1/A/g;
					$vcfline[4]=~s/2/A/g;
					$vcfline[4]=~s/3/T/g;
					$vcfline[4]=~s/4/T/g;
					$vcfline[4]=~s/5/T/g;
					$vcfline[4]=~s/6/C/g;
					$vcfline[4]=~s/7/C/g;
					$vcfline[4]=~s/8/G/g;
					$vcfline[4]=~s/9/G/g;
					$vcfline[4]=$vcfline[3].$vcfline[4];
				}
			}
		}else{#SNV
			$vcfline[3]=substr($ref{$currentCHR},$vcfline[1]-1,1);
			$vcfline[3]=uc($vcfline[3]);
			$iv = 'NNN';
			if($vcfline[3] eq 'A'){
				$iv='TCG';
			}
			elsif($vcfline[3] eq 'T'){
				$iv='AGC';
			}
			elsif($vcfline[3] eq 'C'){
				$iv='GAT';
			}
			elsif($vcfline[3] eq 'G'){
				$iv='CTA';
			}
			else {
				$vcfline[0]='chr0';
			}
			if(rand($TiTvC+1)<$TiTvC && $vcfline[7] eq 'CDSNV' || rand($TiTv+1)<$TiTv && $vcfline[7] eq 'NCNV'){
				$vcfline[4]=substr($iv,0,1);
			}else{
				$vcfline[4]=substr($iv,int(rand(2))+1,1);
			}
		}
		$vcfline[0]='chr0' if ($vcfline[3]=~/N/);
		$currentPOS = $NV{$currentCHR}[$pointers[2]][0]+length($vcfline[3])-1;
	}
	if($vcfline[0] ne 'chr0'){
		$vcfline[3]=uc($vcfline[3]);
		$vcfline[4]=uc($vcfline[4]);
		print VCFOUT join("\t",@vcfline)."\n";
	}
}

close(COM);
close(VAR);
close(VCFOUT);

sub coord_regist {
	my $pointer = -1;
	my $minCHR = 26;
	my $minPOS = 2500000000;

	if($common_cache[0]==$currentCHR && $common_cache[1]<=$currentPOS){
		$pointers[0]=<COM>;
		while($pointers[0]){
			@temp = split(/\t/,$pointers[0]);
			$temp[0]=XYM_trans($temp[0]);
			if(($temp[0] == $currentCHR && $temp[1]<=$currentPOS) || ($temp[0]==24 && $Gender eq 'Female')){
				$pointers[0]=<COM>;
			}else{
				my @info_temp=split(/;/,$temp[7]);
				my @alt_temp=split(/,/,$temp[4]);
				my @afs;
				foreach $info (@info_temp){
					if($info=~/^$Population=/){
						$info=~s/,\./,0/g;
						@afs = split(/[=,]/,$info);
						last;
					}
				}
				if(@afs>0 && sum(@afs[(@afs-@alt_temp)..$#afs]) >= $CVT){
					$common_cache[0]=$temp[0];
					$common_cache[1]=$temp[1];
					delete $common_cache[2];
					delete $common_cache[3];
					$common_cache[2]=\@afs;
					$common_cache[3]=\@alt_temp;
					last;
				}else{
					$pointers[0]=<COM>;
				}
			}
		}	
	}
	if($pointers[0]){
		$minCHR = $common_cache[0];
		$minPOS = $common_cache[1];
		$pointer = 0;
	}else{
		$common_cache[0]=26;
		$common_cache[1]=250000000;
	}
	
	if($allvar_cache[0]==$currentCHR && $allvar_cache[1]<=$currentPOS){
		$pointers[1]=<VAR>;
		while($pointers[1]){
			@temp = split(/\t/,$pointers[1]);
			$temp[0]=XYM_trans($temp[0]);
			if($temp[0] == $currentCHR && $temp[1]<=$currentPOS || $temp[0] == 24 && $Gender eq "Female"){
				$pointers[1]=<VAR>;
			}else{
				$allvar_cache[0]=$temp[0];
				$allvar_cache[1]=$temp[1];
				last;
			}
		}
	}
	if($pointers[1]){
		if($minCHR == $currentCHR && $allvar_cache[0] == $currentCHR && $minPOS > $allvar_cache[1] ||
		   $minCHR != $currentCHR && $allvar_cache[0] == $currentCHR ||
		   $minCHR != $currentCHR && $allvar_cache[0] != $currentCHR &&($minPOS > $allvar_cache[0] || 
	   																	 $minCHR ==$allvar_cache[0] && 
	   																	 $minPOS > $allvar_cache[1])){
			$minCHR = $allvar_cache[0];
			$minPOS = $allvar_cache[1];
			$pointer = 1;
		}
	}else{
		$allvar_cache[0]=26;
		$allvar_cache[1]=250000000;
	}

	while(exists($NV{$currentCHR}) && $pointers[2]<@{$NV{$currentCHR}}){
		if($NV{$currentCHR}[$pointers[2]][0]<=$currentPOS){
			$pointers[2]++;
		}else{
			if($minCHR != $currentCHR || ($minCHR == $currentCHR && $minPOS > $NV{$currentCHR}[$pointers[2]][0])){
				$minCHR = $currentCHR;
				$minPOS = $NV{$currentCHR}[$pointers[2]][0];
				$pointer = 2;
			}
			last;
		}
	}

	while(exists($PV{$currentCHR}) && $pointers[3]<@{$PV{$currentCHR}}){
		if($PV{$currentCHR}[$pointers[3]][0]<=$currentPOS){
			$pointers[3]++;
		}else{
			if($minCHR != $currentCHR || ($minCHR == $currentCHR && $minPOS >= $PV{$currentCHR}[$pointers[3]][0])){
				$minCHR = $currentCHR;
				$minPOS = $PV{$currentCHR}[$pointers[3]][0];
				$pointer = 3;
			}
			last;
		}
	}

	while(exists($SV{$currentCHR}) && $pointers[4]<@{$SV{$currentCHR}}){
		if($SV{$currentCHR}[$pointers[4]][0]<=$currentPOS){
			$pointers[4]++;
		}else{
			if($minCHR != $currentCHR || ($minCHR == $currentCHR && $minPOS >= $SV{$currentCHR}[$pointers[4]][0])){
				$minCHR = $currentCHR;
				$minPOS = $SV{$currentCHR}[$pointers[4]][0];
				$pointer = 4;
			}
			last;
		}
	}

	if($minCHR != $currentCHR){
		delete $NV{$currentCHR} if(exists($NV{$currentCHR}));
		delete $PV{$currentCHR} if(exists($PV{$currentCHR}));
		delete $SV{$currentCHR} if(exists($SV{$currentCHR}));
		print STDERR XYM_trans_rev($currentCHR)." completed.\n" if($currentCHR != -1);

		if($minCHR == 26){
			if(scalar(keys %NV) > 0){
				my @key_temp = keys %NV;
				$minCHR = $key_temp[0];
			}elsif(scalar(keys %PV) > 0){
				my @key_temp = keys %PV;
				$minCHR = $key_temp[0];	
			}elsif(scalar(keys %SV) > 0){
				my @key_temp = keys %SV;
				$minCHR = $key_temp[0];				
			}
		}

		$currentCHR = $minCHR;
		$currentPOS = 0;
		$pointers[2]=0;
		$pointers[3]=0;
		$pointers[4]=0;
		$pointer = -1;
	}
	return $pointer;
}

sub XYM_trans {
	$chr=shift;
	if (($chr eq "X") || ($chr eq "chrX")){
		return 23;
	}
	elsif (($chr eq "Y") || ($chr eq "chrY")){
		return 24;
	}
	elsif (($chr eq "MT") || ($chr eq "chrM") || ($chr eq "M")){
		return 25;
	}
	elsif ($chr=~/^chr/){
		return substr($chr,3);
	}
	elsif ($chr=~/^\d+$/ && $chr >= 1 && $chr <= 25){
		return $chr;
	}
	else{
		return -1;
	}
}

sub XYM_trans_rev {
	$chr=shift;
	if ($chr eq 'X' || $chr eq '23'){
		return 'chrX';
	}
	elsif ($chr eq 'Y' || $chr eq '24'){
		return 'chrY';
	}
	elsif ($chr eq 'MT' || $chr eq 'M' || $chr eq '25'){
		return 'chrM';
	}
	elsif ($chr==-1){
		return 'chr0';
	}
	elsif ($chr=~/^chr/){
		return $chr;
	}
	else {
		return 'chr'.$chr;
	}
}

}__END__