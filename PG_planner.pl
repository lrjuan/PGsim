#!/usr/bin/perl -w
use List::Util qw(sum);

if(@ARGV<2){
	print STDERR "\nUsage:perl PG_simulator.pl user_configuration_file Personal_Genome_ID\n";
}

#Initializing parameters
$Gender="Any";
$Population="CAF";
$Genome="GRCh38.fa.gz";
$GModel="hg38.cds.bed";	
$allvarDB="All_20180418.vcf.gz";
$commonDB="common_all_20180418.vcf.gz";
$svDB="GRCh38.variant_call.vcf.gz";
$pvDB="clinvar_20191001.vcf.gz";

$OVR=0.001;
$KVR=0.9;
$CVR=0.8;
$SVN=0;
$SV_LIM=1000000;
$PVN=0;

$NIDR=0.1;
$NIDL=50;
$CDN=200;
$FSR=0.3;
$TiTv=2;
$TiTvC=2.8;
$PLalpha=1.8;

#Reading Parameter
print STDERR "Initializing parameters...\n";

open(CONF,"$ARGV[0]") || die;
while(<CONF>){
	chomp;
	next if /^#/ || /^\s/;
	s/\t+#.*$//;
	@temp=split(/[=,]/);
	if (@temp<3 || $temp[2]==$temp[1]) {
		${$temp[0]}=$temp[1];
	}
	else {
		${$temp[0]}=$temp[1]+rand($temp[2]-$temp[1]);
		${$temp[0]}=int(${$temp[0]}) unless ("$temp[1]$temp[2]"=~/\./);
	}
}
close(CONF);

if($Gender eq 'Female' || $Gender eq 'female' || $Gender eq 'f' || $Gender eq 'F'){
	$Gender = 'Female';
}elsif($Gender eq 'Male' || $Gender eq 'male' || $Gender eq 'm' || $Gender eq 'M'){
	$Gender = 'Male';
}else{
	if(rand()>0.5){
		$Gender = "Female";
	}else{
		$Gender = "Male";
	}
}

#Analyzing Database, if the analysis had been done and the '.pgstat' files exist, this step will be skipped

print STDERR "Analyzing databases...\n";

my @chrs;
my %chrs;
my $GLN = 0;
my $nGLN = 0;
if(!(-e $Genome)){
	print STDERR "Invalid reference genome!\n";
}elsif(-e "$Genome.pgstat"){
	open(REFSTAT,"<$Genome.pgstat");
	while(<REFSTAT>){
		next if (/^#/);
		chomp;
		if(/^All/){
			@temp=split(/\t/);
			$GLN=$temp[1];
			$nGLN=$temp[2];
		}else{
			push(@chrs,[split(/\t/)]);
			$chrs{$chrs[$#chrs][0]}=$#chrs;
		}
	}
	close(REFSTAT);
}else{
	if($Genome=~/gz$/){
		open(REF,"gzip -dc $Genome|") || die "No such file or directory $Genome\n";
	}else{
		open(REF,"<$Genome") || die "No such file or directory $Genome\n";
	}
	while(<REF>){
		chomp;
		if (/^>/){
			$chr=XYM_trans(substr $_,1);
			push(@chrs,[$chr,0,0]);
			$chrs{$chr}=$#chrs;
		}
		else{
			$chrs[$#chrs][1]+=length;
			$chrs[$#chrs][2]+=tr/N/N/;
		}
	}
	close(REF);
	open(REFSTAT,">$Genome.pgstat") || die;
	print REFSTAT "#Chr\tLength\tNumberOfN\n";
	for($i=0;$i<@chrs;$i++){
		print REFSTAT join("\t",@{$chrs[$i]})."\n";
		$GLN+=$chrs[$i][1];
		$nGLN+=$chrs[$i][2];
	}
	print REFSTAT "All\t$GLN\t$nGLN\n";
	close(REFSTAT);
}
print STDERR "Reference genome analysis done.\n";

my $modLen = 0;
my %mod;

if(!(-e $GModel)){
	print STDERR "Invalid gene model!\n";
}else{
	open(MOD,"<$GModel") || die;
	while(<MOD>){
		chomp;
		@temp=split(/\t/);
		$temp[0]=XYM_trans($temp[0]);
		if(exists($chrs{$temp[0]})){
			$mod{$temp[0]} = [] if(!exists($mod{$temp[0]}));
			push(@{$mod{$temp[0]}},[$temp[1],$temp[2]]);
			$modLen+=$temp[2]-$temp[1];
		}
	}
	close(MOD);
	foreach $chr (@chrs){
		next if(!exists($mod{$chr->[0]}));
		@{$mod{$chr->[0]}}=sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$mod{$chr->[0]}};
	}
}
print STDERR "Gene model analysis done.\n";

my $varNum = 0;

if(!(-e $allvarDB)){
	print STDERR "Invalid known variants database!\n";
}elsif(-e "$allvarDB.pgstat"){
	open(VARSTAT,"<$allvarDB.pgstat") || die;
	while(<VARSTAT>){
		next if (/^#/);
		chomp;
		@temp=split(/\t/);
		$varNum=$temp[1];
	}
	close(VARSTAT);
}else{
	if($allvarDB=~/gz$/){
		open(VAR,"gzip -dc $allvarDB|") || die;
	}
	else{
		open(VAR,"<$allvarDB") || die;
	}
	my %vars;
	while(<VAR>){
		next if(/^#/);
		@temp=split(/\t/);
		$temp[0]=XYM_trans($temp[0]);
		if(exists($chrs{$temp[0]})){
			if(exists($vars{$temp[0]})){
				$vars{$temp[0]}++;
			}else{
				$vars{$temp[0]}=1;
			}
			$varNum++;
		}
	}
	close(VAR);
	open(VARSTAT,">$allvarDB.pgstat") || die;
	print VARSTAT "#Chr\tVarNum\n";
	foreach $chr (@chrs){
		if(exists($vars{$chr->[0]})){
			print VARSTAT "$chr->[0]\t$vars{$chr->[0]}\n";
		}else{
			print VARSTAT "$chr->[0]\t0\n";
		}
	}
	print VARSTAT "All\t$varNum\n";
	close(VARSTAT);
}
print STDERR "All known variants DB analysis done.\n";

my $cvNum = 0;
my $ecvNum = 0;
my $CVT = 0.01;

if(!(-e $commonDB)){
	print STDERR "Invalid common variants database!\n";
}elsif(-e "$commonDB.$Population.pgstat"){
	open(COMSTAT,"<$commonDB.$Population.pgstat") || die;
	while(<COMSTAT>){
		next if (/^#/);
		chomp;
		@temp=split(/\t/);
		$cvNum=$temp[1];
		$ecvNum=$temp[2];
	}
	close(COMSTAT);
}else{
	if($commonDB=~/gz$/){
		open(COM,"gzip -dc $commonDB|") || die;
	}
	else{
		open(COM,"<$commonDB") || die;
	}
	my %cvars;
	foreach $chr (@chrs){
		$cvars{$chr->[0]}=[0,0];
	}
	$comnum=0;
	while(<COM>){
		$comnum++;
		next if(/^#/);
		chomp;
		@temp=split(/\t/);
		$temp[0]=XYM_trans($temp[0]);
		if(exists($chrs{$temp[0]})){
			@temp7=split(/;/,$temp[7]);
			@temp4=split(/,/,$temp[4]);
			foreach $info (@temp7) {
				if($info=~/^$Population=/){
					$info=~s/,\./,0/g;
					@afs = split(/[=,]/,$info);
					$mafs = sum(@afs[(@afs-@temp4)..$#afs]);
					if($mafs>=$CVT){
						$cvars{$temp[0]}->[0]++;
						$cvars{$temp[0]}->[1]+=$mafs;
					}
					last;
				}
			}
		}
	}
	close(COM);

	open(COMSTAT,">$commonDB.$Population.pgstat") || die;
	print COMSTAT "#Chr\tvarNum\texpectedVarNum\n";
	foreach $chr (@chrs){
		print COMSTAT "$chr->[0]\t$cvars{$chr->[0]}->[0]\t$cvars{$chr->[0]}->[1]\n";
		$cvNum+=$cvars{$chr->[0]}->[0];
		$ecvNum+=$cvars{$chr->[0]}->[1];
	}
	print COMSTAT "All\t$cvNum\t$ecvNum\n";
	close(COMSTAT);
}
print STDERR "Common variants DB analysis done.\n";

print STDERR "Calculating other internal parameters...\n";

$eoR = ($OVR*$KVR)/(1-$OVR*$KVR);
$paC = $GLN*$OVR*$KVR*$CVR/$ecvNum;
$paAF = $GLN*$OVR*$KVR*(1-$CVR)/($varNum-$cvNum);
$pnvNum = int($GLN*$OVR*(1-$KVR)*(1+$eoR)*($GLN/($GLN-$nGLN-$modLen)));

#Extracting SVs and PVs

if($SVN>0){
	print STDERR "Extracting structural variations...\n";
	my @svars;
	my %overlen_svars;
	if(!(-e $svDB)){
		print STDERR "Invalid structural variation database!\n";
	}else{
		if($svDB=~/gz$/){
			open(SVAR,"gzip -dc $svDB|");
		}
		else{
			open(SVAR,"<$svDB") || die;
		}
		$i=0;
		while(<SVAR>){
			$i++;
			next if(/^#/);
			chomp;
			@temp=split(/\t/);
			$temp[0]=XYM_trans($temp[0]);
			if(exists($chrs{$temp[0]})){
				@temp7=split(/;/,$temp[7]);
				my $end = -1;
				foreach $info (@temp7) {
					if($info=~/^END=/){
						$end = substr($info,4);
						last;
					}
				}
				$end = $end>$temp[1]?$end:$temp[1];

				if($end-$temp[1]>$SV_LIM){
					$overlen_svars{$i}=1;
					next;
				}

				if(scalar(@svars)==0 || $svars[$#svars][0] != $temp[0] || $svars[$#svars][4]<$temp[1]){
					push(@svars,[$temp[0],$i,$i,$temp[1],$end]);
				}else{
					$svars[$#svars][2]=$i;
					$svars[$#svars][4]=$end if($end > $svars[$#svars][4]);
				}
			}
		}
		close(SVAR);
		my %sv_extract;
		my %sv_extractd_id;
		while(keys %sv_extract < $SVN){
			$tempName=int(rand(@svars));
			$tempId = $svars[$tempName][1]+int(rand($svars[$tempName][2]-$svars[$tempName][1]+1));
			if (!exists($sv_extract{$tempName}) && !exists($overlen_svars{$tempId})){
				$sv_extract{$tempName} = $tempId;
				$sv_extract_id{$tempId} = 1;
			}
		}
		if($svDB=~/gz$/){
			open(SVAR,"gzip -dc $svDB|") || die;
		}
		else{
			open(SVAR,"<$svDB") || die;
		}
		open(SVSIM,">$ARGV[1].$svDB.pgsim.vcf") || die;
		print SVSIM "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
		$i=0;
		while(<SVAR>){
			$i++;
			print SVSIM if(exists($sv_extract_id{$i}));
		}
		close(SVAR);
		close(SVSIM);
	}
}

if($PVN>0){
	print STDERR "Extracting disease-related variants...\n";
	my @pvars;
	if(!(-e $pvDB)){
		print STDERR "Invalid disease-related variants database!\n";
	}else{
		if($pvDB=~/gz$/){
			open(PVAR,"gzip -dc $pvDB|");
		}
		else{
			open(PVAR,"<$pvDB") || die;
		}
		$i=0;
		while(<PVAR>){
			$i++;
			next if(/^#/);
			chomp;
			@temp=split(/\t/);
			$temp[0]=XYM_trans($temp[0]);
			my $sig = "";
			if(exists($chrs{$temp[0]})){
				@temp7=split(/;/,$temp[7]);
				foreach $info (@temp7) {
					if($info=~/^CLNSIG=/){
						$sig = substr($info,7);
						last;
					}
				}
				push(@pvars,$i) if($sig =~ /athogenic/ && $sig !~ /onflict/);
			}
		}
		close(PVAR);
		my %pv_extract;
		my %pv_extract_id;
		while(keys %pv_extract < $PVN){
			$tempName=int(rand(@pvars));
			if(!exists($pv_extract{$tempName})){
				$pv_extract{$tempName}=$pvars[$tempName];
				$pv_extract_id{$pvars[$tempName]}=1;
			}
		}
		if($pvDB=~/gz$/){
			open(PVAR,"gzip -dc $pvDB|");
		}
		else{
			open(PVAR,"<$pvDB") || die;
		}
		open(PVSIM,">$ARGV[1].$pvDB.pgsim.vcf") || die;
		print PVSIM "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
		$i=0;
		while(<PVAR>){
			$i++;
			print PVSIM if(exists($pv_extract_id{$i}));
		}
		close(PVAR);
		close(PVSIM);
	}
}

#Planning Novel Variants

print STDERR "planning novel variants...\n";

my %NovelLoc;
my %CodingLoc;
while(keys %NovelLoc < $pnvNum){
	$tempName=int(rand($GLN));
	$NovelLoc{$tempName}=1 if (!exists($NovelLoc{$tempName}));
}
while(keys %CodingLoc < $CDN){
	$tempName=int(rand($modLen));
	$CodingLoc{$tempName}=1 if (!exists($CodingLoc{$tempName}));
}
my @NovelLoc=sort {$a<=>$b} (keys %NovelLoc);
my @CodingLoc=sort {$a<=>$b} (keys %CodingLoc);

$nvloc = 0;
$cnvloc = 0;
$nvOffset = 0;
$cnvOffset = 0;
open(NVSIM,">$ARGV[1].nv.loc") || die;
foreach $chr (@chrs){
	$modLoc = 0;
	while($nvloc<@NovelLoc || $cnvloc<@CodingLoc){
		if($nvloc<@NovelLoc && (!exists($mod{$chr->[0]}) || $modLoc>=@{$mod{$chr->[0]}}) && $NovelLoc[$nvloc]-$nvOffset<=$chr->[1]){
			print NVSIM "$chr->[0]\t".($NovelLoc[$nvloc]-$nvOffset)."\tNCNV\n";
			$nvloc++;
		}elsif((!exists($mod{$chr->[0]}) || $modLoc>=@{$mod{$chr->[0]}}) && $nvloc<@NovelLoc && $NovelLoc[$nvloc]-$nvOffset>$chr->[1]){
			last;
		}elsif($nvloc<@NovelLoc && $NovelLoc[$nvloc]-$nvOffset<=$chr->[1] && $modLoc<@{$mod{$chr->[0]}} && $NovelLoc[$nvloc]-$nvOffset<=$mod{$chr->[0]}[$modLoc][0]){
			print NVSIM "$chr->[0]\t".($NovelLoc[$nvloc]-$nvOffset)."\tNCNV\n";
			$nvloc++;
		}elsif($nvloc<@NovelLoc && $NovelLoc[$nvloc]-$nvOffset<=$chr->[1] && $modLoc<@{$mod{$chr->[0]}} && $NovelLoc[$nvloc]-$nvOffset>$mod{$chr->[0]}[$modLoc][0] && $NovelLoc[$nvloc]-$nvOffset<=$mod{$chr->[0]}[$modLoc][1]){
			$nvloc++;
		}elsif($cnvloc<@CodingLoc && $modLoc<@{$mod{$chr->[0]}} && $CodingLoc[$cnvloc]-$cnvOffset<=$mod{$chr->[0]}[$modLoc][1]-$mod{$chr->[0]}[$modLoc][0]){
			print NVSIM "$chr->[0]\t".($CodingLoc[$cnvloc]-$cnvOffset+$mod{$chr->[0]}[$modLoc][0])."\tCDSNV\n";
			$cnvloc++;
		}elsif($modLoc<@{$mod{$chr->[0]}} && ($cnvloc>=@CodingLoc || $CodingLoc[$cnvloc]-$cnvOffset>$mod{$chr->[0]}[$modLoc][1]-$mod{$chr->[0]}[$modLoc][0])){
			$cnvOffset+=$mod{$chr->[0]}[$modLoc][1]-$mod{$chr->[0]}[$modLoc][0];
			$modLoc++;
		}else{
			last;
		}
	}
	$nvOffset+=$chr->[1];
}
close(NVSIM);

print STDERR "Rewriting parameters...\n";
open(IPARAM,">$ARGV[1].pgsim.allparams.conf") || die;

print IPARAM "Genome=$Genome\n";
print IPARAM "GModel=$GModel\n";
print IPARAM "allvarDB=$allvarDB\n";
print IPARAM "commonDB=$commonDB\n";
print IPARAM "svDB=$svDB\n";
print IPARAM "pvDB=$pvDB\n";

print IPARAM "Gender=$Gender\n";
print IPARAM "Population=$Population\n";
print IPARAM "OVR=$OVR\n";
print IPARAM "KVR=$KVR\n";
print IPARAM "CVR=$CVR\n";
print IPARAM "SVN=$SVN\n";
print IPARAM "PVN=$PVN\n";

print IPARAM "CDN=$CDN\n";
print IPARAM "NIDL=$NIDL\n";
print IPARAM "NIDR=$NIDR\n";
print IPARAM "FSR=$FSR\n";
print IPARAM "TiTv=$TiTv\n";
print IPARAM "TiTvC=$TiTvC\n";

print IPARAM "GLN=$GLN\n";
print IPARAM "nGLN=$nGLN\n";
print IPARAM "modLen=$modLen\n";
print IPARAM "varNum=$varNum\n";
print IPARAM "CVT=$CVT\n";
print IPARAM "PLalpha=$PLalpha\n";
print IPARAM "cvNum=$cvNum\n";
print IPARAM "ecvNum=$ecvNum\n";
print IPARAM "eoR=$eoR\n";
print IPARAM "paC=$paC\n";
print IPARAM "paAF=$paAF\n";
print IPARAM "pnvNum=$pnvNum\n";

close(IPARAM);

print STDERR "Planning done.\n";

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

__END__