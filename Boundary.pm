package Boundary;
use strict;use warnings;
require Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK=qw(boundary);


#my $fetal=0.1143;
#my $N=2369;
#my $snptype="alpha";
#my $hap1count=1285;
#my @result=boundary($fetal,$snptype,$hap1count,$N);
#print "$hap1count\t$N\t$snptype\t@result\n";


#$hap1count=$N-$hap1count;
#$snptype="beta";
#@result=boundary($fetal,$snptype,$hap1count,$N);
#print "$hap1count\t$N\t$snptype\t@result\n";

sub  boundary{
	my ($fetal,$snptype,$hap1count,$N)=@_;
	my ($q1,$q0)=();
	if ($snptype eq "alpha"){
		$q0=0.5;
		$q1=0.5+$fetal/2;
	}elsif($snptype eq "beta"){
		$q1=0.5;
		$q0=0.5-$fetal/2;
	}else{
		die "wrong";
	}
	my $d=(1-$q1)/(1-$q0);
	my $g=$q1*(1-$q0)/$q0/(1-$q1);
	my $upper=(log(1200)/$N-log($d))/log($g);
	my $lower=(log(1/1200)/$N-log($d))/log($g);
	my $frac=$hap1count/$N;
	my $class;
	if ($frac < $lower){
		$class="Hap2";
	}elsif($frac > $upper){
		$class="Hap1";
	}else{
		$class="Unclassified";		
	}
	#return ($class,$lower,$upper);
	return ($class);
}
