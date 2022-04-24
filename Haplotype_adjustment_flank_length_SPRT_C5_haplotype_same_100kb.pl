#!/usr/bin/perl -w
use strict;
use lib ("./");
use Boundary qw(boundary);


die "<input> <length> <fetal DNA concentration> < cfDNA name in vcf>"  unless @ARGV==4;
my ($combinevcf,$length,$ff,$cfdna_name)=@ARGV;
my (%hash,%C3,%C4,$len,$cfdna_col,%Block);
open(IN1,"$combinevcf ")||die"Cannot open input file:$!\n";
while(<IN1>){
        chomp;
        next if (/^\s*$/);
        if (/^##/){
                 next;
        }elsif(/^#/){
                print $_,"\n";
		my @tmp=split(/\s+/,$_);
                foreach my $i(1..$#tmp){
                        if($tmp[$i] eq $cfdna_name){
                                $cfdna_col=$i;
                        }
                }
		next;
        }
	my @tmp=split(/\s+/,$_);
#	next if($tmp[-1] =~m/unknown/);
	$hash{$tmp[0]}{$tmp[1]}=\@tmp;
	#print join("\t",@{$hash{$tmp[0]}{$tmp[1]}}),"\n";
	if(/C3/){
		if(/:C3:(chr[0-9_X]+):[01]:(Hap[12])/){
			my $b=$1; 
			$C3{$b}{$tmp[1]}=$2;
		}
	}elsif(/C4/){
		if(/:C4:[01]:(chr[0-9_X]+):Hap1M==[01]:(Hap[12])/){	
			my $b=$1; 
			$C4{$b}{$tmp[1]}=$2;
		#	if(! exists $C4{$b}){$C4{$b}{"Hap1"}=0;$C4{$b}{"Hap2"}=0;}
		}
	}	
}
close IN1;

sub  hap{
	my ($block,$pos,$type)=@_;
	my %h;
	my $class="unknown";
	if($type eq "C4"){
		%h=%C4;
	}elsif($type eq "C3"){
		%h=%C3;
	}
	my $min=$length;
	my $len=$length;
	my ($p1,$p2)=(0,0);
	my ($class1,$class2)=("unknown","unknown");
         for (my $i=0;$i<=$len;$i++){
                $p1=$pos+$i;
                if (exists $h{$block}{$p1} ){
                        $class1=$h{$block}{$p1};
                last;
                }
	}
	for (my $i=0;$i<=$len;$i++){
		$p2=$pos-$i;
                if (exists $h{$block}{$p2} ){
                         $class2=$h{$block}{$p2};
                last;
                }

        }

	if($class1 eq $class2 && $class1 ne "unknown"){
		$class=$class1;
	}elsif($class1 ne $class2 && $class1 eq "unknown"){
               $class=$class2;
	}elsif($class1 ne $class2 && $class2 eq "unknown"){
              $class=$class1;
	}

=ped
	for my $p(sort{$a cmp $b} keys %{$h{$block}}){
		my $l=$p-$pos;
		if(abs($l)<1000 && abs($l)<$min){
			$min=abs($l);
			$class=$h{$block}{$p};
		}
	}
=cut
	return ($class) ;

}

for my $chr (sort{$a cmp $b} keys %hash){
	for my $pos (sort {$a<=>$b} keys %{$hash{$chr}}){
		my ($refcount,$altcount);
		if(${$hash{$chr}{$pos}}[$cfdna_col]=~m/^.\/.:([0-9]+),([0-9]+)/){
			($refcount,$altcount)=($1,$2);

		}
		if (${$hash{$chr}{$pos}}[12] =~m/C4/){			
			if (${$hash{$chr}{$pos}}[12] =~m/([01]):C4:[01]:(chr[0-9_X]+):Hap1M==([01]):Unclassified/ && defined $C4{$2}){
		
				my $b=$2;
				my $fa_geno=$1;
				my $hap1m=$3;
=ped
				if ($fa_geno eq 0 && $hap1m eq 0){#alpha
					$Block{$b}{"alpha"}{"hap1"}+=$refcount;
					$Block{$b}{"alpha"}{"hap2"}+=$altcount;
				}elsif($fa_geno eq 1 && $hap1m eq 1){#alpha
					$Block{$b}{"alpha"}{"hap1"}+=$altcount;
					$Block{$b}{"alpha"}{"hap2"}+=$refcount;
				}elsif($fa_geno eq 0 && $hap1m eq 1){#beta
					$Block{$b}{"beta"}{"hap1"}+=$altcount;
					$Block{$b}{"beta"}{"hap2"}+=$refcount;
				}elsif($fa_geno eq 1 && $hap1m eq 0){#beta
					$Block{$b}{"beta"}{"hap1"}+=$refcount;
					$Block{$b}{"beta"}{"hap2"}+=$altcount;
				}
=cut
				my $ma_geno="unknown";
				my $H=&hap($b,$pos,"C4");

				if ($H eq "Hap1" ){
					$ma_geno=$hap1m;
				}elsif ($H eq "Hap2" ){
					$ma_geno=1-$hap1m;
				}
				my $genotype=join ("\/",sort ($fa_geno,$ma_geno));
				my $lab=join(":",$genotype,"C4.new",$b,$H."=".$ma_geno);
				print join("\t",@{$hash{$chr}{$pos}}[0..11],$lab),"\n";
			}else {
				print join("\t",@{$hash{$chr}{$pos}}),"\n";
			}
			
		}elsif (${$hash{$chr}{$pos}}[12] =~m/C5/){
			if (${$hash{$chr}{$pos}}[12] =~m/unknown:C5:(chr[0-9_X]+):(chr[0-9_X]+):Hap1F==([01]):Hap1M==([01])/ && defined $C3{$1} && defined $C4{$2}){#	print "$1\t$2\t";
				my $b1=$1;
				my $b2=$2;
				my $hap1m=$4;
				my ($fa_geno,$ma_geno)=("unknown","unknown");
				my $H1=&hap($b1,$pos,"C3");
				my $H2="unknown";
				#my $H2=&hap($b2,$pos,"C4");
				if($H1 eq "Hap1" ){
					$fa_geno=$3;
				}elsif ($H1 eq "Hap2"){
					$fa_geno=1-$3;
				}
				
				my ($type,$c5);
				if ($fa_geno eq 0 && $hap1m eq 0){#alpha
					$type="alpha";
                                        $Block{$b2}{"alpha"}{"hap1"}+=$refcount;
                                        $Block{$b2}{"alpha"}{"hap2"}+=$altcount;
					    my $total=$Block{$b2}{$type}{"hap1"}+$Block{$b2}{$type}{"hap2"};
                                        $H2=&boundary($ff,$type,$Block{$b2}{$type}{"hap1"},$total);

                                }elsif($fa_geno eq 1 && $hap1m eq 1){#alpha
					$type="alpha";
                                        $Block{$b2}{"alpha"}{"hap1"}+=$altcount;
                                        $Block{$b2}{"alpha"}{"hap2"}+=$refcount;
					    my $total=$Block{$b2}{$type}{"hap1"}+$Block{$b2}{$type}{"hap2"};
                                        $H2=&boundary($ff,$type,$Block{$b2}{$type}{"hap1"},$total);

                                }elsif($fa_geno eq 0 && $hap1m eq 1){#beta
					$type="beta";
                                        $Block{$b2}{"beta"}{"hap1"}+=$altcount;
                                        $Block{$b2}{"beta"}{"hap2"}+=$refcount;   
					my $total=$Block{$b2}{$type}{"hap1"}+$Block{$b2}{$type}{"hap2"};
                                        $H2=&boundary($ff,$type,$Block{$b2}{$type}{"hap1"},$total);

                                }elsif($fa_geno eq 1 && $hap1m eq 0){#beta
					$type="beta";
                                        $Block{$b2}{"beta"}{"hap1"}+=$refcount;
                                        $Block{$b2}{"beta"}{"hap2"}+=$altcount;
					    my $total=$Block{$b2}{$type}{"hap1"}+$Block{$b2}{$type}{"hap2"};
                                        $H2=&boundary($ff,$type,$Block{$b2}{$type}{"hap1"},$total);

                                }
				if($H2 eq "Hap1" ){
					$ma_geno=$4;
					$Block{$b2}{$type}{"hap1"}=0;
					$Block{$b2}{$type}{"hap2"}=0;
					$c5="C5";
				}elsif ($H2 eq "Hap2"){
					$ma_geno=1-$4;
					$Block{$b2}{$type}{"hap1"}=0;
                                        $Block{$b2}{$type}{"hap2"}=0;
					$c5="C5";
				}else{
					$H2=&hap($b2,$pos,"C4");
                                        if($H2 eq "Hap1" ){
                                                $ma_geno=$4;
					}elsif ($H2 eq "Hap2"){
                                                $ma_geno=1-$4;
                                        }
                                        $c5="C5.new";
				}
				my $genotype=join ("\/",sort ($fa_geno,$ma_geno));
				my $lab=join(":",$genotype,$c5,$1,$2,$H1."F=".$fa_geno,$H2."M=".$ma_geno);
				print join("\t",@{$hash{$chr}{$pos}}[0..11],$lab),"\n";
			}else {
				print join("\t",@{$hash{$chr}{$pos}}),"\n";
			}
		}else{
			print join("\t",@{$hash{$chr}{$pos}}[0..12]),"\n";
		}
	} 
} 
