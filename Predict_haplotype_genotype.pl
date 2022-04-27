#!/usr/bin/perl -w
use strict;
use lib ("./");
use Boundary qw(boundary);

die "<fetal DNA concentration> <family merge vcf> <maternal haplotypes constructed by stLFR> <paternal haplotypes constructed by stLFR> <cfDNA name in vcf:TDP1807015083_cfDNA> <maternal name in vcf:TDP1807015084_M_stLFR> <paternal name in vcf:TDP1807015714_P_stLFR> " unless @ARGV==7;

my ($fetal,$combinevcf,$m_hap,$f_hap,$cfdna_name,$m_name,$f_name)=@ARGV;
my ($cfdna_col,$m_col,$f_col);
my (%block,%hap,%hap_f,%C3);
my ($flag,$b,$start,$end)=(0,0);
my ($p,$q,$un,$ex,$un_C2_1,$un_C2_2,$un_C2_3,$GJ)=(0,0.000001,0,0,0,0,0);

open(IN2,"zcat -f $m_hap|")||die"Cannot open input file  $m_hap:$!\n";
while(<IN2>){
        chomp;
	next if (/^\s*$/);
	my @t=split(/\s+/,$_);

	/^#/ and next;
	next if(/HOMVAR/);
        next if(/HETVAR/);
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$mat) = @t[0,1,2,3,4,5,6,7,8,9];
	next if($alt=~m/\,/);
        my %format;
        my @f = split /:/, $t[8];
        @format{@f} = 0..$#f;
	my ($BI_id,$AL1_id,$AL2_id)=($format{"BI"}, $format{"AL1"},$format{"AL2"});
        my ($bi_mat, $al1_mat,$al2_mat) = (split /:/, $t[9])[$BI_id ,$AL1_id, $AL2_id];


	my $gt;
	my $hap1;
	if($al1_mat eq $ref and $al2_mat eq $alt){
		$gt="0/1";
		$hap1=0;
	}elsif($al1_mat eq $alt and $al2_mat eq $ref) {
                $gt="1/0";
		$hap1=1;
        }elsif(length($ref)>1 and length($alt)==1){#Deletion
                if(substr($al1_mat,1) eq substr($ref,1) and $al2_mat eq $alt){
                        $gt="1/0";
			$hap1=1;
                }elsif(substr($al2_mat,1) eq substr($ref,1) and $al1_mat eq $alt){
                        $gt="0/1";
			$hap1=0;
                }
        }elsif(length($alt)>1 and length($ref)==1 ){#Insertion
		if(substr($al1_mat,1) eq substr($alt,1) and $al2_mat eq $ref){
        	        $gt="1/0";
			$hap1=1;
		}elsif(substr($al2_mat,1) eq substr($alt,1) and $al1_mat eq $ref){
			$gt="0/1";
			$hap1=0;
		}
        }
		$hap{$chr."_".$pos}={"block"=>$bi_mat,"hap1"=>$hap1};	
}
close IN2;

($flag,$b)=(0,0);
open(IN3,"zcat -f $f_hap|")||die"Cannot open input file $f_hap:$!\n";
while(<IN3>){
        chomp;
        next if (/^\s*$/);
        my @t=split(/\s+/,$_);

	        /^#/ and next;
        next if(/HOMVAR/);
        next if(/HETVAR/);
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$mat) = @t[0,1,2,3,4,5,6,7,8,9];
        next if($alt=~m/\,/);
        my %format;
        my @f = split /:/, $t[8];
        @format{@f} = 0..$#f;
        my ($BI_id,$AL1_id,$AL2_id)=($format{"BI"}, $format{"AL1"},$format{"AL2"});
        my ($bi_mat, $al1_mat,$al2_mat) = (split /:/, $t[9])[$BI_id ,$AL1_id, $AL2_id];


        my $gt;
        my $hap1;
        if($al1_mat eq $ref and $al2_mat eq $alt){
                $gt="0/1";
                $hap1=0;
        }elsif($al1_mat eq $alt and $al2_mat eq $ref) {
                $gt="1/0";
                $hap1=1;
        }elsif(length($ref)>1 and length($alt)==1){#Deletion
                if(substr($al1_mat,1) eq substr($ref,1) and $al2_mat eq $alt){
                        $gt="1/0";
                        $hap1=1;
                }elsif(substr($al2_mat,1) eq substr($ref,1) and $al1_mat eq $alt){
                        $gt="0/1";
                        $hap1=0;
                }
        }elsif(length($alt)>1 and length($ref)==1 ){#Insertion
                if(substr($al1_mat,1) eq substr($alt,1) and $al2_mat eq $ref){
                        $gt="1/0";
                        $hap1=1;
                }elsif(substr($al2_mat,1) eq substr($alt,1) and $al1_mat eq $ref){
                        $gt="0/1";
                        $hap1=0;
                }
        }
                $hap_f{$chr."_".$pos}={"block"=>$bi_mat,"hap1"=>$hap1};
}
close IN3;

=cut1
open(Inf,"zcat -f $infer|") || die"Cannot open input file $infer:$!\n";
while(<Inf>){
        chomp;
        my @t=split(/\t/,$_);
        
        #if($t[11]=~m/0\/1/ && ( $t[9]=~m/1\/1/ || $t[9]=~m/0\/0/ )){
            my $key_infer="$t[0]_$t[1]";
            my $ref=$t[2];
            my ($ca,$cb)=split(//,$t[13]);
            my ($ma,$mb)=split(/\//,$t[9]);
            my ($fa,$fb)=split(/\//,$t[11]);
            if(($ma eq $mb) && ($fa ne $fb)) {
                if(($ca eq $cb) && ($ca eq $ref)){
                    $C3{$key_infer}="00";
                }elsif(($ca eq $cb) && ($ca ne $ref)){
                    $C3{$key_infer}="11";
                }
                else{
                    $C3{$key_infer}="01";  
                }
           }
}      
close Inf;
=cut

open(IN1,"zcat -f $combinevcf |")||die"Cannot open input file $combinevcf:$!\n";
while(<IN1>){
        chomp;
        next if (/^\s*$/);
	if (/^##/){
	#	print "$_\n";
		next;
	}elsif(/^#/){
		print "$_\tChild\n";
		my @tmp=split(/\s+/,$_);
		foreach my $i(1..$#tmp){
			if($tmp[$i] eq $cfdna_name){
				$cfdna_col=$i;
			}elsif($tmp[$i] eq $m_name){
				$m_col=$i;
			}elsif($tmp[$i] eq $f_name){
				$f_col=$i;
			}
		}
		next;
	}
	my @tmp=split(/\s+/,$_);
	next if($tmp[4]=~m/\,/);
	$tmp[$cfdna_col]=~s/\|/\//;
	$tmp[$f_col]=~s/\|/\//;
	$tmp[$m_col]=~s/\|/\//;
	
if ($tmp[$cfdna_col]=~m/^.\/.:([0-9]+),([0-9]+)/ ){if($1+$2<15){print join("\t",@tmp,"lowdepth"),"\n";next;}}
if ($tmp[$f_col]=~m/^.\/.:([0-9]+),([0-9]+)/ ){if($1+$2<15){print join("\t",@tmp,"lowdepth"),"\n";next;}}
if ($tmp[$m_col]=~m/^.\/.:([0-9]+),([0-9]+)/ ){if($1+$2<15){print join("\t",@tmp,"lowdepth"),"\n";next;}}

if ($tmp[$f_col]=~m/^0\/1:([0-9]+),([0-9]+)/ ){
	if($1/($1+$2)>0.8){
		print join("\t",@tmp,"alleleblance"),"\n";next;
		$tmp[$f_col]=~s/0\/1/0\/0/;
	}elsif($1/($1+$2)<0.2){
		print join("\t",@tmp,"alleleblance"),"\n";next;
		$tmp[$f_col]=~s/0\/1/1\/1/;	
	}
}
if ($tmp[$m_col]=~m/^0\/1:([0-9]+),([0-9]+)/ ){
	if($1/($1+$2)>0.8){
		print join("\t",@tmp,"alleleblance"),"\n";next;
		$tmp[$m_col]=~s/0\/1/0\/0/;
	}elsif( $1/($1+$2)<0.2){
		print join("\t",@tmp,"alleleblance"),"\n";next;
		$tmp[$m_col]=~s/0\/1/1\/1/;
	}
}

	if($tmp[$f_col]=~m/0\/1/ && ( $tmp[$m_col]=~m/1\/1/ || $tmp[$m_col]=~m/0\/0/ )){	# Categories 3
		my $key="$tmp[0]_$tmp[1]";
		my ($b,$Hap);
		my ($f,$genotype)=("NA","NA");
	my $cutoff=2;
               if($tmp[$cfdna_col]=~m/.\/.:([0-9]+),([0-9]+)/){
                        if($1 >$cutoff && $2 >$cutoff){
                                if ($tmp[$m_col]=~m/1\/1/){$f=0;
                                }elsif($tmp[$m_col]=~m/0\/0/){$f=1;
                                }
                                $genotype="0/1";
                        }elsif($1 >$cutoff && $2 <=$cutoff) {
				$f=0; $genotype="0/0";
                        }elsif($1 <=$cutoff && $2 >$cutoff) {
				$f=1; $genotype="1/1";
                        #}else{ print join("\t",@tmp,"unknown:C3"),"\n";next;
			}

			if(exists $hap_f{$key} && $hap_f{$key}{"hap1"} eq $f) {
				$b=$hap_f{$key}{"block"};
				$Hap="Hap1";
			}elsif(exists $hap_f{$key} && $hap_f{$key}{"hap1"} eq abs($f-1)) {
				$Hap="Hap2";
				$b=$hap_f{$key}{"block"};
			}else {	
                              #  $f="NA";
                               # $genotype="NA";
				$Hap="NA";
				$b="not_in_block";
			}
		print join("\t",@tmp,"$genotype:C3:$b:$f:$Hap"),"\n";
               }
	}elsif( ($tmp[$f_col]=~m/0\/0/ || $tmp[$f_col]=~m/1\/1/) && $tmp[$m_col]=~m/0\/1/){# Categories 4 
		my ($genotype,$Hap,$type,$refcount,$altcount);
		my $m="NA";my $f="NA";
		my $key="$tmp[0]_$tmp[1]";
                unless (exists $hap{$key}){
                        print join("\t",@tmp,"unknown:C4:not_in_block"),"\n";   next;
                }

		if($tmp[$cfdna_col]=~m/^.\/.:([0-9]+),([0-9]+)/){
			($refcount,$altcount)=($1,$2);#print "$refcount\t$altcount\t";
			if($1==0 || $2==0){	
				print join("\t",@tmp,"unknown_error:C4"),"\n";next;
			}
		}
		my $b=$hap{$key}{"block"};
		if($tmp[$f_col]=~m/0\/0/ && $hap{$key}{"hap1"} eq 0){#alpha
			$type="alpha";$f=0;$m="Hap1M==0";
			$block{$b}{"alpha"}{"hap1"}+=$refcount;
			$block{$b}{"alpha"}{"hap2"}+=$altcount;
			my $total=$block{$b}{"alpha"}{"hap1"}+$block{$b}{"alpha"}{"hap2"};
			$Hap=&boundary($fetal,$type,$block{$b}{"alpha"}{"hap1"},$total);
			if ($Hap eq "Hap1"){$genotype="0/0";$block{$b}{"alpha"}{"hap1"}=0;$block{$b}{"alpha"}{"hap2"}=0;
			}elsif($Hap eq "Hap2"){$genotype="0/1";$block{$b}{"alpha"}{"hap1"}=0;$block{$b}{"alpha"}{"hap2"}=0;
			}elsif($Hap eq "Unclassified"){$genotype="0";
			}else{$genotype="unknown";}
		}elsif($tmp[$f_col]=~m/1\/1/ && $hap{$key}{"hap1"} eq 1){#alpha
			$type="alpha";$f=1;$m="Hap1M==1";
			$block{$b}{"alpha"}{"hap1"}+=$altcount;
			$block{$b}{"alpha"}{"hap2"}+=$refcount;
			my $total=$block{$b}{"alpha"}{"hap1"}+$block{$b}{"alpha"}{"hap2"};
			$Hap=&boundary($fetal,$type,$block{$b}{"alpha"}{"hap1"},$total);
			if ($Hap eq "Hap1"){$genotype="1/1";$block{$b}{"alpha"}{"hap1"}=0;$block{$b}{"alpha"}{"hap2"}=0;
                        }elsif($Hap eq "Hap2"){$genotype="0/1";$block{$b}{"alpha"}{"hap1"}=0;$block{$b}{"alpha"}{"hap2"}=0;
                        }elsif($Hap eq "Unclassified"){$genotype="1";
                        }else{$genotype="unknown";}
		}elsif($tmp[$f_col]=~m/0\/0/ && $hap{$key}{"hap1"} eq 1){#beta
			$type="beta";$f=0;$m="Hap1M==1";

			$block{$b}{"beta"}{"hap1"}+=$altcount;
			$block{$b}{"beta"}{"hap2"}+=$refcount;
			my $total=$block{$b}{"beta"}{"hap1"}+$block{$b}{"beta"}{"hap2"};
			$Hap=&boundary($fetal,$type,$block{$b}{"beta"}{"hap1"},$total);
			if ($Hap eq "Hap1"){$genotype="0/1";$block{$b}{"beta"}{"hap1"}=0;$block{$b}{"beta"}{"hap2"}=0;
                        }elsif($Hap eq "Hap2"){$genotype="0/0";$block{$b}{"beta"}{"hap1"}=0;$block{$b}{"beta"}{"hap2"}=0;
                        }elsif($Hap eq "Unclassified"){$genotype="0";
                        }else{$genotype="unknown";}
		}elsif( $tmp[$f_col]=~m/1\/1/ && $hap{$key}{"hap1"} eq 0){#beta
			$type="beta";$f=1;$m="Hap1M==0";
			$block{$b}{"beta"}{"hap1"}+=$refcount;
			$block{$b}{"beta"}{"hap2"}+=$altcount;	
			my $total=$block{$b}{"beta"}{"hap1"}+$block{$b}{"beta"}{"hap2"};
			$Hap=&boundary($fetal,$type,$block{$b}{"beta"}{"hap1"},$total);
			if ($Hap eq "Hap1"){$genotype="0/1";$block{$b}{"beta"}{"hap1"}=0;$block{$b}{"beta"}{"hap2"}=0;
                        }elsif($Hap eq "Hap2"){$genotype="1/1";$block{$b}{"beta"}{"hap1"}=0;$block{$b}{"beta"}{"hap2"}=0;
                        }elsif($Hap eq "Unclassified"){$genotype="1";
                        }else{$genotype="unknown";}
		}else{ 	$type="unknown";	
			$genotype="unknown";
			$Hap="Unclassified";
		}	
		print join("\t",@tmp,"$genotype:C4:$f:$b:$m:$Hap"),"\n";
	}elsif($tmp[$f_col]=~m/0\/1/ && $tmp[$m_col]=~m/0\/1/ ){# Categories 5
		my $key="$tmp[0]_$tmp[1]";
		if (exists $hap_f{$key} && exists $hap{$key}){
			my $b_f=$hap_f{$key}{"block"};
			my $b=$hap{$key}{"block"};
			my $f_y=$hap_f{$key}{"hap1"};
			my $m_y=$hap{$key}{"hap1"};
			print join("\t",@tmp,"unknown:C5:$b_f:$b:Hap1F==$f_y:Hap1M==$m_y"),"\n";
		}else{
			print join("\t",@tmp,"unknown:C5:not_in_block"),"\n";
		}
	}else{
		print join("\t",@tmp,"unknown"),"\n";
	}
}
close IN1;

#print "#########################\n$GJ\n";
