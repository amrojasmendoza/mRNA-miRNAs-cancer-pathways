#
#  Saca_grupos
#
#  Created by Eduardo Andrés León on 2016-10-17.
#  Copyright (c) 2016 IPBLN. All rights reserved.
#
#!/usr/bin/perl
$|=1;
use strict;
use Getopt::Long;
use Ubio::Utils::date;
use List::MoreUtils qw(uniq);

my $help;
my $file;
my $type="one";
my $paper=undef;

my $all_genes;
my $all_miRNAs;
my $all_pares;

GetOptions(
	"help" => \$help,
	"file|f=s" => \$file,
	"type|t=s" => \$type,
	"paper"	=>\$paper
);
if($file){
	
	my $names=rename_tumors();
	
	my $data;
	my $all_data;
	my $my_data;
	

	if($type eq "unique"){
		my $pares;
		open(FILE,$file) || die $!;
		while(<FILE>){
			chomp;
			$_=~s/\"//g;
			if($_ !~/^Association/){
				my ($par,undef,undef,undef,undef,undef,$geneUpmiRDown,$geneDownmiRUp)=split(/\t/);
				$pares->{$par}++;
				my($gen,$miRNA)=split(/_/,$par);
				if($geneUpmiRDown){
					if($geneUpmiRDown =~/,/){
						my @tumores_all=split(/,/,$geneUpmiRDown);
						my @tumores=uniq(@tumores_all);
						for (my $i = 0; $i <=$#tumores; $i++) {
							for (my $j = 1; $j <=$#tumores; $j++) {
								if($i !=$j and $i<$j){
									#print "$par $i $tumores[$i]\t$j $tumores[$j]\n";
									$my_data->{$par}->{"$tumores[$i]\t$tumores[$j]"}++;
									#$all_data->{$tumores[$i]}->{$par}++;
									#$all_data->{$tumores[$j]}->{$par}++;
									$all_genes->{$gen}->{$tumores[$i]}++;
									$all_genes->{$gen}->{$tumores[$j]}++;
									
									$all_miRNAs->{$miRNA}->{$tumores[$i]}++;
									$all_miRNAs->{$miRNA}->{$tumores[$j]}++;
									
									$all_pares->{$miRNA}->{$par}++;
									$all_pares->{$gen}->{$par}++;
									
								}
							}
						}
					}else{
						#Aqui esto no es valido, porque buscamos que tengan minimo 2 tumores
						#$all_data->{$geneUpmiRDown}->{$par}++;
						$all_miRNAs->{$miRNA}->{$geneUpmiRDown}++;
						$all_genes->{$gen}->{$geneUpmiRDown}++;
						
					}
				}
				if($geneDownmiRUp){
					if($geneDownmiRUp =~/,/){
						my @tumores_all=split(/,/,$geneDownmiRUp);
						my @tumores=uniq(@tumores_all);
						for (my $i = 0; $i <=$#tumores; $i++) {
							for (my $j = 1; $j <=$#tumores; $j++) {
								if($i !=$j and $i<$j){
									#print "$par $i $tumores[$i]\t$j $tumores[$j]\n";
									$my_data->{$par}->{"$tumores[$i]\t$tumores[$j]"}++;
									#$all_data->{$tumores[$i]}->{$par}++;
									#$all_data->{$tumores[$j]}->{$par}++;
									
									$all_genes->{$gen}->{$tumores[$i]}++;
									$all_genes->{$gen}->{$tumores[$j]}++;
									
									$all_miRNAs->{$miRNA}->{$tumores[$i]}++;
									$all_miRNAs->{$miRNA}->{$tumores[$j]}++;
									
									$all_pares->{$miRNA}->{$par}++;
									$all_pares->{$gen}->{$par}++;
									
								}
							}
						}
					}else{
						#Aqui esto no es valido, porque buscamos que tengan minimo 2 tumores
						#$all_data->{$geneDownmiRUp}->{$par}++;
						$all_miRNAs->{$miRNA}->{$geneDownmiRUp}++;
						$all_genes->{$gen}->{$geneDownmiRUp}++;
						
					}
				}
			}
		}
		close FILE;
		foreach my $par(keys %{$my_data}){
			if(scalar(keys %{$my_data->{$par}})==1){
				foreach my $tumor (keys %{$my_data->{$par}}){

					$data->{$tumor}->{$par}++;

					my($tum1,$tum2)=split(/\t/,$tumor);
					$all_data->{$tum1}->{$par}++;
					$all_data->{$tum2}->{$par}++;
				}
				#print $par ."\t" . scalar(keys %{$data->{$par}}) ."\n" );
			}
		}
		my $luad_lusc;
		
		foreach my $tum (sort {$data->{$b} <=> $data->{$a}} keys %{$data}){
			my($tum1,$tum2)=split(/\t/,$tum);
			my $valorA=scalar(keys %{$data->{$tum}});
			my $valorB=scalar(keys %{$all_data->{$tum1}}) +  scalar(keys %{$all_data->{$tum2}});
			if(!$paper){
				print $tum . "\t" . $names->{$tum1} ."\t". $names->{$tum2} ."\t". ($valorA/$valorB)."\t".$valorA."\t[". scalar(keys %{$all_data->{$tum1}}) ."]\t[". scalar(keys %{$all_data->{$tum2}}) ."]\t" . join(",",keys %{$data->{$tum}}) ."\n";
			}
			if($paper){
				foreach my $par (keys %{$data->{$tum}}){
					if($tum eq "Lung Adenocarcinoma\tLung squamous cell carcinoma" or $tum eq "Lung squamous cell carcinoma\tLung Adenocarcinoma"){
						#print "$tum\t$par\n";
						$luad_lusc->{$tum}->{$par}++;
					}

				}
			}
		}
		if($paper){
			important($luad_lusc,"Lung squamous cell carcinoma");
			
		}
	}
	else{
		open(FILE,$file) || die $!;
		while(<FILE>){
			chomp;
			$_=~s/\"//g;
			if($_ !~/^Association/){
				my ($par,undef,undef,undef,undef,undef,$geneUpmiRDown,$geneDownmiRUp)=split(/\t/);
				if($geneUpmiRDown){
					if($geneUpmiRDown =~/,/){
						my @tumores_all=split(/,/,$geneUpmiRDown);
						my @tumores=uniq(@tumores_all);
						for (my $i = 0; $i <=$#tumores; $i++) {
							for (my $j = 1; $j <=$#tumores; $j++) {
								if($i !=$j and $i<$j){
									#print "$par $i $tumores[$i]\t$j $tumores[$j]\n";
									$data->{"$tumores[$i]\t$tumores[$j]"}->{$par}++;
									$all_data->{$tumores[$i]}->{$par}++;
									$all_data->{$tumores[$j]}->{$par}++;
								}
							}
						}
					}
				}
				if($geneDownmiRUp){
					if($geneDownmiRUp =~/,/){
						my @tumores_all=split(/,/,$geneDownmiRUp);
						my @tumores=uniq(@tumores_all);
						for (my $i = 0; $i <=$#tumores; $i++) {
							for (my $j = 1; $j <=$#tumores; $j++) {
								if($i !=$j and $i<$j){
									#print "$par $i $tumores[$i]\t$j $tumores[$j]\n";
									$data->{"$tumores[$i]\t$tumores[$j]"}->{$par}++;
									$all_data->{$tumores[$i]}->{$par}++;
									$all_data->{$tumores[$j]}->{$par}++;
								}
							}
						}
					}
				}
			}
		}
		close FILE;
		print "Tumor1\tTumor2\tAbrev1\tAbrev2\tParesNorm\tTotalPares\tTotal tumor1\tTotal tumor2\n";
		#print STDERR "Tengo un total de " . scalar(keys %{$data}) . " tumores ($cont)\n";
		foreach my $tum (sort {$data->{$b} <=> $data->{$a}} keys %{$data}){
			my($tum1,$tum2)=split(/\t/,$tum);
			#Defino valor normalizado como el ValorA/ValorB. Siendo ValorA el numero total de pares compartidos entre los dos tumores y
			#ValorB seria la suma del total de pares posibles por cada uno de los 2 tumores implicados.
			my $valorA=scalar(keys %{$data->{$tum}});
			my $valorB=scalar(keys %{$all_data->{$tum1}}) +  scalar(keys %{$all_data->{$tum2}});
			print $tum . "\t" . $names->{$tum1} ."\t". $names->{$tum2} ."\t". ($valorA/$valorB)."\t".$valorA."\t[". scalar(keys %{$all_data->{$tum1}}) ."]\t[". scalar(keys %{$all_data->{$tum2}}) ."]\n";# . join(",",keys %{$data->{$tum}}) ."\n";
		}
	}
}
else{
	print STDERR "Error: Format:\nPares unicos:
	\tperl total_pares_comun_entre_tumores.pl -f TCGA_pairs_DE_pathway_agreement2.xls -t unique [-paper]
Pares compartidos:
	\tperl $0 -f TCGA_pairs_DE_pathway_agreement2.xls\n\n";
}


sub rename_tumors{
	my $hash;
	
	$hash->{"Bladder Urothelial Carcinoma"}="BLCA";
	$hash->{"Breast invasive carcinoma"}="BRCA";
	$hash->{"Cholangiocarcinoma"}="CHOL";
	$hash->{"Esophageal Carcinoma"}="ESCA";
	$hash->{"Head and Neck squamous cell carcinoma"}="HNSC";
	$hash->{"Kidney Chromophobe"}="KICH";
	$hash->{"Kidney renal clear cell carcinoma"}="KIRC";
	$hash->{"Kidney renal papillary cell carcinoma"}="KIRP";
	$hash->{"Liver hepatocellular carcinoma"}="LIHC";
	$hash->{"Lung Adenocarcinoma"}="LUAD";
	$hash->{"Lung squamous cell carcinoma"}="LUSC";
	$hash->{"Prostate adenocarcinoma"}="PRAD";
	$hash->{"Stomach adenocarcinoma"}="STAD";
	$hash->{"Thyroid carcinoma"}="THCA";
	$hash->{"Uterine Corpus Endometrial carcinoma"}="UCEC";
	
	return($hash);
	
}

sub important{
	
	my $chol_lihc=shift;
	my $label=shift;
	#my $pathways=get_pathways();
	my $names=rename_tumors();
	my $expression=get_expression($label);
	
	my $chol_lihc_genes;
	my $chol_lihc_miRNAs;
	my $prueba;
	foreach my $tumor_pair(keys %{$chol_lihc}){
		foreach my $pair(keys %{$chol_lihc->{$tumor_pair}}){
			my($gene,$miRNA)=split(/_/,$pair);
			$chol_lihc_genes->{$gene}->{$pair}++;
			$chol_lihc_miRNAs->{$miRNA}->{$pair}++;
			
			my $miRNA_gene=$miRNA;
			$miRNA_gene=~s/\-3p//g;
			$miRNA_gene=~s/\-5p//g;
			#$chol_lihc_miRNAs->{$miRNA_gene}++;
		}
	}
	
	open(OUT2,">lung.pairs.csv") || die $!;
	foreach my $tumor_pair(keys %{$chol_lihc}){
		foreach my $pair(sort keys %{$chol_lihc->{$tumor_pair}}){
			my($gene,$miRNA)=split(/_/,$pair);
			print OUT2 $gene . "," . $miRNA ."\n";
		}
	}
	close OUT2;
	return();
}

sub get_expression{
	my $tumor=shift;
	my $tumor1=$tumor;
	$tumor1=~s/ /_/g;
	open(FILE,"Differential-expression/$tumor1/$tumor RNASeq\.xls")|| die "Cant find Differential-expression/$tumor1/$tumor RNASeq\.xls: $!\n";
	my $exp="NA";
	my $data;
	while(<FILE>){
		chomp;
		my ($gene,$path,$logFC)=split(/\t/);
		if($logFC>0){
			$exp="+";
		}
		else{
			$exp="-";
		}
		$data->{$gene}=$exp;
	}
	close FILE;
	return($data);
}
