#
#  3.Summary_miRGate
#
#  Created by Eduardo Andrés León on 2015-02-21.
#  Copyright (c) 2015 CNIO. All rights reserved.
#
#!/usr/bin/perl
$|=1;
use strict;
use Getopt::Long;
use Ubio::Utils::date;
use List::MoreUtils qw/ uniq /;

my $help;
my $dirname="";
my $FDR_cutoff=0.05;
my $FC_cutoff=1;
my $FC_cutoff_neg=-1;
my $agreement_cutoff=2;
my $type;
my $data="RNASeq";
my $net;
GetOptions(
	"help" => \$help,
	"data|d=s" => \$data,
	"type|t=s" =>\$type,
	"net|n=s" => \$net,
);
if($type and $data){	
	
	if(lc($type) eq "de"){
		$dirname="Differential-expression/";
	}
	if(lc($type) eq "pairs"){
		$dirname="miRGate_miRNA-mRNA-pairs/";
	}
	my $common_miRNAs;
	my $common_genes;
	my $GenesDE;
	my $miRNAsDE;
	
	my $common_genesUP;
	my $common_genesDown;
	my $common_miRNAsUP;
	my $common_miRNAsDown;

	my $GenesUp;
	my $miRNAsUp;
	my $GenesDown;
	my $miRNAsDown;
	
	my $genes_in_tumor;
	my $miRNAs_in_tumor;
	my $pairs_in_tumor;
	
	my $pairs_genes;
	my $pairs_miRNAs;
	
	my $type_used;
	my $gen_per_tumor;
	
	my $DDR_network;
	my $DDR_paths;
	
	my $common_pairs_trDw_miRUp;
	my $common_pairs_trUp_miRDw;
	my $common_pairs;
	
	my $pairs_trDw_miRUp;
	my $pairs_trUp_miRDw;
	
	my @all_tumors;
	
	my $all_posible_genes;
	open(ALLGENES,"../Pathways/data/genes_pathways.txt") || die "$! Cant' find ../Pathways/data/genes_pathways.txt\n";
	while(<ALLGENES>){
		chomp;
		my($genes,$source,$path)=split(/\t/);
		$all_posible_genes->{$genes}->{$path}++;
	}
	close ALLGENES;
	
	if(lc($type) eq "pairs"){
		my $count_file=1;
		opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
		my @dires= readdir DIR;
		@dires = sort(@dires);
		#print STDERR "Reading " . (scalar(@files)-2) ." files\n";
		foreach my $tumor_folder (@dires){
			if($tumor_folder !~ /^\./ and $tumor_folder ne "README.md"){
				my $saved;
				opendir ( DIR2, $dirname ."/".$tumor_folder ) || die "Error in opening dir $dirname/$tumor_folder\n";
				my @files= readdir DIR2;
				@files = sort(@files);	
				foreach my $filename (@files){
					#Filtering for xls (miRgate output format)
				     if($filename =~ /\.xls$/){
						my $tumor=$filename;
						#Parding filename to gather the tumor name
						$tumor=~s/(.+)_(.+)_miRGate_predictions.xls$/$1/g;
						my $selected_type=$2;
						my $read_file;
						if(!$selected_type and !$net){
							$tumor=~s/(.+)_miRGate_predictions.xls$/$1/g;
							$read_file=("$dirname/$tumor_folder/$tumor\_miRGate_predictions.xls");
							$type_used="AllGenes";
						}
						elsif(lc($selected_type) eq lc($net)){
							$read_file=("$dirname/$tumor_folder/$filename");
							$type_used=$net;
						}
						else{
							next;
						}
						push(@all_tumors,$tumor);
						open(RESULTS, "$read_file") || die "$read_file : $!";
						print STDERR "$count_file) ".date." Reading $tumor data\n";
			#			print STDERR "$count_file) Reading $tumor data ($read_file)\n";
						$count_file++;
						while(<RESULTS>){
							chomp;
							if($_ !~/^Transcript/){
								my @data=split("\t");
			
								#Variable asignment
								my $transcript=$data[0];
								my $miRNA=$data[1];
								my $agreement=$data[6];
								my $methods=$data[7];
								my $ensembl=$data[9];
								my $symbol=$data[10];
								my $tr_FC=$data[14];
								my $tr_FDR=$data[17];
								my $miR_FC=$data[18];
								my $miR_FDR=$data[21];
								my $connections=$data[4];
								my $amplified=$data[13];
								$gen_per_tumor->{$symbol}=$type_used;
								#taking into account targets predicted for att least 2 different methods
								if($agreement>=$agreement_cutoff and $amplified eq "No"){
									#Filtering due FDR and FC for Transripts
									if($tr_FDR<=$FDR_cutoff){
										$common_genes->{$symbol}++;
										$GenesDE->{$tumor}->{$symbol}++;
					
										if(!exists $genes_in_tumor->{$symbol}->{$tumor}){
											$genes_in_tumor->{$symbol}->{$tumor}++;
										}
					
										if($tr_FC>0){
											$common_genesUP->{$symbol}++;
											$GenesUp->{$tumor}->{$symbol}++;
										}
										else{
											$common_genesDown->{$symbol}++;
											$GenesDown->{$tumor}->{$symbol}++;
										}
									}
									#Filtering due FDR and FC for miRNAs
			
									if($miR_FDR<=$FDR_cutoff){
					
										if(!exists $miRNAs_in_tumor->{$miRNA}->{$tumor}){
											$miRNAs_in_tumor->{$miRNA}->{$tumor}++;
										}
					
										$common_miRNAs->{$miRNA}++;
										$miRNAsDE->{$tumor}->{$miRNA}++;
					
										if($miR_FC>0){
											$common_miRNAsUP->{$miRNA}++;
											$miRNAsUp->{$tumor}->{$miRNA}++;
										}
										else{
											$common_miRNAsDown->{$miRNA}++;
											$miRNAsDown->{$tumor}->{$miRNA}++;
										}
									}
									# valid pairs
									if($miR_FDR<=$FDR_cutoff and $tr_FDR<=$FDR_cutoff){
								
										if(!exists $saved->{"$symbol\_$miRNA"}){
											$pairs_genes->{$symbol}+=$agreement;
											$pairs_miRNAs->{$miRNA}+=$agreement;
									
											$saved->{"$symbol\_$miRNA"}++;
										}
							
										if($common_pairs->{"$symbol\_$miRNA"}<$agreement){
											$common_pairs->{"$symbol\_$miRNA"}=$agreement;
										}
										if(!exists $pairs_in_tumor->{"$symbol\_$miRNA"}->{$tumor}){
											$pairs_in_tumor->{"$symbol\_$miRNA"}->{$tumor}=$agreement;
										}else{
											#save the best agreement among all transcripts
											if($pairs_in_tumor->{"$symbol\_$miRNA"}->{$tumor} < $agreement){
												$pairs_in_tumor->{"$symbol\_$miRNA"}->{$tumor}=$agreement;
											}
										}
							
										if($miR_FC>0 and $tr_FC<0){
											my $pair="$symbol\_$miRNA";
											if($common_pairs_trDw_miRUp->{$pair}<$agreement){
												$common_pairs_trDw_miRUp->{$pair}=$agreement;
											}

											if($pairs_trDw_miRUp->{$tumor}->{$pair}<$agreement){
												$pairs_trDw_miRUp->{$tumor}->{$pair}=$agreement;
											}

								
										}
										if($miR_FC<0 and $tr_FC>0){
											my $pair="$symbol\_$miRNA";
											if($pairs_trUp_miRDw->{$tumor}->{$pair}<$agreement){
												$common_pairs_trUp_miRDw->{$pair}=$agreement;
												$pairs_trUp_miRDw->{$tumor}->{$pair}=$agreement;
											}
										}
									}
								}
							}
						}
						close RESULTS;
					}
				}
			}
			closedir(DIR2);
		}
		closedir(DI2);
			
		############# most regulator miRNAs DE in most number of tumors ############
		my $cont=0;
		system("mkdir -p src/DE_data_TCGA");
		open(OUT_MIR_REG,">src/DE_data_TCGA/TCGA_miRNAs_regulators_$type_used\_pathway.xls") || die $!;
		print OUT_MIR_REG "miRNAs\tN of targets\n";
		foreach my $miRNAs (sort {$pairs_miRNAs->{$b} <=> $pairs_miRNAs->{$a}} keys %$pairs_miRNAs){
			print OUT_MIR_REG "$miRNAs\t". $pairs_miRNAs->{$miRNAs} ."\n";
		}
		close OUT_MIR_REG;
		############# most regulated DE gene in most number of tumors ############
		open(OUT_GEN_REG,">src/DE_data_TCGA/TCGA_genes_regulators_$type_used\_pathway.xls") || die $!;
		print OUT_GEN_REG "genes\tN of targets\n";
	
		foreach my $genes (sort {$pairs_genes->{$b} <=> $pairs_genes->{$a}} keys %$pairs_genes){
			print OUT_GEN_REG "$genes\t". $pairs_genes->{$genes} ."\n";
		}
		close OUT_GEN_REG;
		
		if(lc($net) eq "all"){
			print "Junta todo\n";
		}
		############# miRNAs DE in most number of tumors ############
	
		my $cont=0;
		#print '-' x 100;
		#print "\n\nTen miRNAs DE in most number of tumors:\n";
		my $miRNAs_in_tumor_ordered;
		my $miRNAs_Up_in_tumor_ordered;
		my $miRNAs_Down_in_tumor_ordered;
		foreach my $miRNAs (keys %$miRNAs_in_tumor){
			my @tumores=();
			my @tumoresUp=();
			my @tumoresDown=();
			#my @present_tumor=(uniq (keys %{$prueba->{$GeneName}}));
			foreach my $tumor (keys %{$miRNAs_in_tumor->{$miRNAs}}){
				push(@tumores,$tumor);
				push(@tumoresUp,$tumor) if(exists($miRNAsUp->{$tumor}->{$miRNAs}));
				push(@tumoresDown,$tumor)  if(exists($miRNAsDown->{$tumor}->{$miRNAs}));
			
			}
			@tumores=uniq(@tumores);
			$miRNAs_in_tumor_ordered->{$miRNAs}=\@tumores;
			$miRNAs_Up_in_tumor_ordered->{$miRNAs}=\@tumoresUp;
			$miRNAs_Down_in_tumor_ordered->{$miRNAs}=\@tumoresDown;
		}

		open(OUT_MIR,">src/DE_data_TCGA/TCGA_miRNAs_DE_$type_used\_pathway.xls") || die $!;
		print OUT_MIR "miRNAs\tnDE\tnDE Up\tnDE Down\tTumors no DE\tTumors UP\tTumors Down\n";

		foreach my $ff (sort {@{$miRNAs_in_tumor_ordered->{$b}} <=> @{$miRNAs_in_tumor_ordered->{$a}}} keys  %$miRNAs_in_tumor_ordered){
			print OUT_MIR $ff ."\t" . @{$miRNAs_in_tumor_ordered->{$ff}} ."\t". @{$miRNAs_Up_in_tumor_ordered->{$ff}} . "\t". @{$miRNAs_Down_in_tumor_ordered->{$ff}} ."\t". join(",",sort {$a cmp $b} intersec(\@{$miRNAs_in_tumor_ordered->{$ff}},\@all_tumors)). "\t". join(",", sort {$a cmp $b} @{$miRNAs_Up_in_tumor_ordered->{$ff}})."\t". join(",", sort {$a cmp $b} @{$miRNAs_Down_in_tumor_ordered->{$ff}})."\n"; 
		}
		close OUT_MIR;
	
		my $pairs_in_tumor_ordered;
		my $pairs_UpDown_in_tumor_ordered;
		my $pairs_DownUp_in_tumor_ordered;
	
		foreach my $pairs (keys %$pairs_in_tumor){
			my @tumores=();
			my @tumoresUpDown=();
			my @tumoresDownUp=();
			#my @present_tumor=(uniq (keys %{$prueba->{$GeneName}}));
			foreach my $tumor (keys %{$pairs_in_tumor->{$pairs}}){
				push(@tumores,$tumor);
				push(@tumoresUpDown,$tumor) if(exists($pairs_trUp_miRDw->{$tumor}->{$pairs}));
				push(@tumoresDownUp,$tumor) if(exists($pairs_trDw_miRUp->{$tumor}->{$pairs}));
			}
			@tumores=uniq(@tumores);
			$pairs_in_tumor_ordered->{$pairs}=\@tumores;
			$pairs_UpDown_in_tumor_ordered->{$pairs}=\@tumoresUpDown;
			$pairs_DownUp_in_tumor_ordered->{$pairs}=\@tumoresDownUp;
		}
	
		open(OUT_PAIR,">src/DE_data_TCGA/TCGA_pairs_DE_$type_used\_pathway.xls") || die $!;
		print OUT_PAIR "miRNAs\tAgreement\tnDE\tnDE GeneUp-miRNADown\tnDE GeneDown-miRNAUp\tTumors no DE\tTumors GeneUp-miRNADown\tTumors GeneDown-miRNAUp\n";

		foreach my $ff (sort {@{$pairs_in_tumor_ordered->{$b}} <=> @{$pairs_in_tumor_ordered->{$a}}} keys  %$pairs_in_tumor_ordered){
			print OUT_PAIR $ff ."\t" . $common_pairs->{$ff} ."\t". @{$pairs_in_tumor_ordered->{$ff}} ."\t". @{$pairs_UpDown_in_tumor_ordered->{$ff}} . "\t". @{$pairs_DownUp_in_tumor_ordered->{$ff}} ."\t". join(",",sort {$a cmp $b} intersec(\@{$pairs_in_tumor_ordered->{$ff}},\@all_tumors)). "\t". join(",", sort {$a cmp $b} @{$pairs_UpDown_in_tumor_ordered->{$ff}})."\t". join(",", sort {$a cmp $b} @{$pairs_DownUp_in_tumor_ordered->{$ff}})."\n"; 
		}
		close OUT_PAIR;
		
	}
	elsif(lc($type) eq "de"){
		my $count_file=1;
		
		opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
		my @dires= readdir DIR;
		@dires = sort(@dires);
		#print STDERR "Reading " . (scalar(@files)-2) ." files\n";
		foreach my $tumor_folder (@dires){
			#print "Reading $tumor_folder\n";
			if($tumor_folder !~ /^\./ and $tumor_folder ne "README.md"){
				my $saved;
				opendir ( DIR2, $dirname ."/".$tumor_folder ) || die "Error in opening dir $dirname/$tumor_folder\n";
				my @files= readdir DIR2;
				@files = sort(@files);	
		
				foreach my $filename (@files){
					#Filtering for xls (miRgate output format)
				     if($filename =~ /\.xls$/ and $filename =~ / $data/){
				 
						my $tumor=$filename;
						#Parding filename to gather the tumor name
						$tumor=~s/ RNASeq.xls$//g;
						$tumor=~s/ miRNASeq.xls$//g;
						my $selected_type=$2;
						my $read_file=("$dirname/$tumor_folder/$filename");
						#print STDERR "Reading $filename - $read_file\n";				
						push(@all_tumors,$tumor);
						open(RESULTS, "$read_file") || die "$read_file : $!";
						print STDERR "$count_file) ".date." Reading $tumor data\n";
						$count_file++;
						while(<RESULTS>){
							chomp;
							$_=~s/\"//g;
							if($_ !~/^genes/){
								my ($symbol,$tr_FC,undef,undef,$tr_FDR)=split(/\t/,$_);
								if($tr_FDR<=$FDR_cutoff){
									$common_genes->{$symbol}++;
									$GenesDE->{$tumor}->{$symbol}++;
									$gen_per_tumor->{$symbol}->{$selected_type}++;
							
									if(!exists $genes_in_tumor->{$symbol}->{$tumor}){
										$genes_in_tumor->{$symbol}->{$tumor}++;
									}
				
									if($tr_FC>=$FC_cutoff){
										$common_genesUP->{$symbol}++;
										$GenesUp->{$tumor}->{$symbol}++;
									}
									elsif($tr_FC<= $FC_cutoff_neg ){
										$common_genesDown->{$symbol}++;
										$GenesDown->{$tumor}->{$symbol}++;
									}
								}
							}
						}
						close RESULTS;
					}
				}
			}
		}
		close DIR;
		close DIR2;
		my $cont=0;
		my $genes_in_tumor_ordered;
		my $genes_Up_in_tumor_ordered;
		my $genes_Down_in_tumor_ordered;
		my $output;
		my $header;
		if(lc($data) eq "mirnaseq"){
			$output="src/DE_data_TCGA/TCGA_miRNAs_DE_Global.xls";
			$header="miRNA\tnDE\tnDE Up\tnDE Down\tTumors no DE\tTumors UP\tTumors Down\n";
			foreach my $GeneName (keys %{$common_genes}){
				my @tumores=();
				my @tumoresUp=();
				my @tumoresDown=();
				#my @present_tumor=(uniq (keys %{$prueba->{$GeneName}}));
				foreach my $tumor (keys %{$genes_in_tumor->{$GeneName}}){
					push(@tumores,$tumor);
					push(@tumoresUp,$tumor) if(exists($GenesUp->{$tumor}->{$GeneName}));
					push(@tumoresDown,$tumor)  if(exists($GenesDown->{$tumor}->{$GeneName}));
				}
				@tumores=uniq(@tumores);
				$genes_in_tumor_ordered->{$GeneName}=\@tumores;
				$genes_Up_in_tumor_ordered->{$GeneName}=\@tumoresUp;
				$genes_Down_in_tumor_ordered->{$GeneName}=\@tumoresDown;
			}
			system("mkdir -p src/DE_data_TCGA");
			open(OUT_GENES,">$output") || die $!;
			print OUT_GENES "$header";
		
			foreach my $ff (sort {@{$genes_in_tumor_ordered->{$b}} <=> @{$genes_in_tumor_ordered->{$a}}} keys  %$genes_in_tumor_ordered){
				#last if($cont==101 and $type_used eq "AllGenes");
				$cont++;
				if((@{$genes_in_tumor_ordered->{$ff}})>0){
					print OUT_GENES "$ff\t" . uniq(@{$genes_in_tumor_ordered->{$ff}}) ."\t". uniq(@{$genes_Up_in_tumor_ordered->{$ff}}) . "\t". uniq(@{$genes_Down_in_tumor_ordered->{$ff}}) ."\t". join(",", sort {$a cmp $b} intersec(\@{$genes_in_tumor_ordered->{$ff}},\@all_tumors)). "\t". join(",", sort {$a cmp $b} uniq(@{$genes_Up_in_tumor_ordered->{$ff}}))."\t". join(",", sort {$a cmp $b} uniq(@{$genes_Down_in_tumor_ordered->{$ff}}))."\n"; 
				}
			}
			close OUT_GENES;
		}
		else{
			$output="src/DE_data_TCGA/TCGA_genes_DE_Global.xls";
			$header="GeneName\tPathway\tnDE\tnDE Up\tnDE Down\tTumors no DE\tTumors UP\tTumors Down\n";
			foreach my $GeneName (keys %{$all_posible_genes}){
				my @tumores=();
				my @tumoresUp=();
				my @tumoresDown=();
				#my @present_tumor=(uniq (keys %{$prueba->{$GeneName}}));
				foreach my $tumor (keys %{$genes_in_tumor->{$GeneName}}){
					push(@tumores,$tumor);
					push(@tumoresUp,$tumor) if(exists($GenesUp->{$tumor}->{$GeneName}));
					push(@tumoresDown,$tumor)  if(exists($GenesDown->{$tumor}->{$GeneName}));
				}
				@tumores=uniq(@tumores);
				$genes_in_tumor_ordered->{$GeneName}=\@tumores;
				$genes_Up_in_tumor_ordered->{$GeneName}=\@tumoresUp;
				$genes_Down_in_tumor_ordered->{$GeneName}=\@tumoresDown;
			}
			
			system("mkdir -p src/DE_data_TCGA");
			open(OUT_GENES,">$output") || die $!;
			print OUT_GENES "$header";
		
			foreach my $ff (sort {@{$genes_in_tumor_ordered->{$b}} <=> @{$genes_in_tumor_ordered->{$a}}} keys  %$genes_in_tumor_ordered){
				#last if($cont==101 and $type_used eq "AllGenes");
				$cont++;
				if((@{$genes_in_tumor_ordered->{$ff}})>0){
					print OUT_GENES "$ff\t" .join(",", keys %{$all_posible_genes->{$ff}}) ."\t". uniq(@{$genes_in_tumor_ordered->{$ff}}) ."\t". uniq(@{$genes_Up_in_tumor_ordered->{$ff}}) . "\t". uniq(@{$genes_Down_in_tumor_ordered->{$ff}}) ."\t". join(",", sort {$a cmp $b} intersec(\@{$genes_in_tumor_ordered->{$ff}},\@all_tumors)). "\t". join(",", sort {$a cmp $b} uniq(@{$genes_Up_in_tumor_ordered->{$ff}}))."\t". join(",", sort {$a cmp $b} uniq(@{$genes_Down_in_tumor_ordered->{$ff}}))."\n"; 
				}
			}
			close OUT_GENES;
		}
	}
}
elsif($help){
	help();
}
else{
	help();
}
sub help{
        my $usage = qq{
          $0 

            Getting help:
                [--help]

            Needed parameters:
              [dir|d] : Folder with xls files from miRGate             
              [type|t] : Type of analysis: pairs (for miRGate) o DE (for constitutive DE study)
			Optional parameters;
              [net|t] : Network to explore (none for all genes)
			
            Examples:
              perl $0 -d miRGate_Results/ -t pairs -net DDR
                       
};

print STDERR $usage;
exit();
        
}

sub date{
        my $date=("date \"+%D %H:%M:%S\"");
        $date=`$date`;
        chomp $date;
        return("[".$date."]");
}

sub intersec{
	use Array::Utils qw(:all);
	my (@found)=@{$_[0]};
	my (@all)=@{$_[1]};
	
	#Not DE for this gene
	my @minus = array_minus( @all, @found );
	return(uniq(@minus));
}
__END__

head1 NAME

3.Summary_miRGate - Perl Module to ..

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 3.Summary_miRGate

=over 2

=item a

bla bla1

=item b

bla bla2

=back

=head1 EXAMPLE
