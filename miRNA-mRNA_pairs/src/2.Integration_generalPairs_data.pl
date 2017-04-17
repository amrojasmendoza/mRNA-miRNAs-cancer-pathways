#
#  Integrative_pairs_all_tumors.pl
#
#  Created by Eduardo Andrés León on 2016-06-25.
#  Copyright (c) 2016 IPBLN. All rights reserved.
#
#!/usr/bin/perl
$|=1;
use strict;
#use warnings;
use Getopt::Long;
use List::MoreUtils qw/ uniq /;

my $limit=8;
my $miRGate_agreement=2;
my $dirname="src/DE_data_TCGA/";
my $help;
GetOptions(
	"help" => \$help,
	"dir|d=s" => \$dirname,
	"agree|a=s" => \$miRGate_agreement,
	"limit|l=s" => \$limit,
);
if($dirname and $limit){
	opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
	my @files= readdir DIR;
	@files = sort(@files);
	my $data;
	my $all_data;
	foreach my $filename (@files){
		#Filtering for xls (miRgate output format)
	     if($filename =~ /\.xls$/ and $filename =~/^TCGA_pairs_DE_/){
			 open(FILE,$dirname."/".$filename) || die "Cant open $filename\n";
			 while(<FILE>){
				 chomp;
				 my($pairs,$agree,$nDE,$nUp,$nDown)=split(/\t/);
				 my($gene,$miR)=split("_",$pairs);
				 if($agree >=$miRGate_agreement and $nDE >=$limit){
					 $data->{"$gene,$miR,$agree,$nDE,$nUp,$nDown"}++;
				 }
				 $all_data->{$_}++;
			 }
			 close FILE;
		 }
	 }
	 close DIR;
	 
	 my $output="general.pairs.csv";
	 open(OUT,">$output") || die $!;
	 print OUT "GeneName,miRNA,Agreement,nDE,nUp,nDown\n";
	 foreach my $info(sort keys %{$data}){
		 print OUT $info ."\n";
	 }
	 close OUT;
	 
	 my $output2="src/DE_data_TCGA/TCGA_pairs_DE_pathway_agreement2.xls";
	 open(OUT2,">$output2") || die $!;
	 #print OUT2 "GeneName,miRNA,Agreement,nDE,nUp,nDown\n";
	 foreach my $info(sort keys %{$all_data}){
		 print OUT2 $info ."\n";
	 }
	 close OUT2;
	 
}
sub date{
        my $date=("date \"+%D %H:%M:%S\"");
        $date=`$date`;
        chomp $date;
        return("[".$date."]");
}
