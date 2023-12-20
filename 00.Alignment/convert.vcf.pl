#!/usr/bin/perl -w

my ($sample,$infile, $outfile, $mfile) = @ARGV;

open OUT, ">$outfile";
print OUT "chromosomeNumber\tuniqueId\tstart\tend\tref\talleles\tquality\tcaller\n";
my $tag = 1;
open IN, "$infile" or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/);
        my ($chr,$s,$ref,$alt,$qual) = (split)[0,1,3,4,5];
	my $e = $s + length($ref);
	next if($qual < 1);
	$chr =~ s/chr//;
	if($chr =~ /M/){$chr ="MT"}
        print OUT "$chr\t$tag\t$s\t$e\t$ref\t$alt\t$qual\tfreebayes\n";
	$tag += 1;
}
close IN;
my $file = (split /\//,$outfile)[-1];
open OO, ">$mfile";
print OO "[package_infos]\ndescription = SNPs\nmaintainer = wmy\nmaintainer_contact = wmy.cs\nversion = 1\n\n[set_infos]\nspecies = human\nname = $sample\ntype = Agnostic\nsource = freebayesHG38\n\n[snps]\nfilename = $file\n";
close OUT;
close OO;
