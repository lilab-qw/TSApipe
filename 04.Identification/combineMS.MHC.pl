#!/usr/bin/perl -w

my ($mfile,$file1,$ofile) = @ARGV;

my (%hash1);

open IN1, "$file1" or die "$!\n";
my $head1 = <IN1>;
chomp $head1;
while(<IN1>){
        chomp;
        my @temp = split /\t/, $_;
        my $pep = shift @temp;
        my $value = join("\t", @temp);
        $hash1{$pep} = $value;
}
close IN1;


open OUT, ">$ofile";
open IN, "$mfile" or die "$!\n";
my $head = <IN>;
chomp $head;
print OUT "ID\t$head\t$head1\n";
my $id = 1;
while(<IN>){
        chomp;
        my @temp = split /\t/,$_;
        my $pep = $temp[3];
        my ($mhci);
        if(exists $hash1{$pep}){$mhci = $hash1{$pep}}else{$mhci = "-\t-\t-\t-\t-\t-\t-"}
        print OUT "ID$id\t$_\t$mhci\n";
        $id += 1;

}
close IN;
close OUT;
