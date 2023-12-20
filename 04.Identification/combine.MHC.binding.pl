#!/usr/bin/perl -w
use List::Util;

my ($list,$ofile) = @ARGV;
my %hash;
open LI, "$list" or die "$!\n";
while(my $file=<LI>){
        chomp $file;
        open IN, "$file" or die "$!\n";
        while(<IN>){
                chomp;
                s/^\s+//g;
                next if(/^Protein/);
                if(/PEPLIST/){
                        my($hla, $pep, $aff, $rank) = (split /\s+/, $_)[1,2,12,13];
                        my $value = "$hla\t$aff\t$rank";
                        push @{$hash{$pep}}, $value;

                }
        }
        close IN;

}
close LI;

open OUT, ">$ofile";
print OUT "Peptide\tHLA_A_aff\tHLA_B_aff\tHLA_C_aff\tMin_aff\tMin_aff_allele\tMin_rank\tMin_rank_allele\n";
foreach my $key(sort keys %hash){
        my @arrs = @{$hash{$key}};
        my ($af,$bf,$cf,$mf,$mfa,$mr,$mra);
        my (%hash1,%hash2,@affs,@ranks,@A,@B,@C);
        foreach my $arr(@arrs){
                my ($hla,$aff,$rank) = (split /\t/, $arr)[0,1,2];
                push @affs, $aff;
                push @ranks, $rank;
                $hash1{$hla} = $aff;
                $hash2{$hla} = $rank;
                if($hla =~ /^HLA-A/){push @A,$hla}
                if($hla =~ /^HLA-B/){push @B,$hla}
                if($hla =~ /^HLA-C/){push @C,$hla}
                #print "$arr\t";
        }
        @affs = sort { $a <=> $b }(@affs);
        @ranks = sort { $a <=> $b } (@ranks);
        $mf = $affs[0];
        $mr = $ranks[0];
        foreach my $k(sort keys %hash1){if($hash1{$k} == $mf){$mfa = $k}}
        foreach my $k2(sort keys %hash2){if($hash2{$k2} == $mr){$mra = $k2}}
        #print "$mfa:$mf\t$mra:$mr\n";
        if($#A == 0){$af = $hash1{$A[0]}}
        else{
                if($hash1{$A[0]} < $hash1{$A[1]}){$af = $hash1{$A[0]}}
                else{$af = $hash1{$A[1]}}
        }
        if($#B == 0){$bf = $hash1{$B[0]}}
        else{
                if($hash1{$B[0]} < $hash1{$B[1]}){$bf = $hash1{$B[0]}}
                else{$bf = $hash1{$B[1]}}
        }
        if($#C == 0){$cf = $hash1{$C[0]}}
        else{
                if($hash1{$C[0]} < $hash1{$C[1]}){$cf = $hash1{$C[0]}}
                else{$cf = $hash1{$C[1]}}
        }
        print OUT "$key\t$af\t$bf\t$cf\t$mf\t$mfa\t$mr\t$mra\n";

}
close OUT;
