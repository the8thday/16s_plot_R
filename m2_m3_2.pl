#!/usr/bin/perl -w
#usage: perl m2_m3.pl <m2> <m3>
open FH1,"< $ARGV[0]";
open FH2,"> temp.txt";
my $line1 = <FH1>;
chomp $line1;
my @sample = split "\t",$line1;
pop @sample;
shift @sample;
print FH2 "",(join "\t",@sample),"\ttaxonomy\n";
while(<FH1>){
	chomp;
	my @temp = split "\t",$_;
	shift @temp;
	my $taxonomy = pop @temp;
	my $abundance = join "\t",@temp;
	my @temp2 = split "; ",$taxonomy;
	if($taxonomy =~ m/^(k__.*?); (p__.*?); (c__.*?); (o__.*?); (f__.*?); (g__.*?); (s__.*?)$/){
		my $k = $1;
		my $p = "$1; $2";
		my $c = "$1; $2; $3";
		my $o = "$1; $2; $3; $4";
		my $f = "$1; $2; $3; $4; $5";
		my $g = "$1; $2; $3; $4; $5; $6";
		my $s = "$1; $2; $3; $4; $5; $6; $7";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
		print FH2 "$abundance	","$c\n";	
		print FH2 "$abundance	","$o\n";	
		print FH2 "$abundance	","$f\n";	
		print FH2 "$abundance	","$g\n";	
		print FH2 "$abundance	","$s\n";	
	}
	if($taxonomy =~ m/^(k__.*?); (p__.*?); (c__.*?); (o__.*?); (f__.*?); (g__.*?)$/ && ( ! ($taxonomy =~ m/s__/) ) ){
		my $k = $1;
		my $p = "$1; $2";
		my $c = "$1; $2; $3";
		my $o = "$1; $2; $3; $4";
		my $f = "$1; $2; $3; $4; $5";
		my $g = "$1; $2; $3; $4; $5; $6";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
		print FH2 "$abundance	","$c\n";	
		print FH2 "$abundance	","$o\n";	
		print FH2 "$abundance	","$f\n";	
		print FH2 "$abundance	","$g\n";	
}
	if($taxonomy =~ m/^(k__.*?); (p__.*?); (c__.*?); (o__.*?); (f__.*?)$/  && ( ! ($taxonomy =~ m/g__/) ) ){
		my $k = $1;
		my $p = "$1; $2";
		my $c = "$1; $2; $3";
		my $o = "$1; $2; $3; $4";
		my $f = "$1; $2; $3; $4; $5";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
		print FH2 "$abundance	","$c\n";	
		print FH2 "$abundance	","$o\n";	
		print FH2 "$abundance	","$f\n";	
}
	if($taxonomy =~ m/^(k__.*?); (p__.*?); (c__.*?); (o__.*?)$/  && ( ! ($taxonomy =~ m/f__/) ) ){
		my $k = $1;
		my $p = "$1; $2";
		my $c = "$1; $2; $3";
		my $o = "$1; $2; $3; $4";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
		print FH2 "$abundance	","$c\n";	
		print FH2 "$abundance	","$o\n";	
}
	if($taxonomy =~ m/^(k__.*?); (p__.*?); (c__.*?)$/  && ( ! ($taxonomy =~ m/o__/) ) ){
		my $k = $1;
		my $p = "$1; $2";
		my $c = "$1; $2; $3";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
		print FH2 "$abundance	","$c\n";	
}
	if($taxonomy =~ m/^(k__.*?); (p__.*?)$/  && ( ! ($taxonomy =~ m/c__/) ) ){
		my $k = $1;
		my $p = "$1; $2";
		print FH2 "$abundance	","$k\n";	
		print FH2 "$abundance	","$p\n";	
}
	if($taxonomy =~ m/^(k__.*?)$/  && ( ! ($taxonomy =~ m/p__/) ) ){
		my $k = $1;
		print FH2 "$abundance	","$k\n";	
}
}
close FH1;
close FH2;

open FH3,"< temp.txt";
my $head = <FH3>;
chomp $head;
my @temp3 = split "\t",$head;
pop @temp3;
my %hash;
while(<FH3>){
	chomp;
	my @temp = split "\t",$_;
	my $taxonomy = pop @temp;
	#print $temp3[0],"\n";
	#print "$taxonomy\n";
	foreach(0..$#temp3){
		if(exists $hash{$taxonomy}{$temp3[$_]}){
			#print "$temp[$_]\n";
			$hash{$taxonomy}{$temp3[$_]} += $temp[$_];
		}else{
			$hash{$taxonomy}{$temp3[$_]} = $temp[$_];
		}
	}
	
}
close FH3;


open FH4,"> temp.txt";
my @sort_temp3 = sort @temp3;
print FH4 "OTU ID\t",(join "\t",@sort_temp3),"\ttaxonomy\n";

foreach my $taxonomy (sort keys %hash ){
	my @temp;
	foreach my $sample(sort @temp3){
		push(@temp,$hash{$taxonomy}{$sample});
	}
	my $tax = (split "; ",$taxonomy)[-1];
	print FH4 "$tax\t",(join "\t",@temp),"\t$taxonomy\n";
	}
close FH4;


open FH5,"< temp.txt";
open FH6,"> $ARGV[1]";
print FH6 "",scalar <FH5>;
while(<FH5>){
        chomp;
        my @temp = split "\t",$_;
        my $index = $#temp;
        my $taxonomy_short = (split "; ",$temp[0])[-1];
        my @temp3 = @temp[1..$index];
        unshift(@temp3,$taxonomy_short);
        print FH6 "",(join "\t",@temp3),"\n";
}
close FH5;
close FH6;
system("rm temp.txt");
