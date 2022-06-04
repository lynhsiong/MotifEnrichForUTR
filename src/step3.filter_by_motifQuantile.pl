$Quantile=0.1;    ### only keep first Quantile (0.1: top 10%) hits for all fimo sites of a special motif id 

$some_dir="/PATH/to/fimo_results";  ### edit it for your run (output of step2)

$prefix="fimo.SPLIT";  ### (prefix used in your step2)

%hash=();
@line=();

$some_dir_file=$some_dir."/fimo.SPLIT*.out";
if (glob("$some_dir_file")){
opendir($prefix, $some_dir);
@_3rd_cluster = grep { !/^\./ && /^$prefix.+\.out$/ } readdir($prefix);

for ( @_3rd_cluster ){ s/^/$some_dir\//;

$myarr=readfile($_);
@line=(@line, @{$myarr});

}

closedir $prefix;
}



for (@line){

@ln=split("\t", $_);

$randnum=rand(100);
$id=$ln[1]."_".$ln[2]."_".$ln[3]."_".$ln[4]."_".$randnum;
#print $id,"    >===<\n";

$hash{$ln[0]}{$id}=-1*(log($ln[6])/log(10));

}


#$sum=0;
#for $x (keys %hash){
#@y=keys %{$hash{$x}};
#$sum+=scalar(@y);
#}
#print $sum,"   <<<\n";


#$ccc=0;
#%ccc=();
foreach $motif (keys %hash){
   $h=$hash{$motif};
   @k=@v=();
   @k = sort { $h->{$b} <=> $h->{$a} } keys(%$h);
   @v = @{$h}{@k};
#$ccc+=scalar(@k);
#print $motif,"   ",scalar(@k),"   <===\n";
   for(0..int($Quantile*scalar(@k))){
   #for(0..$#k){
       print $motif,"\t",join("\t", @{[split("_", $k[$_])]}),"\t",$v[$_],"\n";

#       $ccc{${[split("_", $k[$_])]}[0]}=1;

   }
   
}

#print scalar(keys %ccc),"\t$ccc\t\tnumber\n";

sub readfile{

$fname=shift;

@ln=();
open TMP, "<$fname";
while(<TMP>){
next if(/^#/);

push @ln, $_;

}
close TMP;

return \@ln;

}
