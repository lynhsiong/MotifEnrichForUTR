
# the input files in the attached config file (you must edit the file "sample.config" to your own config file)
#
# file1: results.0dot1.out
#RBP009_human_eclip_bed3_mw30_num10_loose        ENST00000262291.9       1615    1625    +       52.2674515236865        4.92445303860747
#
# file2: MIIGV.longest.uniqname.lst
#A4GALT  ENST00000401850.5       ENSG00000128274.17      GV      7227
#A4GALT  ENST00000401850.5       ENSG00000128274.17      MII     12157

BEGIN{
unshift @INC, "/perl5/lib/perl5/lib/perl5";
unshift @INC, "/Statistics-Lite-3.62/blib/lib";
} ### modules used in this script, specify according to your situation

use List::Util qw/shuffle/;
use Statistics::Lite qw/mean stddev/;
use Data::Dumper;

$cfg=shift;
$out=shift;
open OUT, ">$out";

%cfg=();
open CFG, "<$cfg";
while(<CFG>){
next if(/^#/);
@cfg=split("\t", $_);
$cfg{$cfg[0]}=$cfg[1];
}
close CFG;

$lib=$cfg{"LIB"};
#$list=shift;

#=head

%hash=();
open LIB, "<$lib";
while(<LIB>){
@l=split("\t", $_);
$hash{$l[1]}{$l[0]}++;
}
close LIB;

#=cut

%h=%id=();
for $file (keys %cfg){
  if($file!~/LIB/){
  $Fdir=$cfg{$file};
  ($h, $id)=readfile($file, $Fdir, \%h, \%id);
  %h=%{$h};
  %id=%{$id};
  }
}


#print OUT Dumper(%h);
#print OUT "================================";
#print OUT Dumper(%id);
#print OUT "-----------===============================";


%sample=();
%rand_tmph=%{$h{"RANDOM"}};
@rand_tmph=keys %rand_tmph;
foreach $key (keys %h){
  
  if($key!~/RANDOM/){
   @j=();
   $tmph=$h{$key};
   @j = keys(%$tmph);
   @{$sample{$key}{"head"}}=@j;

  $oc=scalar(@j);
  for(1..100){
  @k=shuffle(@rand_tmph);
  $seed=int(rand($#k-$oc));
#print OUT scalar(@rand_tmph), " --- ", $#k, " - ", $oc,"\n";
  @{$sample{$key}{"rand"}{$_}}=@k[$seed..($seed+$oc-1)];
  }

  }
}

#print OUT Dumper(%sample);

# Motif_ID	p_values(two_tail)	Motif_occurrence_num	mean and stdev in random(an array)	SeqID1,SeqID2...(Head/tail;sorted_by_number)

foreach $MIIGV (keys %sample){

%sample_num=%{$sample{$MIIGV}};

%rand_motif=();

@head=@{$sample_num{"head"}};
%rand=%{$sample_num{"rand"}};

foreach $rand_n (sort {$a<=>$b} keys %rand){
  for(@{$rand{$rand_n}}){
    if(defined($hash{$_})){
       %gene_hash=%{$hash{$_}};
       foreach $m (keys %gene_hash){
           $rand_motif{$m}{$rand_n}+=$gene_hash{$m};
       }
    }
  }
}

%rand_motif_stat=();
for $stat_m (keys %rand_motif){
@stat_values=values %{$rand_motif{$stat_m}};
$rand_motif_stat{$stat_m}{"mean"}=mean(@stat_values);
$rand_motif_stat{$stat_m}{"stddev"}=stddev(@stat_values);
}

printstat(\@head, \%rand_motif_stat, $MIIGV);

}

close OUT;





sub printstat{

#$label=shift;
#$s_n=shift;
$arr_toTest=shift; @arr_toTest=@{$arr_toTest};
$rand_bkg=shift; %rand_bkg=%{$rand_bkg};
$phase=shift;

%tmph=%motif=();
%gene_pool=();
for(@arr_toTest){
$tmp_gene_id=$_;
if(defined($hash{$tmp_gene_id})){
  %tmph=%{$hash{$tmp_gene_id}};
  foreach $mk (keys %tmph){

    if(defined($id{$tmp_gene_id})){
      $gene_pool{$mk}{$id{$tmp_gene_id}}=$tmph{$mk};
    }

    $motif{$mk}+=$tmph{$mk};

  }
}
}

for $fn_k (keys %motif){
  if(defined($rand_bkg{$fn_k})){
  $fn_m=$rand_bkg{$fn_k}{"mean"}; $fn_sd=$rand_bkg{$fn_k}{"stddev"};
  $fn_p=ttest($motif{$fn_k}, $fn_m, $fn_sd);

  %fn_gene_sort=(); @fn_gene_sort=();
  if(defined($gene_pool{$fn_k})){
     %fn_gene_sort=%{$gene_pool{$fn_k}};
     @fn_gene_sort=sort { $fn_gene_sort{$b} <=> $fn_gene_sort{$a} } keys %fn_gene_sort;
  }
  
  print OUT $phase, "\t", $fn_k, "\t", $fn_p, "\t", $motif{$fn_k}, "\t", $fn_m, ":", $fn_sd, "\t", join(",", @fn_gene_sort), "\n";
  }
}

}

sub ttest{

   $_p=shift;
   $_m=shift;
   $_sd=shift;

   $pi = 3.14159265358979323;
   $_v=$_sd**2;

   $_v=($_v<=0)?0.0000000000000001:$_v;
   
   $_denom=(2*$pi*$_v)**0.5;
   $_num=exp(-1*($_p-$_m)**2/(2*$_v));

   return sprintf("%.16f", $_num/$_denom);

}



sub readfile{

$fname=shift;
$fdir=shift;
$_tmph=shift; %_tmph=%{$_tmph};
$_id=shift; %_id=%{$_id};

open TMP, "<$fdir";
while(<TMP>){
next if(/^#/);

@lll=split("\t", $_);
$_tmph{$fname}{$lll[1]}=$lll[0];
$_id{$lll[1]}=$lll[0];

}
close TMP;

return \%_tmph, \%_id;

}
  





