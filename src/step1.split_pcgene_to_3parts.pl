$in=shift; ### gencode style fasta file for protein-coding gene

%hash_fa=%hash_pos=();
open IN, "<$in";
while(<IN>){
  chomp;
  $ln=$_;
  if($ln=~/^>(\S+)(\|ENSMUSG.+)/){    ### human (ENSG) or mouse prefix
    $id=$1;
    $left=$2;
    if($left=~/.?UTR5:(\d+)-(\d+).?/){
       @{$hash_pos{$id}{"UTR5"}}=@{[$1,$2]};
    }
    if($left=~/.?CDS:(\d+)-(\d+).?/){
       @{$hash_pos{$id}{"CDS"}}=@{[$1,$2]};
    }
    if($left=~/.?UTR3:(\d+)-(\d+).?/){
       @{$hash_pos{$id}{"UTR3"}}=@{[$1,$2]};
    }    
  }else{
    $hash_fa{$id}.=$ln;
  }
}

$outtype="UTR3";  ### UTR5 or UTR3 or CDS
foreach $key (keys %hash_fa){
  if(defined $hash_pos{$key}{$outtype}){
     @tmppos=@{$hash_pos{$key}{$outtype}};
     $seq=substr($hash_fa{$key}, $tmppos[0]-1, $tmppos[1]-$tmppos[0]+1);
     print ">",$key, "\n", $seq, "\n";
  }
}
