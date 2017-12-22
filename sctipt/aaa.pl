#! /usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use threads;
use threads::shared;
my $help=<<USAGE;
USAGE:
-i1: Undetermined_R1        [force]
-i2: Undetermined_R2        [option]
-o:  outdir
-s:  sample and index       [sample1:index1,sample2:index2,.....] 
-m:  mismatch               [default:0]      
USAGE
my ($i1,$o,$i2,$s,$m);
GetOptions
(
"i1:s"  =>\$i1,
"o:s"  =>\$o,
"i2:s"  =>\$i2,
"m:s"  =>\$m,
"s:s"  =>\$s
);
die "$help" unless( defined  $i1 and $s and $o);
$m ||= 0;
my $max_treads=10;
my $num_thread = 0;
my %multi=();
my @a=split/\,/,$s;
foreach my $key (@a){
    my @b=split/\:/,$key;
    if($num_thread < $max_treads){
        $multi{$key}= threads->create(\&abstract,$b[0],$b[1],$o,$i1,$i2,$m);
        $num_thread ++;
    }else{
        while(1){
            my $break = 0;
            foreach my $file (keys %multi){
                if($multi{$file} ->is_joinable()){
                    my $time = $multi{$file} ->join();
                    print STDERR "$file Abstract done.\n";
                    $num_thread --;
                    delete $multi{$file};
                    $break = 1;
                    last;                        				
                }
            }
            last if($break == 1);
            sleep 10;
        }
        redo; 		
    }
}
while(%multi){
    foreach my $file (keys %multi){
        if($multi{$file} ->is_joinable()){
            my $time = $multi{$file} ->join();
            print STDERR "$file Abstract done.\n";
            delete $multi{$file};
        }
    }
}


# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub abstract{
    my ($prefix,$barcode,$out,$Un1,$Un2,$mistach)=@_;
    my $left=$prefix."_R1_001.fastq.gz";
    my $right=$prefix."_R2_001.fastq.gz";	
    my %pair=();
    $pair{$left}= threads->create(\&abstract_single,$barcode,$left,$Un1,$mistach);
    $pair{$right}= threads->create(\&abstract_single,$barcode,$right,$Un2,$mistach);	
    while(%pair){
        foreach my $file (keys %pair){
            if($pair{$file} ->is_joinable()){
               my $time = $pair{$file} ->join();
               print STDERR "$time\n";
               delete $pair{$file};
            }
        }
    }
}	

sub abstract_single{
    my ($barcode,$out,$in,$max_mistach)=@_;	
    open (OUT,"| gzip >$out") or die $!;
    if ($in =~ /gz$/){
        open (IN,"gzip -dc $in|") or die $!;
    }else{
        open (IN,$in) or die $!;
    }    	
    while (<IN>){
        chomp;
        s/\r//;
        next if (/^$/);
        my @tmp=split /\s+/,$_;
        my $index=(split /\:/,$tmp[1])[3];
        my $mismatch=&mismatch($index,$barcode);		
        if ($mismatch <= $max_mistach){
            print OUT $_."\n".<IN>.<IN>.<IN>;
        } 			
    }
    close IN;
    close OUT;
}	
	
sub mismatch{
    my $in1=shift;
    my $in2=shift;
    my $out=0;
    my @tmp1=split //,$in1;
    my @tmp2=split //,$in2;
    for (my $i=0;$i<@tmp1;$i++){
        if ($tmp1[$i] eq $tmp2[$i]){
            $out++;
        }
    }
    my $mismatch=scalar(@tmp1)-$out;
    return $mismatch;
}



