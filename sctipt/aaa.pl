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
-m:  max mismatch           [default:0]
-p： max threads            [default:20]     
USAGE
my ($i1,$o,$i2,$s,$m,$max_treads);
GetOptions
(
"i1:s"  =>\$i1,
"o:s"  =>\$o,
"i2:s"  =>\$i2,
"m:s"  =>\$m,
"p:s"  =>\$max_treads,
"s:s"  =>\$s
);
die "$help" unless( defined  $i1 and $s and $o);
$m ||= 0;
$max_treads||=20;
my $num_thread = 0;
my %multi=();
my @a=split/\,/,$s;
foreach my $key (@a){
    my @b=split/\:/,$key;
    if($num_thread < $max_treads){
        my $left=$b[0]."left";
        my $o1=$b[0]."_R1_001.fastq.gz";		
        $multi{$left}= threads->create(\&abstract_single,$b[1],"$o/$o1",$i1,$m);
        $num_thread ++;
        if (defined $i2){
            my $right=$b[0]."right";		 
            my $o2=$b[0]."_R2_001.fastq.gz";		
            $multi{$right}= threads->create(\&abstract_single,$b[1],"$o/$o2",$i2,$m);		
            $num_thread ++;		
        }		
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
        my $first=$_;		
        my @tmp=split /\s+/,$first;
        my $index=(split /\:/,$tmp[1])[3];
        #print "$index\n$barcode\n";die; 		
        my $m_current=&chongfu($index,$barcode);
        chomp(my $reads=<IN>);
        chomp(my $names=<IN>);
        chomp(my $quality=<IN>);  	
        if ($m_current <= $max_mistach){
            print OUT "$first\n$reads\n$names\n$quality\n";		
        } 			
    }
    close IN;
    close OUT;
}	
	
sub chongfu{
    my ($in1,$in2)=@_;	
    if(length($in1) != length($in2)){
        print length($in1)."\n".length($in2)."\n";	
        die	"the length of two index not equal\n";
    }
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
