#!/usr/bin/perl  -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Math::Complex;
use Cwd;
my $help=<<USAGE;
USAGE:
-i : analysis dir
-o : outdir
-s : Type of DNA. [gDNA|cfDNA|bDNA],default: 'cfDNA'   
USAGE
my ($i,$o,$sampletype);
GetOptions
(
         "i:s"  =>\$i,
         "s:s"  =>\$sampletype,		 
         "o:s"  =>\$o
);
die "$help" unless( defined $i and $o);
$sampletype ||= 'cfDNA';
$sampletype = uc($sampletype);
chomp(our $Script_Home=`dirname $0`);
$Script_Home=&ABSOLUTE_DIR($Script_Home);
my $json="$Script_Home/coverge_web_js/js/echarts.js";
my $json2="$Script_Home/coverge_web_js/js";
$i=&ABSOLUTE_DIR($i);
my $tmp_dir=$i;
$tmp_dir=~s/task_data/fastq_data/;
mkdir $o unless (-d $o);
my $raw_left_qc=glob "$i/trim/fastqc/raw/*_R1_*fastqc.html";
my $raw_right_qc=glob "$i/trim/fastqc/raw/*_R2_*fastqc.html";
my $trim_left_qc= glob "$i/trim/fastqc/trim/*_R1_fastqc.html";
my $trim_right_qc=glob "$i/trim/fastqc/trim/*_R2_fastqc.html";
my $prefix=basename ($trim_left_qc);
$prefix=~m/(.*?)\_R1\_/;
my $s=$1;
my $index=0;
$tmp_dir=~s/$s//;
my $d=$tmp_dir;

my $cmd1="cp -r $d/Reports $o";system $cmd1;
my $cmd2="cp $raw_left_qc  $o/raw_left_qc.html";system $cmd2;
my $cmd3="cp $raw_right_qc  $o/raw_right_qc.html";system $cmd3;
my $cmd4="cp $trim_left_qc  $o/trim_left_qc.html";system $cmd4;
my $cmd5="cp $trim_right_qc  $o/trim_right_qc.html";system $cmd5;
my $cmd6="cp -r $json2 $o";system $cmd6;
my $o_raw_left_qc="raw_left_qc.html";
my $o_raw_right_qc="raw_right_qc.html";
my $o_trim_left_qc= "trim_left_qc.html";
my $o_trim_right_qc="trim_right_qc.html";
my $split_report="Reports/html/index.html";
my $o_json="js/echarts.js";
my $o_json2="js";
#print "$d\n$split_report\n";die;
my @normalized_coverage=();
my @coverge_name=();
my @coverge=();
my @capture=();
my $capture_file="$i/Coverage_result.txt";
my $coverhe_file=glob "$i/bam/Realign/Coverage/*coverage.bed";
open (IN,$capture_file) or die $!;
while (<IN>){
    chomp;
    next if (/^$/);next if(/^\#/);next if (/^Sample/);
    @capture=split /\s+/,$_;
    }
close IN;
my ($more1,$less_1)=(0,0);
my ($depth_50,$depth_300,$depth_500)=(0,0,0);
open (IN,$coverhe_file) or die $!;
while (<IN>){
    chomp;
    next if (/^$/);next if(/^\#/);next if (/^chrom/);
    $index++;
    my @tmp=split /\s+/,$_;
    my $coverge_name="'".$tmp[0].":".$tmp[1]."-"."$tmp[2]"."'";
    push @coverge_name,	$coverge_name;
    push @coverge ,$tmp[6];
    if ($sampletype eq "BDNA" && $tmp[6]<50){$depth_50++;}
    if ($sampletype eq "GDNA" && $tmp[6]<300){$depth_300++;}
    if ($sampletype eq "CFDNA" && $tmp[6]<500){$depth_500++;}	
    my $log;
    if($tmp[7]>0){$log=logn($tmp[7], 2);}else{$log=-2;}	
    if($log>=1){$more1++;}
    if($log<=-1){$less_1++;}	
    my $index_coverge="[".$index.",".$log."]";
    push @normalized_coverage, $index_coverge;
    }
close IN;	
my $data_normalized_coverage=join ",",@normalized_coverage;
my $data_coverge_name=join ",",@coverge_name;
my $data_coverge=join ",",@coverge;
my $more1_per=int(($more1/$index)*10000)/100;
my $less_1_per=int(($less_1/$index)*10000)/100;
my $less_strand;
my $less_strand_abs;
my $strand_depth;
if ($sampletype eq "BDNA"){$less_strand=int(($depth_50/$index)*10000)/100;$less_strand_abs=$depth_50;$strand_depth="50X";}
if ($sampletype eq "GDNA"){$less_strand=int(($depth_300/$index)*10000)/100;$less_strand_abs=$depth_300;$strand_depth="300X";}
if ($sampletype eq "CFDNA"){$less_strand=int(($depth_500/$index)*10000)/100;$less_strand_abs=$depth_500;$strand_depth="500X";}
open (OUT,">$o/$s.threshold.txt")  or die $!;
print OUT "#index\tnumber\tpercent\n";
print OUT "normalized_coverage>1\t$more1\t$more1_per%\n";
print OUT "normalized_coverage<-1\t$less_1\t$less_1_per%\n";
print OUT "lower_than_standard_depth\t$less_strand_abs\t$less_strand%\n";
close OUT;

open (OUT,">$o/$s.html")  or die $!;
print OUT "\<\!DOCTYPE html\>\n\<html lang=\"en\"\>\n\<head\>\n    \<meta charset=\"utf-8\"\>\n    \<title\>$s\<\/title\>\n";
print OUT "\<\/head\>\n\<body\>\n\<h1 align=\"center\"\>$s Statistics Results\<\/h1\>\n\<a href=\"$split_report\"\>bcl2fastq results\<\/a\>\n";
print OUT "\<a href=\"$o_raw_left_qc\" style=\"margin-left:20px\;\"\>Raw R1 FastQC results\<\/a\>\n"; 
print OUT "\<a href=\"$o_raw_right_qc\" style=\"margin-left:20px\;\"\>Raw R2 FastQC results\<\/a\>\n";
print OUT "\<a href=\"$o_trim_left_qc\" style=\"margin-left:20px\;\"\>Trim R1 FastQC results\<\/a\>\n"; 
print OUT "\<a href=\"$o_trim_right_qc\" style=\"margin-left:20px\;\"\>Trim R2 FastQC results\<\/a\>\n";
print OUT "\<div style=\"height:100px\;border:1px solid #ccc;padding:10px\;margin-top:20px\;\"\>\n";
print OUT "\<table border=\"1\" align=\"center\" width=\"100\%\"  style=\"font-family:\'宋体\'\; font-size:12px\; margin-top:5px\;\"\>\n\<tr\>\n";
print OUT "\<th\>Sample\<\/th\>\n\<th\>Target_Size\<\/th\>\n\<th\>Raw_Reads\<\/th\>\n\<th\>Raw_Capture_Rate\<\/th\>\n\<th\>Raw_Target_Coverage\<\/th\>\n";
print OUT "\<th\>Estimate_Percent_Dup\<\/th\>\n\<th\>DeDup_Reads\<\/th\>\n\<th\>DeDup_Capture_Rate\<\/th\>\n\<th\>DeDup_Target_Coverage\<\/th\>\n";
print OUT "\<\/tr\>\n\<tr\>\n\<td\>$capture[0]\<\/td\>\n<td\>$capture[1]\<\/td\>\n<td\>$capture[2]\<\/td\>\n<td\>$capture[3]\<\/td\>\n<td\>$capture[4]\<\/td\>\n";
print OUT "<td\>$capture[5]\<\/td\>\n<td\>$capture[6]\<\/td\>\n<td\>$capture[7]\<\/td\>\n<td\>$capture[8]\<\/td\>\n\<\/tr\>\n\<\/table\>\n\<\/div\>\n";
print OUT "\<div id=\"main\" style=\"height:500px\;border:1px solid #ccc;padding:10px\;\"\>\<\/div\>\n";
print OUT "    \<div id=\"main2\" style=\"height:500px\;border:1px solid #ccc\;padding:10px\;\"\>\<\/div\>\n\n";
print OUT "    \<script src=\"$o_json\"\>\<\/script\>\n\n    \<script type=\"text\/javascript\"\>\n\n";
print OUT "    require\.config({\n        paths: {\n            echarts: \'$o_json2\'\n        }\n    })\;\n\n    require(\n";
print OUT "        [\n           \'echarts\'\,\n            \'echarts\/chart\/bar\'\,\n            \'echarts\/chart\/line\'\,\n";
print OUT "            \'echarts\/chart\/map\'\,\n            \'echarts\/chart\/scatter\'\n        ]\,\n        function (ec) {\n";
print OUT "            var myChart = ec\.init(document\.getElementById(\'main\'))\;\n            ";
print OUT "myChart\.setOption(option = {\n    title : {\n        text: \'Coverage Uniformity\, Sample: $s    (\>1:$more1\,$more1_per\% \; \<-1:$less_1\,$less_1_per\%)\'\,\n    }\,\n";
print OUT "    tooltip : {\n        trigger: \'axis\'\,\n        showDelay : 0\,\n        axisPointer:{\n            show: true\,\n";
print OUT "            type : \'cross\'\,\n            lineStyle: {\n                type : \'dashed\'\,\n                width : 1\n";
print OUT "            }\n        }\n    }\,\n    toolbox: {\n        show : true\,\n        feature : {\n            mark : {show: true}\,\n";
print OUT "            dataZoom : {show: true}\,\n            dataView : {show: true\, readOnly: false}\,\n            restore : {show: true}\,\n";
print OUT "            saveAsImage : {show: true}\n        }\n    }\,\n    xAxis : [\n        {\n            type : \'value\'\,\n";
print OUT "            scale:true\,\n            axisLabel : {\n                formatter: \'{value}\'\n            }\,\n            ";
print OUT "splitLine:{\n                show:false\n            }\n        }\n    ]\,\n    yAxis : [\n        {\n            type : \'value\'\,\n";
print OUT "            splitNumber:4\,\n            min: -2\,\n            max: 2\,\n            splitLine:{\n                show:false\n";
print OUT "            }\n        }\n    ]\,\n    series : [\n        {\n            type:\'scatter\'\,\n            data: [$data_normalized_coverage]\,\n";
print OUT "            markLine : {\n                data : [\n                    [{name: \'base line 1\'\,value: 1\, xAxis: 0\, yAxis: 1}\,{xAxis:10000\,yAxis: 1}]\,\n";
print OUT "                    [{name: \'base line -1\'\,value: -1\, xAxis: 0\, yAxis: -1}\,{xAxis:10000\,yAxis: \-1}]\n";
print OUT "                ]\n           }\n        }\n    ]\n})\;\n            var myChart2 = ec\.init(document\.getElementById(\'main2\'))\;\n";
print OUT "            myChart2\.setOption({\n            title : {\n            text: \'Region Coverage Distribution\, Sample: $s    ($less_strand_abs, $less_strand\% region lower than $strand_depth )\'\,\n";
print OUT "            }\,\n            tooltip : {\n               trigger: \'axis\'\n            },\n            toolbox: {\n                ";
print OUT "show : true\,\n            feature : {\n                 mark : {show: true}\,\n     dataZoom : {show: true}\,\n       dataView : {show: true\, readOnly: false}\,\n";
print OUT "                magicType : {show: true\, type: [\'line\'\,\'bar\']}\,\n              restore : {show: true}\,\n                ";
print OUT "                saveAsImage : {show: true}\n            }\n        }\,\n       calculable : true\,\n        xAxis : [\n       {\n";
print OUT "        type : \'category\'\,\n        data : [$data_coverge_name]\n        }\n        ]\,\n        yAxis : [\n        {\n        type : \'value\'\,\n";
print OUT "        splitArea : {show : true}\n       }\n        ]\,\n        series : [\n           {\n            type:\'bar\'\,\n            data:[$data_coverge]\,\n";
print OUT "        markLine : {\n            data : [\n            {type : \'average\'\, name: \'average\'}\n            ]\n        }\n        }\n";
print OUT "        ]\n        })\;\n    })\;\n\<\/script\>\n\<\/body\>\n\<\/html\>\n";
close OUT;



#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


