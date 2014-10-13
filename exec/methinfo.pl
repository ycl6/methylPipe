# This script extract methylation information from the bismark output sam file both single and paired end
# perl methinfo.pl --file /data/BA/kkishore/DNAmeth_analysis/RRBS/Methylpipe_scripts/LY31_methylpipe.sam --sam_type paired_sam --All /data/BA/kkishore/DNAmeth_analysis/RRBS/Methylpipe_scripts/All

#use warnings;
#use strict;
use Getopt::Long;
use List::Util qw(sum);

my $CpGfile ='';
my $Allfile= '';
my $file;
my $sam_type;
my $no_overlap;

GetOptions(
	   'file=s' =>\$file,
	   'sam_type=s' =>\$sam_type,
	   'no_overlap!' =>\$no_overlap,
	   'CpG=s'=> \$CpGfile,
	   'All=s'=> \$Allfile
);

unless($sam_type){die "sam_type not specified\n";}
unless($file){die "file not specified\n";}

my $IN;
if($file eq "-"){$IN=<STDIN>;}
elsif(! -e $file){
         print "file does not correspond to an existing file\n";
}
else{open ($IN,$file);}

my %sam_types = map { $_ => 1 }  qw(single_sam paired_sam);

if(! exists($sam_types{$sam_type})){
  print "--sam_type argument must be one of the following: 'single_sam','paired_sam' \n";
  
}

#open CpG output file and print the header
# print the methylations

if($sam_type eq "single_sam"){
  parse_sam($IN,$CpGfile,$Allfile,0,0);
}
elsif( $sam_type eq "paired_sam"){
  parse_sam($IN,$CpGfile,$Allfile,$no_overlap,1);
  
}


#########################  Functions ########################

# parse a given CG methlation 
# writes the filter passing CGs to output file

sub parseCGmeth
{
  my($CGmeth,$CGout)=(@_);
  foreach my $key (keys %{$CGmeth})
  {
        my($strand,$chr,$base)=split(/\|/,$key);
	my $noCs=$CGmeth->{$key}->[0];
        my $noTs=$CGmeth->{$key}->[1];
        if ($noCs != 0)
        {
            print $CGout join("\t",($chr,$base,$strand,"CG",$noCs,$noTs)),"\n"; 
        }
   }
  return 1;
}
  

sub parseAllmeth
{
  my($CGmeth,$CHHmeth,$CHGmeth,$Allout)=(@_);
  foreach my $key (keys %{$CGmeth})
  {
        my($strand,$chr,$base)=split(/\|/,$key);
	my $noCs=$CGmeth->{$key}->[0];
	my $noTs=$CGmeth->{$key}->[1];
        if ($noCs != 0)
        {
            print $Allout join("\t",($chr,$base,$strand,"CG",$noCs,$noTs)),"\n";
        }    
   }
  foreach my $key (keys %{$CHHmeth})
  {
        my($strand,$chr,$base)=split(/\|/,$key);
	my $noCs=$CHHmeth->{$key}->[0];
	my $noTs=$CHHmeth->{$key}->[1];
        if ($noCs != 0)
        {
            print $Allout join("\t",($chr,$base,$strand,"CHH",$noCs,$noTs) ),"\n";
        }    
  }
  foreach my $key (keys %{$CHGmeth})
  {
        my($strand,$chr,$base)=split(/\|/,$key);
	my $noCs=$CHGmeth->{$key}->[0];
	my $noTs=$CHGmeth->{$key}->[1];
        if ($noCs != 0)
        {
            print $Allout join("\t",($chr,$base,$strand,"CHG",$noCs,$noTs) ),"\n";
        }   
  }
 return 1;
}

# parse the methylation call string

sub parse_methcall{
	my($methcalls,$i,$key,$CGmeth,$CHHmeth,$CHGmeth)=@_;
	## CpG base
	if( uc($methcalls->[$i]) eq "Z")
		{
		  unless($CGmeth->{$key})
			{
			$CGmeth->{$key}=[0,0];
			} 
		  if( $methcalls->[$i] eq "Z" )
		    {
			$CGmeth->{$key}->[0]++;
		    }          
	  	  elsif( $methcalls->[$i] eq "z")
		    {
			$CGmeth->{$key}->[1]++;
		    }        
	 	}
	## Non CpG base
	else
		{                    
		  if(uc($methcalls->[$i]) eq "X" && (! $CHGmeth->{$key}) )
			{
	    		$CHGmeth->{$key}=[0,0];
	  		}
	  	 elsif(uc($methcalls->[$i]) eq "H" && (! $CHHmeth->{$key}))
			{
	    		$CHHmeth->{$key}=[0,0];
	  		}
		  if( $methcalls->[$i] eq "X" )
	        	{
		   	$CHGmeth->{$key}->[0]++;
	  		}
	  	elsif($methcalls->[$i] eq "H" )
	  		{
	    	     	$CHHmeth->{$key}->[0]++;
	  		}
	  	elsif( $methcalls->[$i] eq "x" )
			{
	    		$CHGmeth->{$key}->[1]++;
			}
		elsif( $methcalls->[$i] eq "h" )
			{
	    		$CHHmeth->{$key}->[1]++;
			}
		  }
    }


# parse the sam file
sub parse_sam{
  my($IN,$CpGfile,$Allfile,$no_overlap,$paired)=(@_);
  my $CGflag=0;
  my $Allflag=0;
  my $CGout;
  my $Allout;

  if($CpGfile ne ''){
    open ($CGout,">".$CpGfile);
    #print $CGout join("\t",qw(chr position strand context noCs noTs)),"\n";
    $CGflag=1;}
  if($Allfile ne '' ){
    open ($Allout,">".$Allfile);
    #print $Allout join("\t",qw(chr position strand context noCs noTs)),"\n";
    $Allflag=1;
  }

  # check if the file looks like sam
  #read-in file to count C bases
  my %CGmeth=(); 
  my %CHHmeth=();
  my %CHGmeth=();
  my $lastPos  =-1;
  my $lastChrom="null";
  my $chrPre;
  my $startPre=-1;

  while(<$IN>)
  {
    if($_=~/Bismark/){next;} # step over the header line
    if($_=~/^@/){next;} # step over the header line

    chomp;
    my @cols   = split(/\t/,$_);
    my $start  = $cols[3];
    my $end    = $start+length($cols[9])-1;
    my $chr    = $cols[2];
    my $methc  = $cols[13]; $methc =~ s/^XM:Z://;
    my @methcalls = split("", $methc); 
    my $rnext   = $cols[6];
    my $pnext   = $cols[7];
    my $tlen  = $cols[8];
    my $slen   = length($cols[9]);
    my $strand;
    if(	  $cols[14] eq "XR:Z:CT" && $cols[15] eq "XG:Z:CT"){$strand="+";} ## OT
    elsif($cols[14] eq "XR:Z:CT" && $cols[15] eq "XG:Z:GA"){$strand="-";} ## OB
    elsif($cols[14] eq "XR:Z:GA" && $cols[15] eq "XG:Z:CT"){$strand="+";} ## CTOT
    elsif($cols[14] eq "XR:Z:GA" && $cols[15] eq "XG:Z:GA"){$strand="-";} ## CTOB

    # if there is no_overlap trim the methcalls and $quals
    # adjust the start
    if($no_overlap && ( ($rnext eq "=") && $paired ) ){
      	if( ($start+$slen-1)>$pnext){
	if(($pnext-$start)>0){
	splice @methcalls,($pnext-$start);
	}
      }

    }
    
    #parses if start-LastPos>100
    if( ($start-$lastPos>100 && $lastPos != -1 ) || ($chr ne $lastChrom && $lastChrom ne "null"  ))
    {
      # if the user wants to write out files write them
      if($CGflag){ parseCGmeth(\%CGmeth,$CGout);}
      if($Allflag){ parseAllmeth(\%CGmeth,\%CHHmeth,\%CHGmeth,$Allout);}
      %CGmeth=();
      %CHHmeth=();
      %CHGmeth=();
    }

    # iterate over the mapped sequence
    for( my $i=0;$i< @methcalls; $i++) 
    {
	if ($methcalls[$i] eq "."){next;}
	my $key;
	if($strand eq "+"){$key=join("|",("+",$chr,$start+$i));}else{$key=join("|",("-",$chr,$start+$i));}
	parse_methcall(\@methcalls,$i,$key,\%CGmeth,\%CHHmeth,\%CHGmeth);
    }
    $lastPos=$end;
    $lastChrom=$chr;
  }
  close $IN;
  if($CGflag){ parseCGmeth(\%CGmeth,$CGout);}
  if($Allflag){ parseAllmeth(\%CGmeth,\%CHHmeth,\%CHGmeth,$Allout);}
  close $CGout;
  close $Allout;
}





