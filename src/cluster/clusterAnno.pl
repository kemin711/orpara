### perl ###
##!/usr/bin/perl -w

# annotate cluster according to the title dump
#@meaningless_wd = qw(a the);

# gprotein like ===> gprotein-like becomes one word
%misspelled= qw(dependant dependent associtaed associated responsiv responsive sesitive sensitive housekeeping house-keeping associted associated bidning binding I 1 II 2 III 3 IV 4);

@complexverbs = read_word_file("complex_verb");
@linkedwords = read_word_file("clusterAnno.conf");
#%lw_hash = wordfile_to_hash("clusterAnno.conf");

#$infile=$ARGV[0];
$infile= $ARGV[0] || "defallcluster";

open IN, "<$infile" or die "Cannot open $infile $!";
open OU, ">wdtoprt";
open CL, ">clwdfreq";

$_ = <IN>;
while (!/Cluster: /) { $_ = <IN>; }
while (1) {
	$cluster = substr($_, 9);
	print CL "\ncluster: $cluster";
	%wdfreq = ();
	%nbfreq = ();
	$_ = <IN>;
	while (!/==========/) {
		/(\w+) \| ([\w\._-]+): (.+)/;
		#$db = $1;
		$id = $2;
		$title = $3;
		#print $title, "\n";  # debug
		if (!$title) {
			die "no tile for $_\n";
		}
		addtohash(\%wdfreq, \%nbfreq, $title);
		#addtohash(\%wdfreq, $title);
		# for pick possible linked words, for initial 
		# setup of the program
		#if ($title =~ /-\s?(\w+)/ && !$lw_hash{$1} ) { 
		#	print $1, "\n"; 
		#}
		$_ = <IN>;
	}
	freqinfo(\%wdfreq, \%nbfreq);
	while ($_ && !/Cluster: /) { $_ = <IN>; }
	if (!$_) { last; }
}

######## sub routines ########################

# add the title's word to the hash table
sub addtohash {
	my $h_r = $_[0];
	my $hr_neighbor = $_[1];  # neighbor hash
	my $tl = preprocess($_[2]);
	my @arr = split /\s+/, $tl;
	my $i;
	my $last_w = "__start";

	for ($i=0; $i<@arr; $i++) {
		$w = $arr[$i];
		print OU "$w\t$id\n";
		$h_r->{$w}++;
		$hr_neighbor->{$last_w}{$w}++;
		$last_w = $w;
	}
	$hr_neighbor->{$last_w}{"__end"}++;
}

sub freqinfo {
	my $h_r = shift;
	my $hr_nbfreq = shift;
	foreach $key (keys %$h_r) {
		print CL $key, " ", $h_r->{$key}, "\n";
	}

	#print "$cluster\n";
	#foreach $key (keys %$hr_nbfreq) {
	#	print $key, "-> ";
	#	foreach $kk (keys %{$hr_nbfreq->{$key} } ) {
	#		print "$kk (",  $hr_nbfreq->{$key}{$kk}, ") ";
	#	}
	#	print "\n";
	#}
	#print "=============\n\n";

	print "\n$cluster\n";
	get_phrase($hr_nbfreq);
}

# preprocessing removes noise
# misspelling
# modifies the whole title, for linked words, add -
sub preprocess {
	my ($k, $v, $lw, $tl);
	$tl = $_[0];
	foreach $k (keys %misspelled) {  #spelling correction
		$v = $misspelled{$k};
		$tl =~ s/\b$k\b/$v/gi;
	}
	# add hyphen -
	my @cwarr;
	foreach $lw (@linkedwords) {
		if (@cwarr = get_complexverb($lw)) {
		   # search linked word
			while ($tl =~ /(\w+) $lw\b/gi) { 
				$with_tail = 0;
				foreach $cw (@cwarr) {
					if ($tl =~/\G$cw/) { 
						$with_tail = 1;
						last; 
					}
				}
				if (!$with_tail) {
					substr($tl, pos($tl)-length($lw)-1, 1) = '-';
				}
			}
		}
		else { # don't have to worry about complex word
			$tl =~ s/(\w+) $lw\b/$1-$lw/gi;
		}
	}
	$tl =~ s/[\(\)]/ /g;
	$tl =~ s/\// /g;
	$tl =~ s/\bG protein/G-protein/i;
	$tl =~ s/,//g;

	#if($_[0] =~ s/(\w+) derived/$1-derived/i) {
	#	print "modified ", $_[0], "\n";
	#}
	#if ($_[0] ne $tl) {
	#	print $_[0], "\n", $tl, " modified\n\n";
	#}
	return $tl;  # modified title
}

## returns a reference to the hash
# shoud sort these in input file
sub read_word_file {
	my $inf = $_[0];
	my @arr;
	open IN, "<$inf" or die "cannot open $inf $!\n";
	while ($_ = <IN>) {
		chomp;
		push(@arr, $_);
	}
	close IN;
	return @arr;
}

## if in the complex word list then return the complexed verb
# else return 0
sub get_complexverb {
	my $wd = $_[0];
	my @warr=();
	foreach $w (@complexverbs) {  # use global variable
		if ($w =~ /$wd (\w+)/) { 
			push(@warr, $1); 
		}
	}
	return @warr;  # at most two words, 99% are single words
}

sub wordfile_to_hash {
	my $inf = $_[0];
	my %h = ();
	open IN, "<$inf" or die "Cannot open $inf $!\n";
	while ($_=<IN>) {
		$h{$_}++;
	}
	close IN;
	return %h;
}

sub get_phrase {
	my $rh = shift;
	my %phrase_scores = ();  # hash of arrays
	my %phrase = ();
	foreach $k ( keys %{$rh->{"__start"}} ) {
		$r_arr_p = [];   # phrases
		$r_arr_ps = [];  # score
		$phrase{$k} = $r_arr_p;
		$phrase_scores{$k} = $r_arr_ps;

		push(@$r_arr_p, $k);
		push(@$r_arr_ps, $rh->{"__start"}{$k});

		#print "keys $k \n";

		$sumS=0;
		$cnt=0;
		while ($k ne "__end") {
			$maxS = 0;
			$maxKK="";
			foreach $kk ( keys %{$rh->{$k}} ) {
				if ($rh->{$k}{$kk} > $maxS) {
					$maxS = $rh->{$k}{$kk};
					$maxKK = $kk;
				}
			}
			$sumS += $maxS;
			$cnt++;
			#print "$k $maxS  ";
			print "$k ";
			$k = $maxKK;
			push(@$r_arr_p, $k);
			push(@$r_arr_ps, $maxS);
		}
		print $sumS/$cnt, "\n";
	}
}
		
