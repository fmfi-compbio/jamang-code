#! /usr/bin/perl -w

use strict;
use Data::Dumper;

my $USAGE = "$0 motifs probs fasta output_prefix

get a list of motifs with names (each line name, them motif)
get parameters (from file?):
  prob. of exiting bg state
  prob. of exiting motif submodel
  prob. of starting deletion
  prob. of extending deletion
  max deletion length
  prob of starting insertion
  prob of extending insertion
  prob of subtitution
get a fasta file
create an HMM
run HMM for each sequence
print bed file with stretches labeled by motifs or background
also for each sequence a summary of the number of times each match state used
as well as detailed viterbi path
    ";

# hmm: (states i,j 0..n, base b 0..3
# state 0 is background always starting
# $hmm->{'states'}[$i]{'emit'}[$b] = emmission prob
# $hmm->{'states'}[$i]{'in'}{$j} = transition prob
# $hmm->{'states'}[$i]{'symbol'} = A/C/G/T or - for background states
# $hmm->{'states'}[$i]{'l_emit'}[$b] = log emmission prob
# $hmm->{'states'}[$i]{'l_in'}{$j} = log transition prob
# $hmm->['states'}[$i]{'group'} = name of motif or 'bg'
# $hmm->['states'}[$i]{'type'} = 'match' or 'insert'
# $hmm->['states'}[$i]{'num'} = order within motif
# $hmm->['states'}[$i]{'key'} = "$group $type $num"
# $hmm->['states'}[$i]{'idx'} = $i
# $hmm->{'map'}[$key] =  state num

die $USAGE unless @ARGV==4;
my ($motifs_f, $probs_f, $fasta_f, $prefix) = @ARGV;
die $motifs_f unless -r $motifs_f;
die $probs_f unless -r $probs_f;
die $fasta_f unless -r $fasta_f;

my ($hmm, $motifs) = get_hmm($motifs_f, $probs_f);

my $in;
open $in, "<", $fasta_f or die;
my ($out1, $out2, $out3);
open $out1, ">", $prefix . ".bed" or die;
open $out2, ">", $prefix . ".stats"  or die;
open $out3, ">", $prefix . ".states"  or die;
while(1) {
    my ($name, $seq) = read_fasta($in);
    last unless defined $name;
    $name =~ s/^>//;
    my $states = viterbi($hmm, $$seq);

    my $n = length($$seq);
    die unless $n == scalar(@$states);

    # print the list of states and bed file,
    # accumulate counts for stats
    my $counts;
    my $prev_group = "";
    my $prev_group_start;
    my $first_bg = 0;
    my $last_bg = 0;
    for(my $i=0; $i<$n; $i++) {
	my $state = $hmm->{'states'}[$states->[$i]];
	my $base = substr($$seq, $i, 1);
	print $out3 join(" ", $name, $i, $state->{'key'}, $base, $state->{'symbol'}), "\n";
	if($state->{'type'} eq 'match') {
	    $counts->{$state->{'group'}}[$state->{'num'}]++;
	}
	my $group = $state->{'group'};
	if($prev_group ne $group) {
	    if(defined $prev_group_start) {
		print $out1 join("\t", $name, $prev_group_start,
				 $i, $prev_group), "\n";
		# if first group is bg, keep its length
		if($prev_group_start == 0 && $prev_group eq "bg") {  
		    $first_bg = $i-$prev_group_start;
		}
	    }
	    $prev_group = $group;
	    $prev_group_start = $i;
	}
    }
    print $out1 join("\t", $name, $prev_group_start,
		     $n, $prev_group), "\n";
    if($prev_group eq "bg") {
	$last_bg = $n-$prev_group_start;
	if($prev_group_start == 0) {
	    $first_bg = $last_bg;
	}
    }
    

    foreach my $motif (@$motifs) {
	print $out2 join(" ", $name, $motif->[0]);
	my $len = length($motif->[1]);
	for(my $i=0; $i<$len; $i++) {
	    my $count = $counts->{$motif->[0]}[$i];
	    $count = 0 unless defined $count;
	    print $out2 " $count";
	}
	print $out2 "\n";
    }
    print $out2 join(" ", $name, "bg", $counts->{"bg"}[0], $first_bg, $last_bg), "\n";
    
}
close $in;
close $out1;
close $out2;
close $out3;

exit(0);

#######
sub viterbi {
    my ($hmm, $seq) = @_;

    my $m = scalar @{$hmm->{'states'}};
    my $n = length($seq);
    my $prev_a = [(-inf) x $m];
    $prev_a->[0] = 0;
    my @b;
    
    for(my $pos=0; $pos<$n; $pos++) {
	my $base = base2num(substr($seq, $pos, 1));
	my @a;
	for(my $state=0; $state<$m; $state++) {
	    my $state_rec = $hmm->{'states'}[$state];
	    my $max = -inf;
	    my $max_state;
	    foreach my $prev (keys %{$state_rec->{'l_in'}}) {
		my $log = $prev_a->[$prev] + $state_rec->{'l_in'}{$prev};
		if($log > $max) {
		    $max = $log;
		    $max_state = $prev
		}
	    }
	    $max += $state_rec->{'l_emit'}[$base];
	    $a[$state] = $max;
	    $b[$pos][$state] = $max_state;
	}
	$prev_a = \@a;
    }

    my @states;
    my $last = 0;
    for(my $pos=$n-1; $pos>=0; $pos--) {
	push @states, $last;
	$last = $b[$pos][$last];
    }
    @states = reverse @states;
    return \@states;
}



#######
sub get_hmm {
    my ($motifs_f, $probs_f) = @_;

    my @motifs;
    my $in;
    open $in, "<", $motifs_f or die;
    my $num_match_states = 0;
    while(my $line = <$in>) {
	my @parts = split " ", $line;
	next if @parts==0;
	die "bad motif line '$line'" unless @parts==2;
	push @motifs, \@parts;
	$num_match_states += length($parts[1]);
    }
    close $in;
    
    my %probs;
    open $in, "<", $probs_f or die;
    while(my $line = <$in>) {
	my @parts = split " ", $line;
	next if @parts==0;
	die "bad prob line '$line'" unless @parts==2;
	die if exists $probs{$parts[0]};
	$probs{$parts[0]}=$parts[1];
    }
    close $in;
    die unless keys %probs==8;
    
    my $hmm = {'states'=>[], 'map'=>{}};
    add_state('bg', 'match', 0, '-', $hmm, \%probs);
    
    foreach my $rec (@motifs) {
	my ($name, $seq) = @$rec;
	for(my $i=0; $i<length($seq); $i++) {
	    my $c = substr($seq, $i, 1);
	    add_state($name, 'match', $i, $c, $hmm, \%probs);
	    add_state($name, 'insert', $i, '-', $hmm, \%probs);
	}
    }

    # transitions to/from bg
    foreach my $state (@{$hmm->{'states'}}) {
	if($state->{'idx'}>0 && $state->{'type'} eq 'match') {
	    add_tr(0, $state->{'idx'},
		   $probs{'exit_bg'}/$num_match_states, $hmm);
	    add_tr($state->{'idx'}, 0,
		   $probs{'exit_motif'}, $hmm);
	}
    }
    add_tr(0, 0, 1-$probs{'exit_bg'}, $hmm);

    # compute probabilities for deletions
    my @skip;
    my $p = $probs{'start_del'} * (1-$probs{'continue_del'});
    my $sum = 0;
    for(my $i=1; $i<=$probs{'max_del'}; $i++) {
	$skip[$i] = $p;
	$sum+=$p;
	$p *= $probs{'continue_del'};
    }
    # correct for only limited deltions considered
    for(my $i=1; $i<=$probs{'max_del'}; $i++) {
	$skip[$i] *= $probs{'start_del'} / $sum;
    }
    
    # transitions to next match, insert and deletions
    foreach my $rec (@motifs) {
	my ($name, $seq) = @$rec;
	my $len = length($seq);
	for(my $i=0; $i<$len; $i++) {
	    my $cur_match = get_state($hmm, $name, 'match', $i);
	    my $cur_insert = get_state($hmm, $name, 'insert', $i);
	    my $next_match = get_state($hmm, $name, 'match', ($i+1)%$len);
	    my $prob_next = 1-$probs{'exit_motif'}
	    -$probs{'start_del'}-$probs{'start_ins'};
	    add_tr($cur_match, $next_match, $prob_next, $hmm);
	    add_tr($cur_match, $cur_insert, $probs{'start_ins'}, $hmm);
	    add_tr($cur_insert, $cur_insert, $probs{'continue_ins'}, $hmm);
	    add_tr($cur_insert, $next_match, 1-$probs{'continue_ins'}, $hmm);

	    # deletions
	    for(my $j=1; $j<@skip; $j++) {
		my $index_after = ($i+$j+1) % $len;
		my $state_after = get_state($hmm, $name, 'match', $index_after);
		add_tr($cur_match, $state_after, $skip[$j], $hmm);		
	    }
	}
    }

    check_tr($hmm);
    
    return ($hmm, \@motifs);
}

#######
sub check_tr {
    my ($hmm) = @_;
    my $n = scalar @{$hmm->{'states'}};
    my @out = (0) x $n;
    my @out2 = "" x $n;
    for(my $i=0; $i<$n; $i++) {
	my $state = $hmm->{'states'}[$i];
	foreach my $j (keys %{$state->{'in'}}) {
	    my $p = $state->{'in'}{$j};
	    $out[$j]+=$p;
	    $out2[$j] .= "   " . $state->{'key'} . " " . $p;
	}
    }
    for(my $i=0; $i<$n; $i++) {
	die join("\t", "bad sum of outgoing for", $hmm->{'states'}[$i]{'key'},
		   $out[$i], $out2[$i])
	    unless abs($out[$i]-1)<1e-6;
    }    
}



#######
sub get_state {
    my ($hmm, $group, $type, $num) = @_;
    my $key = "$group $type $num";
    die unless exists $hmm->{'map'}{$key};
    return $hmm->{'map'}{$key};
}

    
#######
sub add_tr {
    my ($state1, $state2, $prob, $hmm) = @_;
    die "$prob $hmm->{'states'}[$state1]{'key'} $hmm->{'states'}[$state2]{'key'}"
	unless $prob>=0 && $prob<=1;

    $hmm->{'states'}[$state2]{'in'}{$state1} = $prob;
    $hmm->{'states'}[$state2]{'l_in'}{$state1} = log $prob;
}
    
	
#######
sub add_state {
    my ($group, $type, $num, $symbol, $hmm, $probs) = @_;

    my $idx = scalar(@{$hmm->{'states'}});
    my $key = "$group $type $num";
    my $state = {'group'=>$group,
            	 'type'=>$type,
		 'num'=>$num,
		     'key'=>$key,
		     'idx'=>$idx
    };
    if($symbol eq '-') {
	$state->{'emit'} = [0.25, 0.25, 0.25, 0.25];	
    } else {
	my $p = $probs->{'subst'} / 3;
	$state->{'emit'} = [$p, $p, $p, $p];
	my $i = base2num($symbol);
	$state->{'emit'}[$i] = 1-$probs->{'subst'}
    }
    $state->{'symbol'} = uc $symbol;
    $state->{'l_emit'} = [];
    for(my $b=0; $b<4; $b++) {
	$state->{'l_emit'}[$b] = log $state->{'emit'}[$b]
    }
    $hmm->{'states'}[$idx] = $state;
    $hmm->{'map'}{$key} = $idx;	
}


#######
sub base2num {
    my ($b) = @_;
    my $i = index("ACGT", uc $b);
    die "bad base $b" unless $i>=0;
    return $i;
}


############################
sub read_fasta {
    # Read one fasta sequence from the fasta file
    # Return undef, if no sequence found; 
    #        name and reference to sequence otherwise.
    
    my ($input) = @_;
    
    # read name of fasta sequence
    my $name;
    while(1) {
        my $line = <$input>;
        
        # end of file
        return unless defined $line;

        # skip empty lines
        next if $line =~ /^\s*$/;

        # parse the name
        $line =~ s/\s+$//;
        if($line =~ /^>/) {
            $name = $line; 
            last;
        }
        else { die "Incorrect fasta file '$line'."; }
    }

    # read fasta sequence
    my $sequence = "";
    while(1) {
        my $file_pos = tell($input);
        my $line = <$input>;
        
        # end of file
        last unless defined $line; 

        # if there is a beginning of a new sequence
        if($line =~ /^>/) {
            seek($input, $file_pos, 0);
            last;
        }

        # remove all whitespaces from line and append it to the sequence
        $line =~ s/\s//g;
        $sequence .= $line;
    }
    
    return ($name, \$sequence);
}
