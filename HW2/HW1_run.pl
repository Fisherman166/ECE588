#!/usr/bin/perl

use strict;
use warnings;

&main();
exit(0);

###############################################################################
## Functions
###############################################################################
sub main() {
    my $max_count = 1000000000;
    my $max_threads = 16;
    my $output_filename = "output1.txt";
    my $single_thread_time;

    open(my $HW1_output_file, '>', $output_filename) or die "ERROR: Could not open $output_filename for write.";
    
    for(my $num_threads = 1; $num_threads <= $max_threads; $num_threads++) {
        my $cmd = "./sum $max_count $num_threads";
        my $output = `$cmd`;

        my $runtime_secs = &extract_seconds($output);
        my $speedup;
        if($num_threads > 1) {
            $speedup = &calc_speedup($single_thread_time, $runtime_secs);
        }
        else {
            $speedup = 1.0;
            $single_thread_time = $runtime_secs;
        }
        print $HW1_output_file "$num_threads\t$runtime_secs\t$speedup\n";
    }
    close($HW1_output_file);
}

sub extract_seconds() {
    my $output = shift;
    if( $output =~ /(\d+\.\d+)\s+sec/ ) {
        return $1;
    }
    return undef;
}

sub calc_speedup() {
    my ($single_thread_time, $speedup_time) = @_;
    my $speedup = $single_thread_time / $speedup_time;
    return $speedup;
}

