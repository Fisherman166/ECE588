#!/usr/bin/perl

use warnings;
use strict;

my %info;
my %average_info;

&main();
exit(0);


sub main() {

    foreach my $filename (@ARGV) {
        next if($filename =~ /\.pl/);
        print "Opening file $filename\n";
        open(my $filehandle, '<', $filename) or die "Could not open file $filename for read\n";
        &read_file($filehandle);
        close($filehandle);
    }

    &calc_average_time();
    &write_data();
}

sub read_file() {
    my ($filehandle) = @_;

    while( my $line = <$filehandle> ) {
        my @data = extract_data($line);
        push(@{$info{$data[0]}->{'times'}}, $data[1]);
        push(@{$info{$data[0]}->{'speedups'}}, $data[2]);
    }
}

sub calc_average_time() {
    my $total_time = 0;
    my $time_count = 0;
    my $total_speedup = 0;
    my $speedup_count = 0;

    foreach my $thread_num (keys %info) {
        foreach my $runtime ( @{$info{$thread_num}->{'times'}} ) {
            $time_count++;
            $total_time += $runtime;
        }
        foreach my $speedup ( @{$info{$thread_num}->{'speedups'}} ) {
            $speedup_count++;
            $total_speedup += $speedup;
        }
        
        $average_info{$thread_num}->{'time'} = $total_time / $time_count;
        $average_info{$thread_num}->{'speedup'} = $total_speedup / $speedup_count;
        ($total_time, $time_count, $total_speedup, $speedup_count) = (0.0,0,0);
    }
}

sub write_data() {
    open(my $file, '>', "average.txt") or die "Could not open average.txt for write\n";

    foreach my $thread_num (sort { $a <=> $b } keys %average_info) {
        my $time = $average_info{$thread_num}->{'time'};
        my $speedup = $average_info{$thread_num}->{'speedup'};
        print $file "$thread_num\t$time\t$speedup\n";
    }

    close($file);
}

sub extract_data() {
    my $line = shift;
    chomp($line);
    my @values = split(/\t/, $line);
    return @values;
}

