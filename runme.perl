#!/usr/bin/perl -w

while (my $runno = shift @ARGV) {
    $cmd = sprintf("nohup /usr/bin/time ./agana.exe /daq/alpha_agmini/data/run%05dsub*.mid.lz4 >& run%05d.log &", $runno, $runno);
    print $cmd, "\n";
}

#end
