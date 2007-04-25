#!/usr/bin/perl
use strict;
open(MAN,'man.Rd') or die "no man.Rd file";

my $fh = undef;
while ( my $line = <MAN> ) {
  next if $line =~ m#^%#;
  if ( $line =~ m#__CUT__# ) {
    my ( $rd ) = $line =~ m#\{\s*(.+?)\s*\}#;
    $rd =~ s#^\.#0#;
    $rd = "man/$rd.Rd";
    close( $fh ) if defined( $fh );
    open( $fh, ">$rd" ) or die "couldn't open >$rd";
  }
  print $fh $line;
}
close( $fh ) if defined( $fh );
