sub parse_match{
  my $match = pop;
  my @lines = split( /\n/, $match );
  
  foreach( @lines ){
    if( /Query: (\d+) (.*) (\d*)\n/gm ){
      $query_match .= $2;
    }
  }
  
  print "query match = $query_match\n";
}#end parse_match
