#!/usr/bin/perl
use Bio::SeqIO;

my $sequences = Bio::SeqIO->new( 
    -file   => "/dev/stdin",
        -format => "fasta",
        );
        
        while ( my $dna = $sequences->next_seq ){
            my $protein = $dna->translate( 
                    -codontable_id => 1, # standard genetic code
                            -frame         => 0, #reading-frame offset 0
                                );
                                    print ">", $dna->display_id, "\n";
                                        print $protein->seq, "\n\n";
                                        }
