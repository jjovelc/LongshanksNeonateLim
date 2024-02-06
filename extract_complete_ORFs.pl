#!/usr/bin/perl -w -s

my $fasta_file = $ARGV[0];

open my $fasta_fh, "<", $fasta_file or die "Could not open file '$fasta_file' $!";

my ($id, $seq) = ('', '');

while (my $line = <$fasta_fh>){
    chomp $line;
    if ($line =~ />.*$/) {  # Check if the line starts with '>'
        if ($line =~ />.*complete.*$/) {  # Check if 'complete' is in the line
            if (length($seq) > 0) {
                print "$id\n$seq\n";  # Print the previous record if it exists
                $seq = '';
                $id = '';  # Reset the ID
            }
            $id = $line;  # Set the new ID
        }
    } elsif ($id) {
        $seq .= $line;  # Concatenate sequence lines
    }
}

# Add the last sequence
print "$id\n$seq\n" if length($seq) > 0;

close $fasta_fh;

exit;
