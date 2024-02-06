#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# Check for command-line arguments
my ($infile, $ids_list) = @ARGV;
die "Usage: $0 <fastA_file> <ids_file>\n" unless $infile && $ids_list;

# Read IDs into a hash
my %ids;
open my $ID_FILE, '<', $ids_list or die "Can't open $ids_list: $!";
while (my $long_id = <$ID_FILE>) {
    chomp $long_id;
    my @temp = split(' ', $long_id);
    my $id = $temp[0];
    $id =~ s/>//;
    $ids{$id} = 1;
}
close $ID_FILE;

# Process the FASTA file and filter sequences
my $seq_in = Bio::SeqIO->new(
    -format   => 'fasta',
    -file     => $infile,
    -alphabet => 'dna'
);

while (my $seq = $seq_in->next_seq) {
    my $target_id = $seq->id;
    if (exists $ids{$target_id}) {
            print ">", $target_id, " ", "\n", $seq->seq, "\n";
    }
}
