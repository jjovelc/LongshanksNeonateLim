#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# Check for command-line arguments
my ($inpep, $ids_list) = @ARGV;
die "Usage: $0 <pep_file> <ids_file>\n" unless $inpep && $ids_list;

# Read IDs into a hash
my %ids;
open my $ID_FILE, '<', $ids_list or die "Can't open $ids_list: $!";
while (my $long_id = <$ID_FILE>) {
    chomp $long_id;
    my ($id, @temp) = split(' ', $long_id);
    $id =~ s/>//;
    $ids{$id} = join('', @temp);
}
close $ID_FILE;

# Process the FASTA file and filter sequences
my $seq_in = Bio::SeqIO->new(
    -format   => 'fasta',
    -file     => $inpep,
    -alphabet => 'protein'
);

while (my $seq = $seq_in->next_seq) {
    my $target_id = $seq->id;
    my $desc = $seq->desc || '';
    if (exists($ids{$target_id})) {
      print ">", $target_id, " ", $ids{$target_id}, " ", "\n", $seq->seq, "\n";
    }
}

exit;
