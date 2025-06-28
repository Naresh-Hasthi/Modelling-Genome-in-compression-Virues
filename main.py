
from compressor import read_fasta, compress, decompress
from visualize import plot_compression
import os

input_file = "data/sars_cov2.fasta"
compressed_file = "output/encoded.bin"
decompressed_file = "output/decoded_sequence.txt"
stats_file = "output/stats.txt"

dna_sequence = read_fasta(input_file)
codes, bit_length = compress(dna_sequence, compressed_file)

original_size = len(dna_sequence.encode('utf-8'))
compressed_size = os.path.getsize(compressed_file)

decoded = decompress(compressed_file, codes, bit_length)
with open(decompressed_file, 'w') as f:
    f.write(decoded)

assert dna_sequence == decoded, "Mismatch after decompression!"

with open(stats_file, 'w') as f:
    f.write(f"Original Size: {original_size} bytes\n")
    f.write(f"Compressed Size: {compressed_size} bytes\n")
    f.write(f"Compression Ratio: {compressed_size/original_size:.2f}\n")

plot_compression(original_size, compressed_size)
