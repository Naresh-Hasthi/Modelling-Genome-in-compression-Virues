from Bio import SeqIO
from io import StringIO
from collections import Counter, namedtuple
import math
import bz2, gzip, lzma
import matplotlib.pyplot as plt
import heapq

# ==========================
# FASTA Reader
# ==========================
def read_fasta(file):
    content = file.read().decode("utf-8")
    handle = StringIO(content)
    record = next(SeqIO.parse(handle, "fasta"))
    return str(record.seq).upper()

# ==========================
# Entropy Calculation
# ==========================
def compute_entropy(seq, window=100):
    entropy_list = []
    for i in range(0, len(seq) - window, window):
        win = seq[i:i+window]
        probs = [win.count(nuc)/window for nuc in "ACGT"]
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        entropy_list.append(entropy)
    return entropy_list

# ==========================
# Genome Statistics
# ==========================
def compare_genomes(seq1, seq2):
    return {
        "Length 1": len(seq1),
        "Length 2": len(seq2),
        "GC% 1": round((seq1.count('G') + seq1.count('C')) / len(seq1) * 100, 2),
        "GC% 2": round((seq2.count('G') + seq2.count('C')) / len(seq2) * 100, 2)
    }

# ==========================
# Codon Usage Plot
# ==========================
def plot_codon_usage(seq1, seq2):
    def get_codons(seq):
        return [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

    c1 = Counter(get_codons(seq1))
    c2 = Counter(get_codons(seq2))
    codons = sorted(set(c1) | set(c2))
    v1 = [c1.get(c, 0) for c in codons]
    v2 = [c2.get(c, 0) for c in codons]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(codons, v1, label="Genome 1", color='blue')
    ax.plot(codons, v2, label="Genome 2", color='orange')
    ax.set_xticks(range(len(codons)))
    ax.set_xticklabels(codons, rotation=90)
    ax.legend()
    return fig

# ==========================
# Standard Compression
# ==========================
def gzip_compress(seq):
    return gzip.compress(seq.encode())

def bz2_compress(seq):
    return bz2.compress(seq.encode())

def lzma_compress(seq):
    return lzma.compress(seq.encode())

# ==========================
# Run-Length Encoding
# ==========================
def rle_compress(s):
    if not s:
        return ""
    out = []
    count = 1
    prev = s[0]
    for c in s[1:]:
        if c == prev:
            count += 1
        else:
            out.append(f"{prev}{count}")
            prev = c
            count = 1
    out.append(f"{prev}{count}")
    return ''.join(out)

# ==========================
# LZW Compression
# ==========================
def lzw_compress(uncompressed):
    dict_size = 256
    dictionary = {chr(i): i for i in range(dict_size)}
    w = ""
    result = []
    for c in uncompressed:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            result.append(dictionary[w])
            dictionary[wc] = dict_size
            dict_size += 1
            w = c
    if w:
        result.append(dictionary[w])
    return result

# ==========================
# Huffman Compression
# ==========================
class Node(namedtuple("Node", "char freq left right")):
    def __lt__(self, other):
        return self.freq < other.freq

def build_huffman_tree(freqs):
    heap = [Node(char, freq, None, None) for char, freq in freqs.items()]
    heapq.heapify(heap)
    while len(heap) > 1:
        node1 = heapq.heappop(heap)
        node2 = heapq.heappop(heap)
        merged = Node(None, node1.freq + node2.freq, node1, node2)
        heapq.heappush(heap, merged)
    return heap[0]

def generate_codes(node, prefix="", code_map=None):
    if code_map is None:
        code_map = {}
    if node.char is not None:
        code_map[node.char] = prefix
    else:
        generate_codes(node.left, prefix + "0", code_map)
        generate_codes(node.right, prefix + "1", code_map)
    return code_map

def huffman_compress(sequence):
    freqs = Counter(sequence)
    root = build_huffman_tree(freqs)
    huffman_codes = generate_codes(root)
    encoded_sequence = ''.join(huffman_codes[ch] for ch in sequence)
    return encoded_sequence, huffman_codes, len(sequence) * 8, len(encoded_sequence)

def huffman_decompress(encoded_sequence, code_map):
    reverse_map = {v: k for k, v in code_map.items()}
    decoded = ""
    code = ""
    for bit in encoded_sequence:
        code += bit
        if code in reverse_map:
            decoded += reverse_map[code]
            code = ""
    return decoded
