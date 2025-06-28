from collections import Counter
import heapq

# Build Huffman tree and return both encoding map and root
def build_huffman_tree(seq):
    if not seq:
        return None  # return None for empty sequence

    freq = Counter(seq)
    heap = [[wt, [ch, ""]] for ch, wt in freq.items()]
    heapq.heapify(heap)

    while len(heap) > 1:
        lo = heapq.heappop(heap)
        hi = heapq.heappop(heap)
        for pair in lo[1]:
            pair[1] = '0' + pair[1]
        for pair in hi[1]:
            pair[1] = '1' + pair[1]
        heapq.heappush(heap, [lo[0] + hi[0], lo[1] + hi[1]])

    encoding_map = dict(heap[0][1])
    return encoding_map

def huffman_compress(seq, tree):
    if not seq or not tree:
        return ""  # return empty string if no tree or sequence
    return ''.join(tree.get(ch, "") for ch in seq)

def huffman_decompress(encoded, tree):
    if not encoded or not tree:
        return ""
    
    reverse_tree = {v: k for k, v in tree.items()}
    current = ""
    decoded = ""

    for bit in encoded:
        current += bit
        if current in reverse_tree:
            decoded += reverse_tree[current]
            current = ""

    return decoded
