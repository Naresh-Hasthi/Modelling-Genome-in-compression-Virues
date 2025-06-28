
import heapq
from collections import Counter

class Node:
    def __init__(self, symbol, freq):
        self.symbol = symbol
        self.freq = freq
        self.left = None
        self.right = None

    def __lt__(self, other):
        return self.freq < other.freq

def build_huffman_tree(data):
    freq = Counter(data)
    heap = [Node(symbol, freq) for symbol, freq in freq.items()]
    heapq.heapify(heap)
    while len(heap) > 1:
        n1 = heapq.heappop(heap)
        n2 = heapq.heappop(heap)
        merged = Node(None, n1.freq + n2.freq)
        merged.left = n1
        merged.right = n2
        heapq.heappush(heap, merged)
    return heap[0]

def build_codes(root):
    codes = {}
    def _build_codes(node, code=""):
        if node:
            if node.symbol is not None:
                codes[node.symbol] = code
            _build_codes(node.left, code + "0")
            _build_codes(node.right, code + "1")
    _build_codes(root)
    return codes
