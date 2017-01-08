from heapq import heappush, heappop, heapify

def huffman_encode(symb2weights):
    """Huffman encode the given dict mapping symbols to weights
    {'a': 1673754, 'b': 7540662, 'c': 2818538}
    """

    heap = [[weight, [symbol, ""]] for symbol, weight in symb2weights.items()]
    heapify(heap)

    while len(heap) > 1:
        low = heappop(heap)
        high = heappop(heap)

        for pair in low[1:]:
            pair[1] = '0' + pair[1]

        for pair in high[1:]:
            pair[1] = '1' + pair[1]

        heappush(heap, [low[0] + high[0]] + low[1:] + high[1:])

    return sorted(heappop(heap)[1:], key=lambda p: (len(p[-1]), p))

file = [int(l) for l in open('huffman.txt')]
X = {}
i = 0
for element in file[1:]:
    X.update({str(i): element})
    i += 1

print X
huff = huffman_encode(X)
print "Symbol\tWeight\tHuffman Code"
for p in huff:
    print "%s\t\t%s\t%s" % (p[0], X[p[0]], p[1])
