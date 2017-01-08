weights = [int(l) for l in open('mwis.txt')][1:]

def mwis(weights):

    n = len(weights)

    weights = [0] + weights

    maxsetweight = [0, weights[1]]

    for i in range(2, n + 1):
        maxsetweight.append(max(maxsetweight[i - 1], maxsetweight[i - 2] + weights[i] ))

    i = n
    maxset = []

    while i > 1:

        if maxsetweight[i-2] + weights[i] > maxsetweight[i-1]:
            maxset.append(i)
            i -= 2

            if i == 1:
                maxset.append(1)
                break
        else:
            i -= 1

    return (maxsetweight[n], maxset)

a, b = mwis(weights)
print "The weight of the maximum weight independent set of the graph is :", a
print "The vertices that constitute the maximum weight independent set of the path graph are :", b
