import numpy
import pandas
from scipy.stats import norm
from scipy import sparse

eps = numpy.finfo(float).eps
# eps

# archivo a leer y opciones
input = 'test2.txt'
with open(input, 'r') as infile:
    data = pandas.read_csv(infile, delimiter = ' ', index_col = 0).values

cell = range(0, numpy.shape(data[:,:])[1], 1) # all cells; therefore, set cell to a list of cell indexes
alpha = 0.01
boxsize = 0.1
weighted = True

n1, n2 = numpy.shape(data)
upper = numpy.zeros((n1, n2))
lower = numpy.zeros((n1, n2))

for i in range(0, n1):
    s1 = numpy.sort(data[i, :])
    s2 = numpy.argsort(data[i, :])
    n0 = int(n2 - sum(numpy.sign(s1)))
    h = int(numpy.around(boxsize / 2 * sum(numpy.sign(s1))))
    k = 1

    while k <= n2:
        s = 0

        while k + s < n2 and s1[k + s] == s1[k - 1]: # match index with -1
            s = s + 1

        if s >= h:
            upper[i, s2[k-1:k+s]] = data[i, s2[k - 1]]
            lower[i, s2[k-1:k+s]] = data[i, s2[k - 1]]
        else:
            upper[i, s2[k-1:k+s]] = data[i, s2[min(n2, k+s+h) - 1]]
            lower[i, s2[k-1:k+s]] = data[i, s2[max(n0 * int(n0 > h), k - h) - 1]]
        k = k + s + 1

# Construction of cell-specific network
csn = [ numpy.empty((1,1)) for x in range(0, n2) ]
B = numpy.zeros((n1, n2))
p = -norm.ppf(alpha,0,1)

for k in cell:
    for j in range(0, n2):
        B[:, j] = (data[:, j] <= upper[:, k]).astype(int) * (data[:, j] >= lower[:, k]).astype(int)

    a = numpy.sum(B, 1)
    d = (B.dot(B.T) * n2 - numpy.outer(a, a))             / (numpy.sqrt(numpy.outer(a, a) * numpy.outer(n2-a, n2-a) / (n2 - 1) + eps))

    numpy.fill_diagonal(d, 0) # return None as stdout

    if weighted:
        csn[k] = d * (d > 0).astype(int) # carefull, these networks are different from the sparse ones
    else:
        csn[k] = sparse.csr_matrix((d > p).astype(int)).T # carefull with the interpretation of these networks

# write reports
if weighted:
    for index, network in enumerate(csn):
        with open('network_cell{:0d}'.format(index), 'w') as outfile:
            pandas.DataFrame(network).to_csv(outfile, sep = ' ', header = False, index = False)
else:
    for index, network in enumerate(csn):
        print(network.toarray())
        with open('network_cell{:0d}'.format(index), 'w') as outfile:
            pandas.DataFrame(network.toarray()).to_csv(outfile, sep = ' ', header = False, index = False)
