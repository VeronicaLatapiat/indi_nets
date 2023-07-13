import numpy as np
from scipy.stats import norm
import pandas as pd

def csnet(data, c=None, alpha=0.01, boxsize=0.1, weighted=False, textdata=None):
    """
    Construct the cell-specific network (csn) from a gene expression matrix.
    
    Args:
    data: Gene expression matrix, rows = genes, columns = cells
    c: Construct the CSNs for all cells, set c = None (Default);
       Construct the CSN for cell k, set c = k
    alpha: Significance level (e.g., 0.001, 0.01, 0.05 ...)
           Larger alpha leads to more edges, Default = 0.01
    boxsize: Size of neighborhood, Default = 0.1
    weighted: True if edge is weighted, False (Default) if edge is not weighted
    textdata: Header names of the columns (genes) in the data matrix
    
    Returns:
    csn: Cell-specific network, the kth CSN is in csn[k]
         Rows = genes, columns = genes
    """
    
    if weighted is None:
        weighted = False
    
    if boxsize is None:
        boxsize = 0.1
    
    if alpha is None:
        alpha = 0.01
    
    n1, n2 = data.shape
    
    if c is None:
        c = np.arange(1, n2 + 1)
    
    upper = np.zeros((n1, n2))
    lower = np.zeros((n1, n2))
    
    for i in range(n1):
        s1 = np.sort(data[i, :])
        n3 = n2 - np.sum(np.sign(s1))
        h = round(boxsize / 2 * np.sum(np.sign(s1)))
        k = 0
        while k < n2:
            s = 0
            while k + s + 1 < n2 and s1[k + s + 1] == s1[k]:
                s += 1
            if s >= h:
                upper[i, k:k+s+1] = data[i, k]
                lower[i, k:k+s+1] = data[i, k]
            else:
                upper[i, k:k+s+1] = data[i, min(n2, k+s+h)]
                lower[i, k:k+s+1] = data[i, max(n3 * (n3 > h) + 1, k - h)]
            k = k + s + 1
    
    csn = {}
    B = np.zeros((n1, n2))
    p = norm.ppf(alpha, 0, 1)
    
    for k in c:
        for j in range(n2):
            B[:, j] = (data[:, j] <= upper[:, k-1]) & (data[:, j] >= lower[:, k-1])
        
        a = np.sum(B, axis=1)
        d = (np.dot(B, B.T) * n2 - np.outer(a, a)) / np.sqrt((np.outer(a, a)) * ((n2 - a) * (n2 - a)) / (n2 - 1) + np.finfo(float).eps)
        np.fill_diagonal(d, 0)
        
        if weighted:
            csn[k] = d * (d > 0)
        else:
            csn[k] = sparse.csr_matrix(d > p)
        
        print('Cell', k, 'is completed')
    
    if textdata is not None:
        dim = n1 - 1
        for l in range(len(csn)):
            for m in range(len(csn[l+1])):
                hij = csn[l+1][m, :].todense()
                A3 = textdata[m+1]
                A1 = m + 1
                str_filename = A3 + '_' + str(A1) + '.txt'
                pd.DataFrame(hij).to_csv(str_filename, sep='\t', header=False, index=False)
        
        fileID = open('Gene_names.txt', 'w')
        for gene_name in textdata:
            fileID.write(gene_name + ';\n')
        fileID.close()
    
    return csn
