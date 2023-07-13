
import numpy as np
from scipy.stats import norm

def SSN(sample, ref):
    """
    Construct the SSN network structure.
    
    Args:
    sample: calculated sample
    ref: reference samples
    
    Returns:
    index_R: network structure
    p: p-values for significance testing
    
    Example usage:
    sample = new_T[:, 0]
    ref = new_N
    index_R, p = SSN(sample, ref)
    """
    
    R, P = np.corrcoef(ref.T)
    final_R0 = R.copy()
    final_R0[np.isnan(final_R0)] = 0
    
    NEW_data = np.hstack((ref, sample.reshape(-1, 1)))
    R1, P1 = np.corrcoef(NEW_data.T)
    final_R1 = R1.copy()
    final_R1[np.isnan(final_R1)] = 0
    
    index_R = final_R1 - final_R0
    
    m, n = ref.shape
    Z = index_R / ((1 - final_R0**2) / (n - 1))
    Z[np.isinf(Z)] = np.max(Z)
    Z[np.isneginf(Z)] = -np.max(Z)
    Z[np.isnan(Z)] = 0
    
    p = 1 - norm.cdf(np.abs(Z))
    
    return index_R, p
