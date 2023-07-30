from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

import numpy as np

A1 = np.array([0] * 100)
A2 = np.array([0.868812657252926] * 100)
A3 = np.array([0.025] * 100)

A2[0] = 1
A2[-1] = 1

A3[0] = 0
A3[-1] = 0

Bij = np.array([0.474509609679611] * 100)
Bij[0] = 0
Bij[-7] = 0.474501171553038
Bij[-6] = 0.474526485932756
Bij[-5] = 0.474222713376145
Bij[-4] = 0.477522020866004
Bij[-3] = 0.448823952392832
Bij[-2] = 0.630150854309913
Bij[-1] = 10


A = lil_matrix((100, 100))
A.setdiag(A3, k=1)
A.setdiag(A2, k=0)
A.setdiag(A1, k=-1)

B = np.asarray(Bij)
A = A.tocsr()

C = spsolve(A, B).transpose()
