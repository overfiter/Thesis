import numpy as np
from scipy.linalg import solve
from numpy.linalg import cholesky    # Cholesky 分解：A = L * L^T

### 生成数据： y = beta * X + f(U) + error
### 其中 X 为多元正态分布，均值为0，相邻变量相关性为 rou ,间隔变量相关性乘指数衰减； 
### 其中 U 为均匀分布
### 其中 f() 为 正弦函数

def generater(n,p,k,rou):
    # generate real beta
    beta = np.append(5*np.ones(k) * (-1)**(np.random.randn(k)>0), np.zeros(p-k))
    np.random.shuffle(beta)

    # generate sigma
    index_p = np.array(range(p))+1
    sigma = np.abs(index_p.reshape((p,1)) - index_p)
    sigma = np.exp(sigma*np.log(rou))

    # generate X
    X = np.dot(np.random.randn(n, p), cholesky(sigma).T)
    X_norm = X / np.linalg.norm(X,axis=0)
    
    # generate U
    U = np.random.uniform(0, 2*np.pi, n)

    # generate error
    error = np.random.randn(n)

    # generate y
    y = np.dot(X,beta) + 10*np.sin(U) + error
    return y,X,U


def ortho(X):
    p = np.shape(X)[1]
    n = np.shape(X)[0]
    if np.shape(X)[0] < 2*np.shape(X)[1]:
        print("ERROR: n < 2p")
        return 0

    U = np.ones(np.shape(X))
    XU = np.hstack((X,U))

    for i in range(p):
        A = XU[n-p-i:,:p+i]
        b = -np.dot(XU[:n-p-i,:p+i].T,XU[:n-p-i,p+i])
        XU[n-p-i:,p+i] = solve(A.T, b) 
    return XU[:,p:]  

def augment(X1 ,X0 = None):
    if X0 is None:
        print("when augumenting X0 is None, so your input should be all X, or only X1 with appling spliding method")
    else:
        print("when augumenting X0 is not None, so \"data recycling\" method would be applied")
        
    p = np.shape(X1)[1]
    X = X1 if X0 is None else np.vstack((X0,X1))
    norm_X = np.linalg.norm(X,axis=0)
    X1_norm = X1 / norm_X
    X0_norm = None if X0 is None else X0 / norm_X
    U = ortho(X1)
    U_norm = U / np.linalg.norm(U,axis=0)
    
    # get X_tilde
    sigma = np.dot(X1_norm.T,X1_norm)
    sigma_inv = np.mat(sigma).I
    s = min(1,2*np.linalg.eig(sigma)[0].min() )* np.ones(p) - 0.0001   # 这里减去一个很小的数，否则可能导致下一行的分解由于非半正定而报错
    C = cholesky(2*np.diag(s) - np.diag(s)*sigma_inv*np.diag(s))
    X1_tilde = np.dot(X1_norm, np.eye(p) - sigma_inv*np.diag(s)) + np.dot(U_norm,C.T)
    
    XX = np.hstack((X1_norm,X1_tilde))
    
    if X0 is not None:
        XX0 = np.hstack((X0_norm,X0_norm))
        XX = np.vstack((XX0,XX))
    
    return np.array(XX) 