import numpy as np
from nonparam import KernalSmooth
from sklearn.linear_model import Lasso


class WeightFunction:
    beta = None
    linear_model = None

    
    def fit(self, X, U, y, h, alpha = 0, kernal = "gaussian"):
        p = np.shape(X)[1]
        n = np.shape(X)[0]
        ksm = KernalSmooth()
        g1 = ksm.fitpredict(U,y,U,0.3,"gaussian")
        g2 = ksm.fitpredict(U,X,U,0.3,"gaussian")
        
        self.linear_model = Lasso(alpha = alpha).fit(X-g2, y-g1)
        self.beta = self.linear_model.coef_ 
        
        