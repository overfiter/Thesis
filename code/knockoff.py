import numpy as np
from scipy.linalg import solve
from numpy.linalg import cholesky    # Cholesky 分解：A = L * L^T

from sklearn.linear_model import Lasso


class LinearInferenceMachine:
    precision = None
    XX = None
    y = None
    record = None
    alpha = None
    feature_index = None
    W = None
    p = None
    
    def __init__(self,XX,y,precision):
        self.precision = precision
        self.XX = XX
        self.y = y
        self.record = set()
        self.alpha = []
        self.feature_index = []
        self.p = int(np.shape(XX)[1]/2)
        self.W = -np.ones(self.p)

    def dialectical(self, alpha_low, alpha_high):
        model_high = Lasso(alpha=alpha_high,fit_intercept=0).fit(self.XX,self.y)
        model_low = Lasso(alpha=alpha_low,fit_intercept=0).fit(self.XX,self.y)
        new_element = set(model_low.coef_.nonzero()[0]) - set(model_high.coef_.nonzero()[0]) - self.record

        if len(new_element) == 0:
            pass

        elif len(new_element) == 1:

            local_alpha_low = alpha_low
            local_alpha_high = alpha_high
            while local_alpha_high - local_alpha_low > self.precision:
                alpha_mid = (local_alpha_high + local_alpha_low)/2
                model_mid = Lasso(alpha=alpha_mid,fit_intercept=0).fit(self.XX,self.y)
                model_low = Lasso(alpha=local_alpha_low,fit_intercept=0).fit(self.XX,self.y)
                find_mid = set(model_mid.coef_.nonzero()[0])
                find_low = set(model_low.coef_.nonzero()[0])
                if len(find_low - find_mid - self.record) == 0:
                    local_alpha_low = alpha_mid
                else:
                    local_alpha_high = alpha_mid

            self.record.add(list(new_element)[0])
            self.feature_index.append(list(new_element)[0])
            self.alpha.append(local_alpha_low)

        elif len(new_element) > 1:
            alpha_mid = (alpha_low + alpha_high)/2
            self.dialectical(alpha_mid,alpha_high)
            self.dialectical(alpha_low,alpha_mid)
    
    
    # 计算统计量W
    def compute_W(self):
        # 确认下界
        alpha_low = 0  # 下界
        if len(Lasso(alpha=alpha_low,fit_intercept=0).fit(self.XX,self.y).coef_.nonzero()[0]) != 2 * self.p:
            return "ERROR: some coef is equal to zero strictly"

        # 找到上界
        alpha_high = alpha_low
        while len(Lasso(alpha=alpha_high,fit_intercept=0).fit(self.XX,self.y).coef_.nonzero()[0]) != 0:
            alpha_high += 1

        # 递归二分法, 招到每个特征第一次进入lasso时的alpha
        # 结果保存到全局变量alpha、feature_index中
        self.dialectical(alpha_low,alpha_high)

        # 根据alpha、feature_index 计算每个特征的 Wj
        for j in range(self.p):
            wj = self.alpha[self.feature_index.index(j)]
            wj_tilde = self.alpha[self.feature_index.index(j+self.p)]
            if wj > wj_tilde:
                self.W[j] = wj
            else:
                self.W[j] = -wj_tilde
        