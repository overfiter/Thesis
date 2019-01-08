import numpy as np


class KernalSmooth:
    
    
    def kernal_function(self, x, kernal):
        if kernal == "gaussian":
            return 1/(np.sqrt(2*np.pi)) * np.exp(-(np.array(x)**2)/2)
        elif kernal == "uniform":
            return 0.5 * (np.abs(x) <= 1)
        elif kernal == "epanechnikov":
            return 0.75 * (1-x**2) * (np.abs(x) <= 1)
        elif kernal == "biquadrate":
            return (15/16) * (1-x**2)**2 * (np.abs(x) <= 1)
        else :
            print("your kernal function is not avaliable")
            
            
    def fitpredict(self, U, y, x, h, kernal = "gaussian"):
        y_hat = []
        for each_x in x:
            w = self.kernal_function((U - each_x)/h, kernal)/h
            w = w / w.sum()

            # 如果之预测一个变量
            if len(np.shape(y)) == 1:
                y_hat += [(w * y).sum()]

            # 如果需要预测多个变量
            elif len(np.shape(y)) == 2:
                multi_w = w
                for i in range( np.shape(y)[1]-1 ):
                    multi_w = np.vstack((multi_w,w))
                multi_w = multi_w.T
                y_hat.append((multi_w * y).sum(axis=0))

        return np.array(y_hat)