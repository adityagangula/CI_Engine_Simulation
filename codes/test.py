import numpy as np
import matplotlib.pyplot as plt

z = np.array([ 9.27185689e-30, -2.69196802e-25,  3.33277639e-21, -2.29855128e-17, 9.68752486e-14, -2.57578506e-10,  4.28980804e-07, -4.25098468e-04, 2.31965986e-01 , 7.85901406e-03])
def f(N):
    x = 0
    for i in range(0,10):
        x = x + z[i]*(N**(9-i))
    return x

x = np.array([0,1,2,3,4,5,6,7,8,9])
index = np.where(x == np.amax(x))
print(index[0][0])

