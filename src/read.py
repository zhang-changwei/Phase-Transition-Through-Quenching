import numpy as np
from matplotlib import pyplot as plt

TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
ZM = 2.17675
ZW = 0.55
TI = 1.5 * TC

rf = np.load("WFig11L12_rfave.npy")
rf2= np.load("WFig11L12_rfave.npy")

v = np.load("WFig11L12_v.npy")
v2= np.load("WFig11L24_v.npy")


# m2ave = np.load("fig1L8Tauq4000.npy")
# print (m2ave)
# plt.figure()
# plt.plot(np.arange(m2ave.size), m2ave)
# plt.show()

#v0= np.load("data\\MlinL8_v.npy")
# v = np.load("data\\MlinL12_v.npy")
# v2= np.load("data\\MlinL24_v.npy")
# v3= np.load("data\\MlinL48_v.npy")
# vn= np.load("data\\MnonlinL12_v.npy")
# L0=8
# L=12
# L2=24
# L3=48
# #m2ave0= np.load("data\\MlinL8_m2ave.npy")
# m2ave = np.load("data\\MlinL12_m2ave.npy")
# m2ave2= np.load("data\\MlinL24_m2ave.npy")
# m2ave3= np.load("data\\MlinL48_m2ave.npy")
# m2aven= np.load("data\\MnonlinL12_m2ave.npy")

plt.figure()
plt.plot(v*12**(ZW+1/NU), rf, label="L=12")
plt.plot(v2*24**(ZW+1/NU), rf2, label="L=12")


#plt.plot(t3, cn3, label="$\\tau_q=4000$")
#plt.plot(1/(v0), m2ave0*L0**2, label="L=8")
# plt.plot(1/(v), m2ave*L**2, label="L=12")
# plt.plot(1/(v3), m2ave3*L3**2, label="L=48")
# plt.plot(1/(v2), m2ave2*L2**2, label="L=24")
# plt.plot(1/(v), m2aven*L**2, label="L=12n")
# plt.plot(v*L**(ZM+1/NU), m2ave*L**(2*BETA*NU), label="L=12")
# plt.plot(v2*L2**(ZM+1/NU), m2ave2*L2**(2*BETA*NU), label="L=24")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()