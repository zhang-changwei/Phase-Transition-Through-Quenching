import numpy as np
from matplotlib import pyplot as plt
import MC

rng = np.random.default_rng()
L = 12
TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
TI = 1.5*TC
R = 2
MCBIN = L**2
REPEAT = 100

# tauq = np.logspace(0, 4, num=41, base=10, endpoint=True)
tauq = np.logspace(0, 3, num=31, base=10, endpoint=True)
v = 0.5/tauq**R
m2ave = np.array([], dtype=np.float64)
m2std = np.array([], dtype=np.float64)

for j in range(tauq.size):
    m20 = np.array([], dtype=np.float64)
    for i in range(REPEAT):
        m2 = MC.quench_metropolis2(L, TC, MCBIN, tauq[j], R)
        m20 = np.append(m20, m2)
        print (j, i, m2)
    m2ave = np.append(m2ave, np.average(m20))
    m2std = np.append(m2std, np.std(m20)) 

np.save("data\\MnonlinL12_m2ave.npy", m2ave)
np.save("data\\MnonlinL12_m2std.npy", m2std)
np.save("data\\MnonlinL12_tauq.npy", tauq)
np.save("data\\MnonlinL12_v.npy", v)


plt.figure()
plt.plot(1/v, m2ave*L**2, label="L=12")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()