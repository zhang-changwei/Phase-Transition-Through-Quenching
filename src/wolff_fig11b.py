import numpy as np
from matplotlib import pyplot as plt
import MC

rng = np.random.default_rng()
L = 64
TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
TI = 1.5*TC
R = 1
ZW = 0.55
MCBIN = L**2
REPEAT = 100

cluster = MC.quench_wolff2(L, TC, tauq=2000, rf=False)
for i in range(REPEAT):
    c = MC.quench_wolff2(L, TC, tauq=2000, rf=False)
    cluster = np.vstack((cluster, c))
    print (i)
clusterave = np.average(cluster, axis=0)

np.save("fig11bL64Tauq2000.npy", clusterave)

plt.figure()
t = np.linspace(TI, TC, num=clusterave.size, endpoint=True)
plt.plot(t, clusterave)
plt.legend()
plt.yscale("log")
plt.show()