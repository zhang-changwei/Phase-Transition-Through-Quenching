import numpy as np
from matplotlib import pyplot as plt
import MC

rng = np.random.default_rng()
L = 48
TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
TI = 1.5*TC
R = 1
ZW = 0.55
MCBIN = L**2
REPEAT = 100

tauq = np.logspace(0, 3, num=31, base=10, endpoint=True)
v = 0.5/tauq**R
rfave = np.array([], dtype=np.float64)
rfstd = np.array([], dtype=np.float64)

for j in range(tauq.size):
    rf0 = np.array([], dtype=np.float64)
    for i in range(REPEAT):
        rf = MC.quench_wolff2(L, TC, tauq[j], rf=True)
        rf0 = np.append(rf0, rf)
        print (j, i, rf)
    rfave = np.append(rfave, np.average(rf0))
    rfstd = np.append(rfstd, np.std(rf0)) 

np.save("WFig11L48_rfave.npy", rfave)
np.save("WFig11L48_rfstd.npy", rfstd)
np.save("WFig11L48_tauq.npy", tauq)
np.save("WFig11L48_v.npy", v)


plt.figure()
plt.plot(v*L**(ZW+1/NU), rfave, label="L=12")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()