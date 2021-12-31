import numpy as np
from matplotlib import pyplot as plt
import MC

rng = np.random.default_rng()
L = 8
TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
TI = 1.5*TC
MCBIN = L**2
REPEAT = 100

m20 = MC.quench_metropolis(L, TC, MCBIN, 4000)
for i in range(REPEAT):
    m2 = MC.quench_metropolis(L, TC, MCBIN, 4000)
    m20 = np.vstack((m20,m2))
    print (i)
m2ave, m2std = np.average(m20, axis=0), np.std(m20, axis=0)

np.save("fig1L8Tauq4000.npy", m2ave)


plt.figure()
plt.plot(np.arange(m2ave.size), m2ave)
plt.show()