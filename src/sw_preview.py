import numpy as np
from matplotlib import pyplot as plt
import MC

rng = np.random.default_rng()
L = 12
TC, BETA, NU = 2/(np.log(1+np.sqrt(2))), 1/8, 1
TI = 1.5*TC
R = 1
MCBIN = L**2
REPEAT = 100

#lat = np.ones((L,L), dtype=np.int8)
lat = rng.choice([1,-1], size=(L,L))
mag = MC.Mag(lat)
log_mag2 = np.array([], dtype=np.float64)

for i in range(100):
    dmag = MC.SwendsenWang(lat, L, 1/TI)
    mag += dmag 
    log_mag2 = np.append(log_mag2, (mag/L**2)**2)

plt.figure()
plt.plot(np.arange(log_mag2.size), log_mag2, label="L=12")
plt.legend()
plt.show()