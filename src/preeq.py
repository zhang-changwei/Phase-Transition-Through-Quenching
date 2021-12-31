import numpy as np
from matplotlib import pyplot as plt
import MC

def Pre_Eq(lat, L, beta=1):
    mag = Mag(lat)
    E = 0
    log_mag2 = np.array([], dtype=np.float64) 
    log_E = np.array([], dtype=np.float64) 
    for i in range(L**2*150): 
        dmag, dE = Metropolis(lat, L, beta)
        mag += dmag
        E += dE
        if i%(L**2)==0: 
            log_mag2 = np.append(log_mag2, mag)
            log_E = np.append(log_E, E)
    return log_mag2, log_E

def Mag(lat):
    return lat.sum()

def MagDelta(cluster, updown):
    return -2*updown*len(cluster)

def EDelta(lat, L, x, y):
    return 2*lat[x][y]*(lat[(x-1)%L][y]+ lat[(x+1)%L][y]+ lat[x][(y-1)%L]+ lat[x][(y+1)%L])

def SpinFlip(lat, cluster, updown):
    for i in cluster:
        lat[i[0]][i[1]] = -updown

def Metropolis(lat, L, beta=1):
    posx,posy = rng.integers(0,L,endpoint=False), rng.integers(0,L,endpoint=False)
    dE = EDelta(lat, L, posx, posy)
    rand = rng.random()
    dmag = 0
    if dE<=0 or rand<np.exp(-beta*dE):
        dmag = -2*lat[posx][posy]
        lat[posx][posy] = -lat[posx][posy]
    else:
        dE = 0
    return dmag, dE

def main(L, ti):
    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))


    logmag2, logE = Pre_Eq(lat, L, 1/ti)

    return logmag2, logE

rng = np.random.default_rng()
m,E = main(32, 3.4)
print (E)

plt.figure()
plt.plot(np.arange(E.size), E, label="E")
plt.plot(np.arange(m.size), m, label="$m^2$")
plt.show()

# PRE equilibrium: 50 steps is enough