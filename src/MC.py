import numpy as np
from copy import deepcopy

# assuming J = 1, kB = 1

def Mag(lat):
    return lat.sum()

def MagDelta(cluster, updown):
    return -2*updown*len(cluster)

def EDelta(lat, L, x, y):
    return 2*lat[x][y]*(lat[(x-1)%L][y]+ lat[(x+1)%L][y]+ lat[x][(y-1)%L]+ lat[x][(y+1)%L])

def SpinFlip(lat, cluster, updown):
    for i in cluster:
        lat[i[0]][i[1]] = -updown

def Cluster(lat, L, seed, beta=1):
    # seed: tuple
    tx,ty = seed[0],seed[1]
    updown = lat[tx][ty]
    lat[tx][ty] = 0
    cluster = [] # cluster index
    cbound = [seed]
    cout = []
    while True:
        for b in cbound:
            x,y = b[0],b[1]
            # left
            rand = rng.random()
            if lat[(x-1)%L][y]==updown and rand<1-np.exp(-2*beta):
                lat[(x-1)%L][y] = 0
                cout.append(((x-1)%L,y))
            # right
            rand = rng.random()
            if lat[(x+1)%L][y]==updown and rand<1-np.exp(-2*beta):
                lat[(x+1)%L][y] = 0
                cout.append(((x+1)%L,y))
            # top
            rand = rng.random()
            if lat[x][(y-1)%L]==updown and rand<1-np.exp(-2*beta):
                lat[x][(y-1)%L] = 0
                cout.append((x,(y-1)%L))
            # bottom
            rand = rng.random()
            if lat[x][(y+1)%L]==updown and rand<1-np.exp(-2*beta):
                lat[x][(y+1)%L] = 0
                cout.append((x,(y+1)%L))
        if len(cout)==0:
            cluster = cluster + cbound
            break
        else:
            cluster = cluster + cbound
            cbound = cout
            cout = []

        if len(cluster)>=L**2:
            break
    return cluster, updown

def Wolff(lat, L, beta=1):
    posx,posy = rng.integers(0,L,endpoint=False), rng.integers(0,L,endpoint=False)
    seed = (posx,posy)
    cluster, updown = Cluster(lat, L, seed, beta)
    # SpinFlip
    SpinFlip(lat, cluster, updown)
    # MagDelta
    dmag = -2*updown*len(cluster)
    return dmag, len(cluster)

def Metropolis(lat, L, beta=1):
    posx,posy = rng.integers(0,L,endpoint=False), rng.integers(0,L,endpoint=False)
    dE = EDelta(lat, L, posx, posy)
    rand = rng.random()
    dmag = 0
    if dE<=0 or rand<np.exp(-beta*dE):
        dmag = -2*lat[posx][posy]
        lat[posx][posy] = -lat[posx][posy]
    return dmag

def bondSW(lat, L, beta):
    bondx = np.zeros((L,L), dtype=np.bool8)
    bondy = np.zeros((L,L), dtype=np.bool8)
    p = 1-np.exp(-2*beta) # bond possibility
    for i in range(L):
        for j in range(L):
            rand1,rand2 = rng.random(),rng.random()
            if lat[i][j]==lat[(i+1)%L][j] and rand1<p:
                bondx[i][j] = True
            else:
                bondx[i][j] = False
            if lat[i][j]==lat[i][(j+1)%L] and rand2<p:
                bondy[i][j] = True
            else:
                bondy[i][j] = False
    return bondx, bondy

def indexSW(L):
    index = []
    for i in range(L):
        for j in range(L):
            index.append((i,j))
    return index

def SwendsenWang(lat, L, beta=1):
    dmag = 0
    bondx, bondy = bondSW(lat, L, beta)
    index = indexSW(L)
    cluster = []
    clusterupdate = []
    clusterboundry = []
    while len(index)>0:
        cluster.append(index[0])
        ix, iy = index[0][0],index[0][1]
        index.remove((ix,iy))
        updown = lat[ix][iy]
        if bondx[ix][iy]==True and index.count(((ix+1)%L,iy))>0: 
            clusterupdate.append(((ix+1)%L,iy))
            index.remove(((ix+1)%L,iy))
        if bondx[ix-1][iy]==True and index.count(((ix-1)%L,iy))>0: 
            clusterupdate.append(((ix-1)%L,iy))
            index.remove(((ix-1)%L,iy))
        if bondy[ix][iy]==True and index.count((ix,(iy+1)%L))>0: 
            clusterupdate.append((ix,(iy+1)%L))
            index.remove((ix,(iy+1)%L))
        if bondy[ix][iy-1]==True and index.count((ix,(iy-1)%L))>0: 
            clusterupdate.append((ix,(iy-1)%L))
            index.remove((ix,(iy-1)%L))

        while len(clusterupdate)>0:
            clusterboundry = clusterupdate
            cluster = cluster + clusterupdate
            clusterupdate = []

            for i in clusterboundry:
                ix, iy = i[0],i[1]
                if bondx[ix][iy]==True and index.count(((ix+1)%L,iy))>0: 
                    clusterupdate.append(((ix+1)%L,iy))
                    index.remove(((ix+1)%L,iy))
                if bondx[ix-1][iy]==True and index.count(((ix-1)%L,iy))>0: 
                    clusterupdate.append(((ix-1)%L,iy))
                    index.remove(((ix-1)%L,iy))
                if bondy[ix][iy]==True and index.count((ix,(iy+1)%L))>0: 
                    clusterupdate.append((ix,(iy+1)%L))
                    index.remove((ix,(iy+1)%L))
                if bondy[ix][iy-1]==True and index.count((ix,(iy-1)%L))>0: 
                    clusterupdate.append((ix,(iy-1)%L))
                    index.remove((ix,(iy-1)%L))
        # online form one cluster
        # print (cluster)

        rand = rng.random()
        if rand < 1/2:
            dmag += MagDelta(cluster, updown)
            SpinFlip(lat, cluster, updown)
        cluster = []
        clusterupdate = []
        clusterboundry = []
    return dmag

def Pre_Eq(lat, L, beta=1):
    mag = Mag(lat)
    for i in range(L**2*50): 
        dmag = Metropolis(lat, L, beta)
        mag += dmag
    return mag

def quench_metropolis(L, tc, bin, tauq, r=1):
    ti = 1.5*tc
    v = 0.5/tauq**r

    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))
    mag = Pre_Eq(lat, L, 1/ti)

    log_mag2 = np.array([], dtype=np.float64) 

    mcstep, rcd = 0, 0
    t = (v*(tauq-mcstep)**r+1)*tc # temperature
    while mcstep <= tauq:
        beta = 1/t

        dmag = Metropolis(lat, L, beta)
        mag += dmag

        rcd += 1
        # record
        if rcd>=bin :
            log_mag2 = np.append(log_mag2, (mag/L**2)**2)
            rcd = 0
            mcstep = mcstep + 1
            # print (t)
            t = (v*(tauq-mcstep)**r+1)*tc
    
    return log_mag2

def quench_metropolis2(L, tc, bin, tauq, r=1):
    ti = 1.5*tc
    v = 0.5/tauq**r

    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))
    mag = Pre_Eq(lat, L, 1/ti)

    mcstep, rcd = 0, 0
    t = (v*(tauq-mcstep)**r+1)*tc # temperature
    while mcstep <= tauq:
        beta = 1/t

        dmag = Metropolis(lat, L, beta)
        mag += dmag

        rcd += 1
        # record
        if rcd>=bin :
            rcd = 0
            mcstep = mcstep + 1
            # print (t)
            t = (v*(tauq-mcstep)**r+1)*tc
    
    return (mag/L**2)**2

def quench_swendsenwang(L, tc, tauq):
    ti = 1.5*tc
    v = 0.5/tauq

    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))
    mag = Pre_Eq(lat, L, 1/ti)

    mcstep = 0
    while mcstep <= tauq:
        t = (v*(tauq-mcstep)+1)*tc # temperature
        beta = 1/t

        dmag = SwendsenWang(lat, L, beta)
        mag += dmag

        mcstep += 1
        # print (t)
    
    return (mag/L**2)**2

def quench_wolff(L, tc, tauq):
    ti = 1.5*tc
    v = 0.5/tauq

    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))
    mag = Pre_Eq(lat, L, 1/ti)

    mcstep = 0
    while mcstep <= tauq:
        t = (v*(tauq-mcstep)+1)*tc # temperature
        beta = 1/t

        dmag, clustersize = Wolff(lat, L, beta)
        mag += dmag

        mcstep += 1
        # print (t)
    
    return (mag/L**2)**2

def Rf(L, latinit, latend):
    return 1 - np.count_nonzero(latinit + latend)/L**2

def quench_wolff2(L, tc, tauq, rf=False):
    ti = 1.5*tc
    v = 0.5/tauq

    #lat = np.ones((L,L), dtype=np.int8)
    lat = rng.choice([1,-1], size=(L,L))
    mag = Pre_Eq(lat, L, 1/ti)

    log_c    = np.array([], dtype=np.float64)
    latinit = deepcopy(lat)

    mcstep = 0
    while mcstep <= tauq:
        t = (v*(tauq-mcstep)+1)*tc # temperature
        beta = 1/t

        dmag, clustersize = Wolff(lat, L, beta)
        mag += dmag

        # rcd
        log_c    = np.append(log_c, clustersize/L**2)

        mcstep += 1

    if rf==False:    
        return log_c
    else:
        latend = deepcopy(lat)
        rf = Rf(L, latinit, latend)
        return rf

rng = np.random.default_rng()
