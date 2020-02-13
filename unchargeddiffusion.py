# 2D gas diffusion using Metropolis-Hastings Algorithm

import random
import matplotlib.pyplot as plt
from math import e, pi


def ljpotential(sigma, epsilon, r):  # Lennard Jones potential
    A = 4*epsilon*(sigma**12)
    B = 4*epsilon*(sigma**6)

    V = (A/(r**12)) - (B/(r**6))
    
    return V

def distance(r1, r2):
    x1, y1 = r1
    x2, y2 = r2
    return ((x2 - x1)**2 + (y2 - y1)**2)**0.5

def boltzmannprob(E, kT):   
    return e**(-E/kT)

def averagev(T, M):
    v = ((8*8314*T)/(pi*M))**0.5

    return v

def display(points):
    n = len(points)
    x = []
    y = []

    for i in range(n):
        x.append(points[i][0])
        y.append(points[i][1])

    plt.scatter(x, y, s = 1.5)
    plt.show()

n = 500   # Iterations

boxx1 = 0    # Box size where particles are contained, units are 10000 pm
boxx2 = 30
boxy1 = 0
boxy2 = 30

plt.xlim(boxx1, boxx2)
plt.ylim(boxy1, boxy2)

startx1 = 14  # Starting region from where particles will diffuse
startx2 = 15
starty1 = 0
starty2 = 1

particles = 50  # Number of atoms simulated

particleinfo = [[] for i in range(particles)]   # Stores coordinates

for i in range(particles):
    xrandom = random.uniform(startx1, startx2)
    yrandom = random.uniform(starty1, starty2)
    while [xrandom, yrandom] in particleinfo:
        xrandom = random.uniform(startx1, startx2)
        yrandom = random.uniform(starty1, starty2)

    particleinfo[i] = [xrandom, yrandom]

display(particleinfo)

v = averagev(300, 4)      # m/s
timescale = 5            # ps
movementlimit = v*timescale/10000
sigma = 0.0256
epsilon = 0.0849
kT = 0.02587       # Assume T = 300K

for i in range(n):
    
    energytotal = 0
    
    for j in range(particles):
        particle = particleinfo[j]
        energyinitial = 0

        for k in range(particles):
            if (k == j):
                continue

            energyinitial += ljpotential(sigma, epsilon, distance(particle, particleinfo[k]))


        #print(energyinitial)

        xincr = random.uniform(-movementlimit, movementlimit)
        yincr = random.uniform(-movementlimit, movementlimit)
        xnew = particle[0] + xincr
        ynew = particle[1] + yincr
        
        while (((xnew > boxx2) or (xnew < boxx1)) or ((ynew > boxy2) or (ynew < boxy1))):
            xincr = random.uniform(-movementlimit, movementlimit)
            yincr = random.uniform(-movementlimit, movementlimit)
            xnew = particle[0] + xincr
            ynew = particle[1] + yincr

        while [xnew, ynew] in particleinfo:
            xincr = random.uniform(-movementlimit, movementlimit)
            yincr = random.uniform(-movementlimit, movementlimit)
            xnew = particle[0] + xincr
            ynew = particle[1] + yincr

            while (((xnew > boxx2) or (xnew < boxx1)) or ((ynew > boxy2) or (ynew < boxy1))):
                xincr = random.uniform(-movementlimit, movementlimit)
                yincr = random.uniform(-movementlimit, movementlimit)
                xnew = particle[0] + xincr
                ynew = particle[1] + yincr
            

        particlenew = [xnew, ynew]
        
        energyfinal = 0

        for k in range(particles):
            if (k == j):
                continue

            energyfinal += ljpotential(sigma, epsilon, distance(particlenew, particleinfo[k]))

        if (energyfinal <= energyinitial):
            particleinfo[j] = particlenew
            energytotal += energyfinal
            continue

        else:
            probability = boltzmannprob(energyfinal, kT)/boltzmannprob(energyinitial, kT)
            rand = random.uniform(0, 1)

            if (rand <= probability):
                particleinfo[j] = particlenew
                energytotal += energyfinal

            else:
                particleinfo[j] = particle
                energytotal += energyinitial

    if (((i+1)%100 == 0) or (i == 0)):
        print(energytotal*0.01037, "eV")
        plt.xlim(boxx1, boxx2)
        plt.ylim(boxy1, boxy2)
        display(particleinfo)

#plt.xlim(boxx1, boxx2)
#plt.ylim(boxy1, boxy2)
#display(particleinfo)

        

        


