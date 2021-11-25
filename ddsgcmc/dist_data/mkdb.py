#!/usr/bin/python

import numpy as np
import sys

Erange=0.0

print 'preprocessing...'
for a in xrange(1,int(sys.argv[1])+1):
    for b in xrange(1,int(sys.argv[1])+1):
        if a == b:
            continue
        rawAB = "sample%dto%d.dat" % (a,b)
        dEAB = []

        with open(rawAB,'r') as f:
            for line in f:
                data = line.strip().split()
                dEAB.append(float(data[0]))

        tmp = np.max(dEAB)-np.min(dEAB)
        print 'range = %f' % (tmp)
        if tmp > Erange:
            Erange = tmp
# Set database energy increment here
dE=0.001

for a in xrange(1,int(sys.argv[1])+1):
    for b in xrange(1,int(sys.argv[1])+1):
        if a == b:
            continue

#        print '%d to %d ...' % (a,b)
        rawAB = "sample%dto%d.dat" % (a,b)
        outAB = "dErsv_%dto%d.dat" % (a,b)

        print 'convert %s to %s ...' % (rawAB,outAB)

        dEAB = []

        with open(rawAB,'r') as f:
            for line in f:
                data = line.strip().split()
                dEAB.append(float(data[0]))

        eAB_min = np.min(dEAB)-dE
        eAB_max = eAB_min+Erange+dE
        nbinAB = int(np.ceil((eAB_max-eAB_min)/dE))
        hAB, binAB = np.histogram(dEAB,nbinAB,range=(eAB_min,eAB_max))

        eAB = []
        for i in xrange(1,len(binAB)):
            eAB.append((binAB[i-1]+binAB[i])/2.0)

        nAB = 1.0*np.array(hAB)
        normAB = np.sum(nAB)
        nAB /= float(normAB)
        cum=0.0
        with open(outAB,'w') as f:
            for i in xrange(len(eAB)):
                cum += nAB[i]
                f.writelines("%16.8e%16.8e%16.8e\n" % (eAB[i],nAB[i],cum))

