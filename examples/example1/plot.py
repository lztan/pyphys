import matplotlib.pyplot as plt
import bandtools as bt
import units as un
import numpy as np

datx = bt.readmat("ramantot.dat")
#expt = bt.readmat("/home/lztan/science/feti/halide_perov/expt/cspbi3_300k.csv")

wmin = -200
wmax = 200
nw = 400
ws = np.linspace(wmin,wmax,nw)
thy = bt.broaden(datx[:,0]*1000/0.1241,datx[:,1],wmin,wmax,nw,2.0)

plt.plot(thy[:,0],thy[:,1]/10,"ro-",label="Theory (300K)")
#plt.plot(expt[:,0],expt[:,1]*2000,"k-",label="Expt. (300K)")

plt.legend(fontsize=14)


plt.xlim((-200,200))
plt.ylim((0,5000))

plt.xlabel(r"$\omega$ (cm$^{-1}$)")
plt.ylabel(r"Raman Spectrum")

plt.show()
