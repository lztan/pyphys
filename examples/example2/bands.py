import qetools as qe
import bandtools as bt
import numpy as np

spa = qe.spaghetti("out_nscf")

nk2 = spa.bands.shape[1]
print(nk2)

nk = 4096

nv = 392
nb = 10

#bt.printmat("kpts", spa.kpts)

upbands = spa.bands[nv-nb:nv+nb,0:nk]
dnbands = spa.bands[nv-nb:nv+nb,nk:]

bt.printmat("up.dat",np.transpose(upbands))
bt.printmat("dn.dat",np.transpose(dnbands))
