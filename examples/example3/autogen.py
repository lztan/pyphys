import numpy as np
import solidstate as ss

lattdim =np.array([[8.8489999771,0,0],[0,8.8489999771,0],[0,0,12.6420001984]]) 
dims = [2,1,1]

mol = []
mol.append(ss.Atom("C",np.array([4.424499989,0.000000000,4.449984046])))
mol.append(ss.Atom("N",np.array([4.424499989,0.000000000,3.059364031])))
mol.append(ss.Atom("H",np.array([4.424499989,7.7889999771-8.8489999771,4.683952399])))
mol.append(ss.Atom("H",np.array([3.524499989,0.470000000,4.822151273])))
mol.append(ss.Atom("H",np.array([5.324499989,0.470000000,4.822151273])))
mol.append(ss.Atom("H",np.array([4.424499989,7.7889999771-8.8489999771,2.825395693])))
mol.append(ss.Atom("H",np.array([3.524499989,0.470000000,2.687196819])))
mol.append(ss.Atom("H",np.array([5.324499989,0.470000000,2.687196819])))

zpos = 0.5*(4.449984046+3.059364031)

#translate the molecule to the origin
mol0 = ss.trans(np.array([-4.424499989,0.0,-zpos]),mol)

cn = mol0[0:2]
ch3 = mol0[2:5]
nh3 = mol0[5:8]

#PbI3 lattice
pbi3 = []
pbi3.append(ss.Atom("Pb",np.array([0,0,0])))
pbi3.append(ss.Atom("Pb",np.array([0.000000000,0.000000000,6.321000099])))
pbi3.append(ss.Atom("Pb",np.array([4.424499989,4.424499989,6.321000099])))
pbi3.append(ss.Atom("Pb",np.array([4.424499989,4.424499989,0.000000000])))
pbi3.append(ss.Atom("I",np.array([0.000000000,0.000000000,3.125102415])))
pbi3.append(ss.Atom("I",np.array([0.000000000,0.000000000,9.446102703])))
pbi3.append(ss.Atom("I",np.array([4.424499989,4.424499989,9.446102703])))
pbi3.append(ss.Atom("I",np.array([4.424499989,4.424499989,3.125102415])))
pbi3.append(ss.Atom("I",np.array([1.895190272,6.319690129,0.058153202])))
pbi3.append(ss.Atom("I",np.array([6.953809837,2.529309848,0.058153202])))
pbi3.append(ss.Atom("I",np.array([2.529309848,1.895190272,0.058153202])))
pbi3.append(ss.Atom("I",np.array([6.319690129,6.953809837,0.058153202])))
pbi3.append(ss.Atom("I",np.array([1.895190272,2.529309848,6.379153154])))
pbi3.append(ss.Atom("I",np.array([6.953809837,6.319690129,6.379153154])))
pbi3.append(ss.Atom("I",np.array([2.529309848,6.953809837,6.379153154])))
pbi3.append(ss.Atom("I",np.array([6.319690129,1.895190272,6.379153154])))

# dims = [4,4,4]
# angles = [:,:,:,:]
# angles[x,y,z,j,0]= polar angle 
# angles[x,y,z,j,1]= azimuthal angle 
# angles[x,y,z,j,2]= H angle (C) 
# angles[x,y,z,j,3]= H angle (N)   
def makestruct(dims, angles):
    result = []
    zdir = np.array([0,0,1.0])
    xdir = np.array([1.0,0,0])
    for ix in range(dims[0]):
        for iy in range(dims[1]):
            for iz in range(dims[2]):
                #first rotate H
                ch31 = ss.rot(zdir,angles[ix,iy,iz,0,2],ch3)
                ch32 = ss.rot(zdir,angles[ix,iy,iz,1,2],ch3)
                ch33 = ss.rot(zdir,angles[ix,iy,iz,2,2],ch3)
                ch34 = ss.rot(zdir,angles[ix,iy,iz,3,2],ch3)

                nh31 = ss.rot(zdir,angles[ix,iy,iz,0,3],nh3)
                nh32 = ss.rot(zdir,angles[ix,iy,iz,1,3],nh3)
                nh33 = ss.rot(zdir,angles[ix,iy,iz,2,3],nh3)
                nh34 = ss.rot(zdir,angles[ix,iy,iz,3,3],nh3)

                mol1 = ch31+nh31+cn
                mol2 = ch32+nh32+cn
                mol3 = ch33+nh33+cn
                mol4 = ch34+nh34+cn
            
                #next, rotate CN axis
                mol1r = ss.rot(xdir,angles[ix,iy,iz,0,0],mol1)
                mol2r = ss.rot(xdir,angles[ix,iy,iz,1,0],mol2)
                mol3r = ss.rot(xdir,angles[ix,iy,iz,2,0],mol3)
                mol4r = ss.rot(xdir,angles[ix,iy,iz,3,0],mol4)

                mol1r = ss.rot(zdir,angles[ix,iy,iz,0,1],mol1r)
                mol2r = ss.rot(zdir,angles[ix,iy,iz,1,1],mol2r)
                mol3r = ss.rot(zdir,angles[ix,iy,iz,2,1],mol3r)
                mol4r = ss.rot(zdir,angles[ix,iy,iz,3,1],mol4r)

                # finally, translate CN to new position
                t1 = 0.5*lattdim[0,:]+0.0*lattdim[1,:]+0.25*lattdim[2,:] 
                t2 = 0.0*lattdim[0,:]+0.5*lattdim[1,:]+0.25*lattdim[2,:]
                t3 = 0.5*lattdim[0,:]+0.0*lattdim[1,:]+0.75*lattdim[2,:] 
                t4 = 0.0*lattdim[0,:]+0.5*lattdim[1,:]+0.75*lattdim[2,:]

                d = ix*lattdim[0,:] + iy*lattdim[1,:] + iz*lattdim[2,:]
                
                t1 = t1 + d
                t2 = t2 + d
                t3 = t3 + d
                t4 = t4 + d

                mol1t = ss.trans(t1,mol1r)
                mol2t = ss.trans(t2,mol2r)
                mol3t = ss.trans(t3,mol3r)
                mol4t = ss.trans(t4,mol4r)

                #make PbI3 lattice 

                pbi3t = ss.trans(d,pbi3)

                result = result + mol1t
                result = result + mol2t
                result = result + mol3t
                result = result + mol4t
                result = result + pbi3t

    return result


angles = np.zeros((dims[0],dims[1],dims[2],4,4))       

for ix in range(dims[0]):
    for iy in range(dims[1]):
        for iz in range(dims[2]):
            for j in range(4):
                angles[ix,iy,iz,j,0]= rand()*np.pi 
                angles[ix,iy,iz,j,1]=0
                angles[ix,iy,iz,j,2]=0
                angles[ix,iy,iz,j,3]=0            

structure =  makestruct(dims,angles) 
ss.printatoms("test",structure) 
