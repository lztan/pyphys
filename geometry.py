import numpy as np

def unitfy(v):
  """normalize a vector"""
  norm2 = np.vdot(v,v)
  return v / np.sqrt(norm2)

def rotmat(u,th):
  """rotation matrix for angle th, axis u"""
  ans = np.zeros((3,3))
  nu = unitfy(u)
  ux = nu[0]
  uy = nu[1]
  uz = nu[2]
  costh = np.cos(th)
  sinth = np.sin(th)
  ans[0,0] = costh + ux**2 * (1-costh)
  ans[0,1] = ux*uy*(1-costh)-uz*sinth
  ans[0,2] = ux*uz*(1-costh)+uy*sinth
  ans[1,0] = uy*ux*(1-costh)+uz*sinth
  ans[1,1] = costh + uy**2 * (1-costh)
  ans[1,2] = uy*uz*(1-costh)-ux*sinth
  ans[2,0] = uz*ux*(1-costh)-uy*sinth
  ans[2,1] = uz*uy*(1-costh)+ux*sinth
  ans[2,2] = costh + uz**2 * (1-costh)
  return ans


