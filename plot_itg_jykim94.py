import numpy as np
import matplotlib.pyplot as plt
import gpdf as gp
from scipy.optimize import root

etai=2.5
LnbyR=0.2
rbyR=0.18
kpar=0.1
tau=1.0

def epsfun(v):
    om=v[0]+1j*v[1]
    omsi=-ky
    omdi=2*omsi*LnbyR
    za=-om/omdi
    zb=-np.sqrt(2)*kpar/omdi
    b=ky**2
    i10=gp.Inm(za,zb,b,1,0)
    i12=gp.Inm(za,zb,b,1,2)
    i30=gp.Inm(za,zb,b,3,0)
    eps=1+1/tau+(i10*(om-omsi*(1-1.5*etai))-omsi*etai*(i12+i30))/omdi
    res=[float(np.real(eps)),float(np.imag(eps))]
    return res

kys=np.arange(0.01,2.4,0.01)
ky0=0.7
om0=[-0.2,0.2]
omky=np.zeros(np.shape(kys))*(1+1j)
ind0=np.argmin((kys-ky0)**2)
numk=len(kys)
inds=np.concatenate((np.arange(ind0,numk), np.arange(ind0-1,-1,-1)))
for l in inds:
    ky=kys[l]
    res=root(epsfun,om0,tol=1e-8,method='hybr')
    omky[l]=res.x[0]+1j*res.x[1]
    if(l==numk-1):
        om0=[float(np.real(omky[ind0])),float(np.imag(omky[ind0]))]
    else:
        om0=res.x
phi=np.loadtxt('jykim94a.dat')
ax1=plt.subplot(1,2,1)
plt.plot(kys,np.imag(omky),color='k',linewidth=2)
plt.plot(phi[:,0],phi[:,1],'--',color='r')
plt.xlabel('$k_y$',fontsize=18)
plt.ylabel('$\gamma$',fontsize=18)
plt.tight_layout()
#plt.axis([0.0,2.0,-0.1,0.2])
#plt.plot(kys,-np.real(omky),'r')
#plt.show()
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()
ax1.set_aspect((x1-x0)/(y1-y0))
ax1=plt.subplot(1,2,2)
phi=np.loadtxt('jykim94a_re.dat')
plt.plot(kys,np.real(omky),color='k',linewidth=2)
plt.plot(phi[:,0],-phi[:,1],'--',color='r')
plt.xlabel('$k_y$',fontsize=18)
plt.ylabel('$\omega_r$',fontsize=18)
plt.tight_layout()
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()
ax1.set_aspect((x1-x0)/(y1-y0))
plt.savefig("fig_kim94a.pdf")
gam=np.real(omky)
np.savetxt('test.txt',gam, delimiter=" ")
