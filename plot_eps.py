import numpy as np
import matplotlib.pyplot as plt
import gpdf as gp

etai=2.5
LnbyR=0.2
rbyR=0.18
kpar=0.1
tau=1.0
ky=0.06

def epsfun(v):
    om=v
    omsi=-ky
    omdi=2*omsi*LnbyR
    za=-om/omdi
    zb=-kpar/omdi*np.sqrt(2)
    b=ky**2
    i10=gp.Inm(za,zb,b,1,0)
    i12=gp.Inm(za,zb,b,1,2)
    i30=gp.Inm(za,zb,b,3,0)
    eps=1+1/tau+(i10*(om-omsi*(1-1.5*etai))-omsi*etai*(i12+i30))/omdi
    return eps
xx,yy=np.meshgrid(np.arange(-1,1,0.01),np.arange(-1,1,0.01))
om=xx+1j*yy;
cnts=np.arange(-2,2.2,0.2)
wdts=np.ones(cnts.shape);
wdts[0]=2.0
wdts[5]=2.0
wdts[10]=4.0
wdts[15]=2.0
wdts[20]=2.0
eps=epsfun(om)
#gp.Inm(om,0.0,0.09,1,0)
plt.pcolormesh(np.real(om),np.imag(om),np.imag(eps),shading='gouraud',rasterized=True)
plt.clim(-0.01,0.01)
plt.colorbar()
plt.contour(np.real(om),np.imag(om),np.real(eps),[0.0],colors='k')
plt.contour(np.real(om),np.imag(om),np.imag(eps),[0.0],colors='r')
plt.show()
