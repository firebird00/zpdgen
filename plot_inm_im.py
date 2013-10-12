import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.colors import LogNorm
nlist=[[0,0],[1,0],[2,0],[0,2],[1,2],[2,2],[0,4],[1,4],[2,4]]
cnts=np.arange(-2,2.2,0.2)
wdts=np.ones(cnts.shape);
wdts[0]=2.0
wdts[5]=2.0
wdts[10]=4.0
wdts[15]=2.0
wdts[20]=2.0
vcb = np.arange(-12,14,2)
ii=0;
for ns in nlist:
    ii=ii+1;
    [n,m]=ns;
    filename='out_I'+str(n)+str(m)+'.h5'
    f = h5py.File(filename, "r")
    za=f['fields/za']
    eps=f['fields/eps']
    plt.subplot(3,3,ii)
#    plt.contourf(za['real'],za['imag'],eps['real'],500)
    plt.pcolormesh(za['real'],za['imag'],eps['imag'],shading='gouraud',rasterized=True)
#    ax=plt.gca()
#    ax.set_rasterization_zorder(-10)
    plt.clim(-10,10)
    cs=plt.contour(za['real'],za['imag'],eps['imag'],cnts,colors='k',linewidths=wdts*0.5)
#    plt.clabel(cs,[cnts[0],cnts[5],cnts[10],cnts[15],cnts[20]],inline=1,fmt="%1.1f",fontsize=10)
    aa=r'$Im[I_{'+str(n)+str(m)+r'}(\zeta_\alpha,\zeta_\beta,b)]$'
    plt.text(-1,3,aa,fontsize=14)
    if (m!=4 and n!=0):
        plt.tick_params(labelbottom='off',labelleft='off')
    elif(n==0 and m==4) :
        plt.text(2,-7.5,r'$Re[\zeta_\alpha]$',fontsize=14)
        plt.text(-7.5,3,r'$Im[\zeta_\alpha]$',fontsize=14,rotation=90)
        plt.tick_params(axis='both', which='major', labelsize=9,labelbottom='on',labelleft='on')
    elif(m==4) :
        plt.text(2,-7.5,r'$Re[\zeta_\alpha]$',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=9,labelbottom='on',labelleft='off')
    elif(n==0): 
        plt.text(-7.5,3,r'$Im[\zeta_\alpha]$',fontsize=14,rotation=90)
        plt.tick_params(axis='both', which='major', labelsize=9,labelleft='on',labelbottom='off')

plt.subplots_adjust(wspace=0.02,hspace=0.02,left=0.05,right=0.92,top=0.95,bottom=0.05)
fig=plt.gcf()
cax = fig.add_axes([0.93, 0.1, 0.02, 0.80])
plt.figure(2)
s=plt.contourf(za['real'],za['imag'],eps['imag']*(eps['imag']<14)*(eps['imag']>-14),100)
plt.clim(-10,10)
plt.close()
plt.figure(1)
cbar=fig.colorbar(s, cax=cax,ticks=vcb)
cbar.ax.tick_params(labelsize=8) 
cbar.solids.set_rasterized(True) 
#cbar.ax.set_ylim([cbar.norm(-10), cbar.norm(10)])
#cbar.outline.set_ydata([cbar.norm(-10)] * 2 + [cbar.norm(10)] * 4 + [cbar.norm(-10)] * 3)
fig.set_size_inches(8,8)
#plt.show()
plt.savefig("patates.pdf")
