import numpy as np
import matplotlib.pyplot as plt
import gpdf as gp
nlist=[[1,0],[2,0],[3,0],[1,2],[2,2],[3,2],[1,4],[2,4],[3,4]]
cnts=np.arange(-2,2.2,0.2)
wdts=np.ones(cnts.shape);
wdts[0]=2.0
wdts[5]=2.0
wdts[10]=4.0
wdts[15]=2.0
wdts[20]=2.0
vcb = np.arange(-12,14,2)
ii=0;
xx,yy=np.meshgrid(np.arange(-6,6,0.1),np.arange(-6,6,0.1))
za=xx+1j*yy;
zb=0.0
b=0.09
for ns in nlist:
    ii=ii+1;
    [n,m]=ns;
    print('computing I'+str(n)+str(m)+' ...\n')
    inm=gp.Inm(za,zb,b,n,m)
    plt.subplot(3,3,ii)
    plt.pcolormesh(np.real(za),np.imag(za),np.real(inm),shading='gouraud',rasterized=True)
    plt.clim(-10,10)
    cs=plt.contour(np.real(za),np.imag(za),np.real(inm),cnts,colors='k',linewidths=wdts*0.5)
    aa=r'$Re[I_{'+str(n)+str(m)+r'}(\zeta_\alpha,\zeta_\beta,b)]$'
    plt.text(-1,3,aa,fontsize=14)
    if (m!=4 and n!=1):
        plt.tick_params(labelbottom='off',labelleft='off')
    elif(n==1 and m==4) :
        plt.text(2,-7.5,r'$Re[\zeta_\alpha]$',fontsize=14)
        plt.text(-7.5,3,r'$Im[\zeta_\alpha]$',fontsize=14,rotation=90)
        plt.tick_params(axis='both', which='major', labelsize=9,labelbottom='on',labelleft='on')
    elif(m==4) :
        plt.text(2,-7.5,r'$Re[\zeta_\alpha]$',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=9,labelbottom='on',labelleft='off')
    elif(n==1): 
        plt.text(-7.5,3,r'$Im[\zeta_\alpha]$',fontsize=14,rotation=90)
        plt.tick_params(axis='both', which='major', labelsize=9,labelleft='on',labelbottom='off')
plt.subplots_adjust(wspace=0.02,hspace=0.02,left=0.05,right=0.92,top=0.95,bottom=0.05)
fig=plt.gcf()
cax = fig.add_axes([0.93, 0.1, 0.02, 0.80])
plt.figure(2)
s=plt.contourf(np.real(za),np.imag(za),np.real(inm)*(np.real(inm)<14)*(np.real(inm)>-14),100)
plt.clim(-10,10)
plt.close()
plt.figure(1)
cbar=plt.colorbar(s, cax=cax,ticks=vcb)
cbar.ax.tick_params(labelsize=8) 
cbar.solids.set_rasterized(True) 
fig.set_size_inches(8,8)
plt.savefig("figures_re.pdf")
