import numpy as np
import matplotlib.pyplot as plt
import gpdf as gp
import time
nlist=[[1,0],[1,2],[3,0]]
ii=0;
xx,yy=np.meshgrid(np.arange(-6,6,0.1),np.arange(-6,6,0.1))
za=xx+1j*yy;
zb=0.0
b=0.09
for ns in nlist:
    ii=ii+1;
    [n,m]=ns;
    print('computing I'+str(n)+str(m)+' ...')
    t0 = time.clock()
    inm=gp.Inm(za,zb,b,n,m)
    print(time.clock()-t0,"seconds")
