from numpy import reshape,shape,transpose
from inmzpd import inmweid

def Inm(za,zb,b,n,m,mw=16):
    res=inmweid(za,zb,b,n,m,mw)
    if (shape(za) ==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(za)))
    res=transpose(res)
    return res

#def Gm(z1,z2,m):
#    if(shape(z1)!=shape(z2)):
#        print("error: shape(z1)!=shape(z2)!\n");
#        res=()
#        return res;
#    res=gmweid(z1,z2,m)
#    if (shape(z1) ==() and shape(z2)==()):
#        res=res[0]
#        return res
#    res=reshape(res,shape(transpose(z1)))
#    res=transpose(res)
#    return res
