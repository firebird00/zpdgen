from numpy import reshape,shape,transpose
from inmzpd import inmweid

def inm(za,zb,b,n,m):
    res=inmweid(za,zb,b,n,m)
    if (shape(za) ==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(za)))
    res=transpose(res)
    return res
