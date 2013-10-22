      program timetest
      double precision b,zb,xi,yi,u,v
      double complex za,res,i,zamin,zamax,dza,ttf(1000,1000)
      external inmzpd
      integer n,m,numx,numy,lx,ly,k,ns(3,2)
      logical flag
      real finish,start
      data ns/1,1,3,0,2,0/
      i=cmplx(0,1)
      zamin=-6.0-6.0*i
      zamax=6.0+6.0*i
      dza=0.1+0.1*i
      zb=0.0
      b=0.09

      numx=nint(dble(zamax-zamin)/dble(dza))
      numy=nint(dimag(zamax-zamin)/dimag(dza))

      do 20 k=1,3
         n=ns(k,1)
         m=ns(k,2);
         print '("computing I",i1,i1," for ", i3," x",i3,
     *" points...")',n,m,numx+1,numy+1
         call cpu_time(start)
         do 10 lx=0,numx
            do 10 ly=0,numy
               za=zamin+dble(dza)*lx+i*dimag(dza)*ly
               xi=dble(za)
               yi=dimag(za)
               call inmzpd(xi,yi,zb,b,n,m,u,v,flag)
               if(flag) goto 100
               ttf(lx+1,ly+1)=cmplx(u,v);
 10         continue
         call cpu_time(finish)
         print '("Time = ",f12.6," seconds.")',finish-start
 20   continue
      stop

 100  write(*,*) "error!"
      stop
      end

