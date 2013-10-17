      program timetest
      double precision b,zb,xi,yi,u,v,omdi,etai,tau,ky,kpar,pars(5)
      double precision omsi
      double complex za,res,i,zamin,zamax,dza,ttf(1000,1000)
      double complex i10,i12,i30,om
      external inmzpd
      external epszpd
      external epskur
      integer n,m,numx,numy,lx,ly,k,ns(3,2),numt,nws(4),nw,l
      logical flag
      real finish,start
      data ns/1,1,3,0,2,0/
      data nws/1,4,8,16/
      i=cmplx(0,1)
      zamin=-6.0-6.05*i
      zamax=6.0+5.95*i
      dza=0.1+0.1*i
      zb=0.0
      b=0.09
      numx=nint(dble(zamax-zamin)/dble(dza))
      numy=nint(dimag(zamax-zamin)/dimag(dza))
      omdi=-0.5
      etai=3.4
      tau=1.0
      ky=0.3
      kpar=0.0
      pars=(/omdi,etai,tau,ky,kpar/)
      omsi=-ky

      do 50 l=1,4
         nw=nws(l)
         do 20 k=1,5
            if(k.LT.4) then
               n=ns(k,1)
               m=ns(k,2);
               print '("computing I",i1,i1," for ", i3," x",i3,
     *" points, nw=",i2,"...")',n,m,numx+1,numy+1,nw
            else if (k.EQ.4) then
               print '("computing eps(om,k) [using Inms] for "
     *, i3," x",i3," points, nw=",i2,"...")',numx+1,numy+1,nw
            else if (k.EQ.5) then
               print '("computing eps(om,k) [combined] for ", i3," x",i3,
     *" points, nw=",i2,"...")',numx+1,numy+1,nw
            endif
            call cpu_time(start)
            
            do 10 lx=0,numx
               do 10 ly=0,numy
                  za=zamin+dble(dza)*lx+i*dimag(dza)*ly
                  xi=dble(za)
                  yi=dimag(za)
                  if(k.LT.4) then
                     call inmzpd(xi,yi,zb,b,n,m,nw,u,v,flag)
                     ttf(lx+1,ly+1)=cmplx(u,v);
                  else if (k.EQ.4) then
                     call inmzpd(xi,yi,zb,b,1,0,nw,u,v,flag)
                     i10=cmplx(u,v)
                     call inmzpd(xi,yi,zb,b,1,2,nw,u,v,flag)
                     i12=cmplx(u,v)
                     call inmzpd(xi,yi,zb,b,3,0,nw,u,v,flag)
                     i30=cmplx(u,v)
                     om=-2*za*omdi
                     ttf(lx+1,ly+1)=1.0+1.0/tau+0.5/omdi*(
     *                    i10*(om-omsi*(1.0-1.5*etai))
     *                    -omsi*etai*(i30+i12))
                  else if (k.EQ.5) then
                     call epszpd(xi,yi,pars,nw,u,v,flag)
                     ttf(lx+1,ly+1)=cmplx(u,v);
                  endif
 10       continue
          call cpu_time(finish)
          print '("Time = ",f12.6," seconds.")',finish-start
 20    continue
 50    continue

      do 40 k=1,5
         if(k.LT.4) then
            n=ns(k,1)
            m=ns(k,2);
            print '("computing I",i1,i1," for ", i3," x",i3,
     *" points...")',n,m,numx+1,numy+1
         else if (k.EQ.4) then
            print '("computing eps(om,k) [using Inms] for ", i3," x",i3,
     *" points...")',numx+1,numy+1
         else if (k.EQ.5) then
            print '("computing eps(om,k) [combined] for ", i3," x",i3,
     *" points...")',numx+1,numy+1
         endif
         call cpu_time(start)
         do 30 lx=0,numx
            do 30 ly=0,numy
               za=zamin+dble(dza)*lx+i*dimag(dza)*ly
               xi=dble(za)
               yi=dimag(za)
               if(k.LT.4) then
                  call inmkur(xi,yi,zb,b,n,m,u,v,flag)
                  ttf(lx+1,ly+1)=cmplx(u,v);
               else if (k.EQ.4) then
                  call inmkur(xi,yi,zb,b,1,0,u,v,flag)
                  i10=cmplx(u,v)
                  call inmkur(xi,yi,zb,b,1,2,u,v,flag)
                  i12=cmplx(u,v)
                  call inmkur(xi,yi,zb,b,3,0,u,v,flag)
                  i30=cmplx(u,v)
                  om=-2*za*omdi
                  ttf(lx+1,ly+1)=1.0+1.0/tau+0.5/omdi*(
     *                 i10*(om-omsi*(1.0-1.5*etai))
     *                 -omsi*etai*(i30+i12))
               else if (k.EQ.5) then
                  call epskur(xi,yi,pars,u,v,flag)
                  ttf(lx+1,ly+1)=cmplx(u,v);
               endif
 30         continue
      call cpu_time(finish)
      print '("Time = ",f12.6," seconds.")',finish-start
 40         continue
      stop
      end
