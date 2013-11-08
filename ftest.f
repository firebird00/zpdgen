      program timetest
      double precision b,zb,xi,yi,u,v,pars(5),omdi,etai,tau,
     *     ky,kpar,omsi,LnbyR
      double complex za,res,i,zamin,zamax,dza,ttf(1000,1000)
      external inmzpd,sigmazpd
      integer n,m,numx,numy,lx,ly,k,ns(3,2),numn,numm
      logical flag
      real finish,start
      double complex anm(4,3),om
      data ns/1,1,3,0,2,0/
      data anm/ 12 * 0.0/
      i=cmplx(0,1)
      zamin=-6.0-6.0*i
      zamax=6.0+6.0*i
      dza=0.1+0.1*i
      numx=nint(dble(zamax-zamin)/dble(dza))
      numy=nint(dimag(zamax-zamin)/dimag(dza))

      numn=4
      numm=3

      LnbyR=0.2
      etai=2.5
      tau=1.0
      ky=0.3
      kpar=0.1
      omsi=-ky
      omdi=2*omsi*LnbyR

      pars(1)=omdi
      pars(2)=etai
      pars(3)=tau
      pars(4)=ky
      pars(5)=kpar

      b=ky**2
      zb=-dsqrt(2.0d0)*kpar/omdi

c      om = cmplx(0.20400002598762512,1.0728836485895954E-008)
c      za=-om/omdi
c      za=dcmplx( -3.0000000000000107, -2.7000000000000117)  
c      zb=0.0000000000000000
c      b=8.9999999999999997d-002
c      xi=dble(za)
c      yi=dimag(za)
c      n=0
c      m=0
c      call inmzpd(xi,yi,zb,b,n,m,u,v,flag)
c      stop
      
c      om = cmplx(0.20400002598762512,1.0728836485895954E-008)
c      xi=dble(om)
c      yi=dimag(om)
c      call epszpd(xi,yi,pars,u,v,flag)
c      stop

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

      print '("computing epscombined for ", i3," x",i3,
     *" points...")',numx+1,numy+1
      call cpu_time(start)
      do 30 lx=0,numx
         do 30 ly=0,numy
            za=zamin+dble(dza)*lx+i*dimag(dza)*ly
            om=-za*omdi
            xi=dble(om)
            yi=dimag(om)
            call epszpd(xi,yi,pars,u,v,flag)
            if(flag) goto 100
            ttf(lx+1,ly+1)=cmplx(u,v);
 30   continue
      call cpu_time(finish)
      print '("Time = ",f12.6," seconds.")',finish-start

      print '("computing epscombined (mat) for ", i3," x",i3,
     *" points...")',numx+1,numy+1
      anm(4,1)=-omsi*etai/omdi
      anm(2,3)=-omsi*etai/omdi
      call cpu_time(start)
      do 40 lx=0,numx
         do 40 ly=0,numy
            za=zamin+dble(dza)*lx+i*dimag(dza)*ly
            om=-za*omdi
            xi=dble(om)
            yi=dimag(om)
c     Note the fortran arrays start from 1 not zero...
            anm(2,1)=(om-omsi*(1-1.5*etai))/omdi
            call sigmazpd(xi,yi,zb,b,anm,numn,numm,u,v,flag)
            if(flag) goto 100
            ttf(lx+1,ly+1)=cmplx(u,v);
 40   continue
      call cpu_time(finish)
      print '("Time = ",f12.6," seconds.")',finish-start
      stop

 100  write(*,*) "error!"
      stop
      end

