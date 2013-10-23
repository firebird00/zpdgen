      program test
      double precision b,zb,xi,yi,u,v
      double complex za,res,i
      external inmzpd
      integer n,m
      logical flag
      real finish,start
      i=cmplx(0,1)
      za=-6.0+6.0*i
      zb=0.0
      b=0.09
      n=1
      m=0
      print '("computing I",i1,i1,".")',n,m
      xi=dble(za)
      yi=dimag(za)
      call cpu_time(start)
      call inmzpd(xi,yi,zb,b,n,m,u,v,flag)
      if(flag) goto 100
      call cpu_time(finish)
      res=cmplx(u,v)
      print '("Time = ",f12.6," seconds.")',finish-start
      write (*,*),"res:",res
      stop

 100  write(*,*) "error!"
      stop
      end

