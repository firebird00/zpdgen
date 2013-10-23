C      INMZPD
      subroutine inmzpd (xi, yi, zb, bi, n, m, nw, u, v, flag)
Cf2py double precision intent(out) u,v
Cf2py integer intent(out) flag
C
C  PARAMETER LIST
C     XI     = REAL      PART OF ZA
C     YI     = IMAGINARY PART OF ZA
C     BI     = BI
C     ZB    = ZB
C     N     = N
C     M     = M
C     U      = REAL      PART OF INM(Z)
C     V      = IMAGINARY PART OF INM(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI, BI, ZB, N,M      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C

      implicit none
      DOUBLE PRECISION xi,yi,bi,zb,z1,z2,u,v, 
     *     epsabs,epsrel,alim,blim,abserr,
     *     work(40000),resr,resi,zbb,bbi,limsingsm
      double complex zaa,za,i,w
      INTEGER n,m,np1,nu,j,l,iwork(10000),nweid,nw,neval
      LOGICAL A, B, FLAG
      integer nlimit,mf,nf,last,npts2,spoints(3),ier
      PARAMETER (nlimit=10000,
     *     epsrel=1.0e-2,epsabs=1.0e-6,npts2=3,
     *     limsingsm=1.0e-8)
      double precision fpd_re,fpd_im
      EXTERNAL fpd_re,fpd_im,resfpd_im,resfpd_re
      external dqagi,dqagp,dqag,prerr
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      common /nweid/ nweid
      FLAG = .FALSE.
      i=cmplx(0,1)
      nweid=nw
      za=XI+i*YI
      zaa=za
      zbb=zb
      bbi=bi
      mf=m
      nf=n
      w=zbb**2/4-zaa
      if(dabs(dimag(w)).GT.limsingsm.OR.dble(w).LT.0) then
      Alim=0.0
         CALL DQAGI(Fpd_re,alim,1,epsabs,epsrel,resr,abserr,neval,ier,
     *        nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAGI(Fpd_im,alim,1,epsabs,epsrel,resi,abserr,neval,ier,
     *        nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=resr;
         v=resi;
      else
         alim=0.0
         blim=dsqrt(2.01*dble(w))
         spoints(1)=dsqrt(2.0*dble(w))
         CALL DQAGP(Fpd_re,alim,blim,npts2,spoints,epsabs,epsrel,resr,
     *        abserr,neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAGP(Fpd_im,alim,blim,npts2,spoints,epsabs,epsrel,resi,
     *        abserr,neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=resr;
         v=resi;
         alim=blim;
         CALL DQAGI(Fpd_re,alim,1,epsabs,epsrel,resr,abserr,neval,ier,
     *        nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAGI(Fpd_im,alim,1,epsabs,epsrel,resi,abserr,neval,ier,
     *        nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u+resr;
         v=v+resi;
      endif
      if(dimag(zaa).LT.0.AND.dble(w).GT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAG(resFpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u-resr
         v=v-resi
      else if(dimag(zaa).EQ.0.AND.dble(w).GT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAG(resFpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u-0.5*resr
         v=v-0.5*resi
      endif
      RETURN
*
 100  FLAG = .TRUE.
      call prerr(zaa,zbb,bbi)
      RETURN
*
      END

      double precision function Fpd_re(s)
      double complex Fpd
      external Fpd
      double precision s
      Fpd_re=dble(Fpd(s))
      return
      end function Fpd_re

      double precision function Fpd_im(s)
      double complex Fpd
      external Fpd
      double precision s
      Fpd_im=dimag(Fpd(s))
      return
      end function Fpd_im

      double complex function Fpd(s)
      double precision s,limsingsm,xbr,Jr0,Ji0,zbb,bbi
      integer mf,nf,ierr,nz
      double complex z1,z2,zaa,Gm,weidGm,w
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      external weidGm
      parameter (limsingsm = 1.0e-8)
      z1=0.5*(zbb+zsqrt(zbb**2-2.0*(s**2+2.0*zaa)))
      z2=0.5*(zbb-zsqrt(zbb**2-2.0*(s**2+2.0*zaa)))
      Gm=weidGm(z1,z2,mf)
      if ((zabs(z1).LT.limsingsm).and.(zabs(z2).LT.limsingsm)) then
         Fpd=0.0
      else
         xbr=dble(bbi*2.0)**(0.5)*s
c         call zbesj(xbr,0.0,0,1,1,Jr0,Ji0,nz,ierr)
         JR0=DBESJ0(XBR)
         Fpd=2.0*dexp(-s**2)*Jr0**2*Gm*s**nf
      endif
      return
      end function Fpd

      double complex function weidGm(z1,z2,m)
      double complex weidZm,z1,z2
      double precision adelz,limsingsm
      integer m
      parameter (limsingsm = 1.0e-8)
      external weidZm
      adelz=zabs(z1-z2)
      if(adelz.LT.limsingsm) then
         if(m.gt.0) then
            weidGm=m*weidZm(z1,m-1)-2*weidZm(z1,m+1)
         else
            if (m.eq.0) then
               weidGm=-2*weidZm(z1,m+1)
            endif
         endif
      else
         weidGm=(weidZm(z1,m)-weidZm(z2,m))/(z1-z2)
      endif
      return
      end function weidGm

      double complex function weidZm(z,m)
      double complex z,i,Z0
      double precision xi,yi,u,v,flag,dgamma,sqrtpi
      integer m,k,nweid
      external wofzwh16,wofzwh8,wofzwh4,wofzh,dgamma
      parameter (sqrtpi = 1.77245385090552)
      common /nweid/ nweid
      i=cmplx(0,1)
      xi=dble(z)
      yi=dimag(z)
      select case (nweid)
      case(1)
         call wofzh(xi,yi,u,v,flag)
      case(4)
         call wofzwh4(xi,yi,u,v,flag)
      case(8)
         call wofzwh8(xi,yi,u,v,flag)
      case(16)
         call wofzwh16(xi,yi,u,v,flag)
      case default
         call wofzwh16(xi,yi,u,v,flag)
      end select
      Z0=u+i*v
      weidZm=i*sqrtpi*Z0*z**m
      if (m.gt.0) then
         do 20 k=1,m
            weidZm=weidZm+1.0/sqrtpi*dgamma((m-k+1)*0.5d0)*z**(k-1)
 20      continue
      endif 
      return
      end function weidZm

      double precision function resFpd_re(s)
      double complex resFpd
      external resFpd
      double precision s
      resFpd_re=dble(resFpd(s))
      return
      end function resFpd_re

      double precision function resFpd_im(s)
      double complex resFpd
      external resFpd
      double precision s
      resFpd_im=dimag(resFpd(s))
      return
      end function resFpd_im

      double complex function resFpd(mu)
      double precision mu,xbr,xbi,Jr0,Ji0,zbb,bbi,sqrtpi
      integer mf,nf,ierr,nz
      double complex zaa,i,w,xb,J0
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      parameter (sqrtpi = 1.77245385090552)
      xb=2.0*zsqrt(bbi*(1-mu**2)*w)
      xbr=dble(xb)
      xbi=dimag(xb)
      i=cmplx(0,1)
      call zbesj(xbr,xbi,0,1,1,Jr0,Ji0,nz,ierr)
      J0=cmplx(Jr0,Ji0)
      resFpd=i*2**(0.5*(nf+3))*J0**2*sqrtpi*w**(nf*0.5)*
     *     (1-mu**2)**((nf-1)*0.5)*(mu*zsqrt(w)+0.5*zbb)**mf*
     *     zexp(-(mu*zsqrt(w)+0.5*zbb)**2-2.0*(1.0-mu*mu)*w)
c      resFpd=i*2**(0.5*(nf+3))*J0**2*sqrtpi*w**(nf*0.5)*
c     *     zexp(-(mu*zsqrt(w)+0.5*zbb)**2-2.0*(1.0-mu*mu)*w)
c      if(nf.gt.1) then
c         resFpd=resFpd*(1-mu**2)**(0.5*(nf-1))
c      endif
c      if(mf.gt.0) then
c         resFpd=resFpd*(mu*zsqrt(w)+0.5*zbb)**mf
c      endif
      return
      end function resFpd

C     inm
      subroutine Inmweid(za,zb,b,n,m,nw,numza,res)
Cf2py double complex dimension(numza) :: za
Cf2py double precision :: b
Cf2py double precision :: zb
Cf2py  integer :: n
Cf2py  integer :: m
Cf2py  integer :: nw
Cf2py  integer, optional,depend(za) :: numza=len(za)
Cf2py double complex dimension(numza) :: res
      double precision b,zb,xi,yi,u,v
      double complex za(numza),res(numza)
      integer n,m,k,numza,nw
      logical flag
      external inmzpd,prerr
      double precision largestxi,largestyi,largeval
      parameter (largestxi=1E4,largestyi=1E4,largeval=1d120)
      do 10 k=1,numza
         xi=dble(za(k))
         yi=dimag(za(k))
         if (dabs(xi).gt.largestxi.or.dabs(yi).gt.largestyi) goto 120
         call inmzpd(xi,yi,zb,b,n,m,nw,u,v,flag)
         if(flag) goto 110
         res(k)=cmplx(u,v)
 10   continue
      return
 110  call prerr(za(k),zb,b)
      return
 120  res(k)=largeval
      continue
      end

C     inm
      subroutine Gmweid(z1,z2,m,numza,res)
Cf2py double complex dimension(numza) :: z1
Cf2py double complex dimension(numza) :: z2
Cf2py  integer :: m
Cf2py  integer, optional,depend(z1) :: numza=len(z1)
Cf2py double complex dimension(numza) :: res
      double complex z1(numza),z2(numza),res(numza)
      external weidGm
      double complex weidGm
      integer m,k,numza
      do 30 k=1,numza
         res(k)=weidGm(z1(k),z2(k),m)
 30   continue
      return
      end
