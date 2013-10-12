C      INMZPD
      subroutine inmkur (xi, yi, zb, bi, n, m, u, v, flag)
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
      DOUBLE PRECISION xi,yi,bi,zb,z1,z2,c,daux,x,y,xabs,yabs,u,v, 
     *     factor,rmaxreal, rmaxexp, rmaxgoni, h,h2, kapn, qlambda,
     *     qrho, rx,ry,sx,sy,tx,ty, u1, u2, v1,v2,xabsq,xaux,xquad,xsum,
     *     yquad, ysum,epsabs,epsrel,alim,blim,key,result,abserr,
     *     limit,lenw, work(40000),resr,resi,zbb,bbi
      double complex zaa,za,i,w
      INTEGER n,m,np1,nu,j,l,iwork(10000),neval,ier
      LOGICAL A, B, FLAG
      integer nlimit,mf,nf,last
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *     RMAXREAL = 0.5D+154,
     *     RMAXEXP  = 708.503061461606D0,
     *     RMAXGONI = 3.53711887601422D+15,
     *     nlimit=10000)
      double precision fkur_re,fkur_im
      EXTERNAL DQAG, fkur_re,fkur_im,resfpd_im,resfpd_re
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      FLAG = .FALSE.
      i=cmplx(0,1)
      za=XI+i*YI
      zaa=za
      zbb=zb
      bbi=bi
      mf=m
      nf=n
      w=zbb**2/4-zaa
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      Alim=0.0
c      CALL DQAG(Fpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,neval,ier,
c     *     nlimit,40000,last,iwork,work)
      CALL DQAGI(Fkur_re,alim,1,epsabs,epsrel,resr,abserr,neval,ier,
     *     nlimit,40000,last,iwork,work)
c      write (*,*) resr

c      CALL DQAG(Fpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,neval,ier,
c     *     nlimit,40000,last,iwork,work)
      CALL DQAGI(Fkur_im,alim,1,epsabs,epsrel,resi,abserr,neval,ier,
     *     nlimit,40000,last,iwork,work)

      u=resr;
      v=resi;

      if(dimag(zaa).LT.0.AND.dble(zaa).LT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         CALL DQAG(resFpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         u=u-resr
         v=v-resi
      else if(dimag(zaa).EQ.0.AND.dble(zaa).LT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         CALL DQAG(resFpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         u=u-resr*0.5
         v=v-resi*0.5
      endif
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END

      double precision function Fkur_re(r)
      double complex Fkur_int_re
      external Fkur_int_re
      double precision r,rf
      common /fkurcom/ rf
      integer neval,ier,nlimit,last,iwork(10000)
      double precision alim,blim,epsabs,epsrel,resr,abserr,work(40000)
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      Alim=-1.0
      Blim=1.0
      rf=r
       nlimit=10000
      CALL DQAG(Fkur_int_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *     neval,ier, nlimit,40000,last,iwork,work)
      Fkur_re=resr
      return
      end function Fkur_re

      double precision function Fkur_int_re(mu)
      double complex Fkur,res
      external Fkur
      double precision rf,mu
      common /fkurcom/ rf
      res=Fkur(mu)
c      write(*,*) res
      Fkur_int_re=dble(res)
      return
      end function Fkur_int_re

      double precision function Fkur_im(r)
      double complex Fkur_int_im
      external Fkur_int_im
      double precision r,rf
      common /fkurcom/ rf
      integer neval,ier,nlimit,last,iwork(10000)
      double precision alim,blim,epsabs,epsrel,resi,abserr,work(40000)
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      Alim=-1.0
      Blim=1.0
      rf=r
       nlimit=10000
      CALL DQAG(Fkur_int_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *     neval,ier, nlimit,40000,last,iwork,work)
      Fkur_im=resi
      return
      end function Fkur_im

      double precision function Fkur_int_im(mu)
      double complex Fkur,res
      external Fkur
      double precision rf,mu
      common /fkurcom/ rf
      res=Fkur(mu)
c      write (*,*) res
c      Fkur_int_im=dimag(Fkur(mu))
      Fkur_int_im=dimag(res)
      return
      end function Fkur_int_im

      double complex function Fkur(mu)
      double precision limsinglg,xbr,xbi,Jr0,Ji0,zbb,bbi,sqrtpi
      integer mf,nf,ierr,nz
      double complex zaa,i,w,xb,J0,res
      double precision rf,mu
      common /fkurcom/ rf
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      parameter (limsinglg = 1.0e-5, sqrtpi = 1.77245385090552)
      xb=2.0*dsqrt(bbi*(1-mu**2))*rf
      xbr=dble(xb)
      xbi=dimag(xb)
      i=cmplx(0,1)
      call zbesj(xbr,xbi,0,1,1,Jr0,Ji0,nz,ierr)
      J0=cmplx(Jr0,Ji0)
      If(dble(w).ne.0) then
         res=2.0**(0.5*(nf+3))*J0**2/sqrtpi*rf**(nf+1)*
     *        (mu*rf+zbb*0.5)**mf*
     *        (1-mu**2)**(0.5*(nf-1))/(rf-zsqrt(w))/(rf+zsqrt(w))*
     *        dexp(-(mu*rf+0.5*zbb)**2-2.0*(1.0-mu**2)*rf**2)
      else
         res=0.0
      end if
c      write (*,*) res
      Fkur=res;
      return
      end function Fkur
