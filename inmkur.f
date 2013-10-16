
      subroutine epskur(omr,omi,pars,u,v,flag)
      double precision omr,omi,ky,kpar,omdi,tau,etai,u,v,
     *     minomdlim,sqrttwo,omsi,epsabs,epsrel,alim,blim,
     *     abserr,bbi,work(40000),resr,resi,zbb
      double complex zaa,i,w,om
      double precision resFepspd_im,resFepspd_re,pars(5),
     *     Fepskur_im,Fepskur_re
      integer nlimit,last,neval,ier,iwork(10000)
      logical flag
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      PARAMETER (MINOMDLIM = -1e-6,
     *     nlimit=10000,
     *     sqrttwo =1.414213562373095 )
      FLAG=.FALSE.
      om=cmplx(omr,omi)
      omdi=pars(1)
      etai=pars(2)
      tau=pars(3)
      ky=pars(4)
      kpar=pars(5)

      if(omdi.LT.MINOMDLIM) then
         zaa=-0.5*cmplx(omr,omi)/omdi
         zbb=kpar/sqrttwo/omdi
         w=zbb**2/4-zaa
         bbi=ky**2;
         omsi=-ky;
         epsrel = 1.0e-2
         epsabs = 1.0e-6
         Alim=-1.0
         Blim=1.0
         CALL DQAG(Fepskur_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         CALL DQAG(Fepskur_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         u=resr+1.0+1.0/tau;
         v=resi;
      endif

      RETURN
      end subroutine epskur

      double precision function Fepskur_re(mu)
      double complex Fepskur_int_re
      external Fepskur_int_re
      double precision mu,muf
      integer neval,ier,nlimit,last,iwork(10000)
      double precision alim,blim,epsabs,epsrel,resr,abserr,work(40000)
      common /fkurcom/ muf
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      muf=mu
      nlimit=10000
      Alim=0.0
      CALL DQAGI(Fepskur_int_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *     ier,nlimit,40000,last,iwork,work)
      Fepskur_re=resr
      return
      end function Fepskur_re

      double precision function Fepskur_im(mu)
      double complex Fepskur_int_im
      external Fepskur_int_im
      double precision mu,muf
      integer neval,ier,nlimit,last,iwork(10000)
      double precision alim,blim,epsabs,epsrel,resi,abserr,work(40000)
      common /fkurcom/ muf
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      muf=mu
      nlimit=10000
      Alim=0.0
      CALL DQAGI(Fepskur_int_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *     ier,nlimit,40000,last,iwork,work)
      Fepskur_im=resi
      return
      end function Fepskur_im

      double precision function Fepskur_int_im(r)
      double complex Fepskur,res
      external Fepskur
      double precision r,muf
      common /fkurcom/ muf
      res=Fepskur(r)
      Fepskur_int_im=dimag(res)
      return
      end function Fepskur_int_im

      double precision function Fepskur_int_re(r)
      double complex Fepskur,res
      external Fepskur
      double precision muf,r
      common /fkurcom/ muf
      res=Fepskur(r)
      Fepskur_int_re=dble(res)
      return
      end function Fepskur_int_re

      double complex function Fepskur(r)
      double precision limsinglg,xbr,xbi,Jr0,Ji0,zbb,bbi,sqrtpi,
     *     omdi,omsi,tau,etai,ky,kpar
      integer mf,nf,ierr,nz
      double complex zaa,i,w,xb,J0,res,om
      double precision r,muf
      common /fkurcom/ muf
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      parameter (sqrtpi = 1.77245385090552)
      xb=2.0*dsqrt(bbi*(1-muf**2))*r
      xbr=dble(xb)
      xbi=dimag(xb)
      i=cmplx(0,1)
      call zbesj(xbr,xbi,0,1,1,Jr0,Ji0,nz,ierr)
      J0=cmplx(Jr0,Ji0)
      If(dble(w).ne.0) then
         res=2.0/omdi*J0**2/sqrtpi*r**2/(r-zsqrt(w))/(r+zsqrt(w))*
     *        dexp(-(muf*r+0.5*zbb)**2-2.0*(1.0-muf**2)*r**2)*
     *        (om-omsi*(1.0-1.5*etai)-2*(1-muf**2)**2*r**2*omsi*etai
     *        -(muf*r+zbb*0.5)**2*omsi*etai)
      else
         res=0.0
      end if
c      write (*,*) res
      Fepskur=res;
      return
      end function Fepskur


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
      Alim=-1.0
      Blim=1.0
      CALL DQAG(Fkur_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *     neval,ier, nlimit,40000,last,iwork,work)
      CALL DQAG(Fkur_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *     neval,ier, nlimit,40000,last,iwork,work)
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

      double precision function Fkur_re(mu)
      double complex Fkur_int_re,zaa,w
      external Fkur_int_re
      double precision mu,muf,zbb,bbi,limsinglg
      integer neval,ier,nlimit,last,iwork(10000),mf,nf,npts2,spoints(3)
      double precision alim,blim,epsabs,epsrel,resr,abserr,work(40000)
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      common /fkurcom/ muf
      PARAMETER (limsinglg=1.0e-5)
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      muf=mu
      nlimit=10000
      if(dabs(dimag(w)).GT.limsinglg.OR.dble(w).LT.0) then
         Alim=0.0
         CALL DQAGI(Fkur_int_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         Fkur_re=resr
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_Re (in if):"
            write (*,*) "w:", w, "ier:", ier
         end if
      else
         alim=0.0
         blim=dsqrt(1.01*dble(w))
         spoints(1)=dsqrt(dble(w))
         npts2=3
         CALL DQAGP(Fkur_int_re,alim,blim,npts2,spoints,epsabs,epsrel,
     *        resr,abserr,neval,ier,nlimit,40000,last,iwork,work)
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_Re (in else, DQAGP):"
            write (*,*) "w:", w, "ier:", ier
         end if
         Fkur_re=resr;
         alim=blim;
         CALL DQAGI(Fkur_int_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         Fkur_re=Fkur_re+resr;
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_Re (in else, DQAGI):"
            write (*,*) "w:", w, "ier:", ier
         end if
      endif
      return
      end function Fkur_re

      double precision function Fkur_int_re(r)
      double complex Fkur,res
      external Fkur
      double precision muf,r
      common /fkurcom/ muf
      res=Fkur(r)
      Fkur_int_re=dble(res)
      return
      end function Fkur_int_re

      double precision function Fkur_im(mu)
      double complex Fkur_int_im,zaa,w
      external Fkur_int_im
      double precision mu,muf,zbb,bbi,limsinglg
      common /fkurcom/ muf
      integer neval,ier,nlimit,last,iwork(10000),mf,nf,npts2,spoints(3)
      double precision alim,blim,epsabs,epsrel,resi,abserr,work(40000)
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      PARAMETER (limsinglg=1.0e-5)
      epsrel = 1.0e-2
      epsabs = 1.0e-6
      muf=mu
      nlimit=10000
      Alim=0.0
      if(dabs(dimag(w)).GT.limsinglg.OR.dble(w).LT.0) then
         CALL DQAGI(Fkur_int_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         Fkur_im=resi
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_im (in if):"
            write (*,*) "w:", w, "ier:", ier
         end if
      else
         alim=0.0
         blim=dsqrt(1.01*dble(w))
         npts2=3
         spoints(1)=dsqrt(dble(w))
         epsrel = 1.0e-4
         epsabs = 1.0e-12
         CALL DQAGP(Fkur_int_im,alim,blim,npts2,spoints,epsabs,epsrel,
     *        resi,abserr,neval,ier,nlimit,40000,last,iwork,work)
         Fkur_im=resi
         alim=blim;
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_im (in else DQAGP):"
            write (*,*) "w:", w, "ier:", ier
         end if
         CALL DQAGI(Fkur_int_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         if (ier.GT.0) then
            write(*,*) "problem in Fkur_im (in else DQAGI):"
            write (*,*) "w:", w, "ier:", ier
         end if
         Fkur_im=Fkur_im+resi;
      endif
      return
      end function Fkur_im

      double precision function Fkur_int_im(r)
      double complex Fkur,res
      external Fkur
      double precision r,muf
      common /fkurcom/ muf
      res=Fkur(r)
      Fkur_int_im=dimag(res)
      return
      end function Fkur_int_im

      double complex function Fkur(r)
      double precision limsinglg,xbr,xbi,Jr0,Ji0,zbb,bbi,sqrtpi
      integer mf,nf,ierr,nz
      double complex zaa,i,w,xb,J0,res
      double precision r,muf
      common /fkurcom/ muf
      common /inmcom/ mf,nf,zbb,bbi,zaa,w
      parameter (limsinglg = 1.0e-5, sqrtpi = 1.77245385090552)
      xb=2.0*dsqrt(bbi*(1-muf**2))*r
      xbr=dble(xb)
      xbi=dimag(xb)
      i=cmplx(0,1)
      call zbesj(xbr,xbi,0,1,1,Jr0,Ji0,nz,ierr)
      J0=cmplx(Jr0,Ji0)
      If(dble(w).ne.0) then
         res=2.0**(0.5*(nf+3))*J0**2/sqrtpi*r**(nf+1)*
     *        (muf*r+zbb*0.5)**mf*
     *        (1-muf**2)**(0.5*(nf-1))/(r-zsqrt(w))/(r+zsqrt(w))*
     *        dexp(-(muf*r+0.5*zbb)**2-2.0*(1.0-muf**2)*r**2)
      else
         res=0.0
      end if
c      write (*,*) res
      Fkur=res;
      return
      end function Fkur
