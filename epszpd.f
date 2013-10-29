C      EPSZPD
      subroutine epszpd (omr, omi, pars, nw, u, v, flag)
C  PARAMETER LIST
C     omr     = REAL      PART OF OM
C     omi     = IMAGINARY PART OF OM
C     U      = REAL      PART OF eps(om)
C     V      = IMAGINARY PART OF eps(om)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  OMR, OMI     ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C
      implicit none
      DOUBLE PRECISION omr,omi,ky,kpar,omdi,etai,u,v,
     *     epsabs,epsrel,alim,blim,abserr,
     *     work(40000),resr,resi,zbb,bbi,minomdlim,tau,
     *     pars(5),sqrttwo,omsi
      double complex zaa,i,w,om
      INTEGER n,m,np1,nu,j,l,iwork(10000),limit
      LOGICAL A, B, FLAG
      integer nlimit,mf,nf,last,neval,ier,nw,nweid
      PARAMETER (MINOMDLIM = -1e-6,
     *     nlimit=10000,epsrel=1.0e-2,epsabs=1.0e-6,
     *     sqrttwo =1.414213562373095 )
      double precision fpd_re,fpd_im
      EXTERNAL DQAG, fepspd_re,fepspd_im,resfepspd_im,resfepspd_re
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      common /nweid/ nweid
      FLAG = .FALSE.
      i=cmplx(0,1)
      om=cmplx(omr,omi)
      omdi=pars(1)
      etai=pars(2)
      tau=pars(3)
      ky=pars(4)
      kpar=pars(5)
      nweid=nw
      if(omdi.LT.MINOMDLIM) then
         zaa=-cmplx(omr,omi)/omdi
         zbb=-kpar*sqrttwo/omdi
         w=zbb**2/4-zaa
         bbi=ky**2;
         omsi=-ky;
         Alim=0.0
         CALL DQAGI(Fepspd_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAGI(Fepspd_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=resr+1.0+1.0/tau;
         v=resi;
      end if

      if(omi.LT.0.AND.dble(w).GT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFepspd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAG(resFepspd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u-resr
         v=v-resi
      endif
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END

      double precision function Fepspd_re(s)
      double complex Fepspd
      external Fepspd
      double precision s
      Fepspd_re=dble(Fepspd(s))
      return
      end function Fepspd_re

      double precision function Fepspd_im(s)
      double complex Fepspd
      external Fepspd
      double precision s
      Fepspd_im=dimag(Fepspd(s))
      return
      end function Fepspd_im

      double complex function Fepspd(s)
      double precision s,limsinglg,xbr,Jr0,Ji0,zbb
      double precision bbi,omdi,omsi,etai,tau,ky,kpar
      integer ierr,nz
      double complex z1,z2,zaa,G0,G2,weidGm,w,res,om
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      external weidGm
      parameter (limsinglg = 1.0e-5)
      z1=0.5D0*(zbb+zsqrt(zbb**2-2.0*(s**2+2.0*zaa)))
      z2=0.5D0*(zbb-zsqrt(zbb**2-2.0*(s**2+2.0*zaa)))
      if ((zabs(z1).LT.limsinglg).and.(zabs(z2).LT.limsinglg)) then
         Fepspd=0.0
      else
         G0=weidGm(z1,z2,0)
         G2=weidGm(z1,z2,2)
         xbr=dsqrt(bbi*2.0)*s
         jr0=dbesj0(xbr)
c         call zbesj(xbr,0.0,0,1,1,Jr0,Ji0,nz,ierr)
         Fepspd=2.0d0/omdi*
     *        Jr0**2*dexp(-s**2)*(
     *        (om-omsi*(1.0+(s**2-1.5)*etai))*G0-omsi*etai*G2
     *        )*s
      end if
      return
      end function Fepspd

      double precision function resFepspd_re(s)
      double complex resFepspd
      external resFepspd
      double precision s
      resFepspd_re=dble(resFepspd(s))
      return
      end function resFepspd_re

      double precision function resFepspd_im(s)
      double complex resFepspd
      external resFepspd
      double precision s
      resFepspd_im=dimag(resFepspd(s))
      return
      end function resFepspd_im

      double complex function resFepspd(mu)
      double precision mu,limsinglg,xbr,xbi,Jr0,Ji0,zbb
      double precision bbi,sqrtpi,omdi,etai,tau,ky,kpar,omsi
      integer ierr,nz
      double complex zaa,i,w,xb,J0,om
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      parameter (limsinglg = 1.0e-5, sqrtpi = 1.77245385090552)
      xb=2.0*(bbi*(1.0d0-mu**2)*w)**(0.5)
      xbr=dble(xb)
      xbi=dimag(xb)
      i=cmplx(0,1)
      call zbesj(xbr,xbi,0,1,1,Jr0,Ji0,nz,ierr)
      J0=cmplx(Jr0,Ji0)
      resFepspd=i*4.0*J0**2*sqrtpi/omdi*zsqrt(w)*
     *     zexp(-(mu*zsqrt(w)+0.5*zbb)**2-2.0*(1.0-mu**2)*w)*
     *     (om-omsi*(1.0+etai*(2.0*w*(1-mu**2)
     *     +(mu*zsqrt(w)+0.5*zbb)**2-1.5)))
      return
      end function resFepspd
