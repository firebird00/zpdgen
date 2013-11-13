C       COMPUTATION OF SPECIAL FUNCTIONS
C 
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs. 
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
C      Changed according to ERRATA also.
C      
C      Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
C
C       **********************************

        SUBROUTINE CBJ0(Z,CJ0)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z), 
C                Y1(z), and their derivatives for a complex
C                argument
C       Input :  z --- Complex argument
C       Output:  CBJ0 --- J0(z)
C                CDJ0 --- J0'(z)
C                CBJ1 --- J1(z)
C                CDJ1 --- J1'(z)
C                CBY0 --- Y0(z)
C                CDY0 --- Y0'(z)
C                CBY1 --- Y1(z)
C                CDY1 --- Y1'(z)
C       =======================================================
C
        DOUBLE PRECISION A,B,A1,B1,A0,RP2,PI,EL
        DOUBLE COMPLEX CJ0,Z,CI,Z1,Z2,CR,CT1,CT0,CQ0,CP0,CU
        DIMENSION A(12),B(12),A1(12),B1(12)
        INTEGER K,K0
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CJ0=(1.0D0,0.0D0)
           RETURN
        ENDIF
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CJ0=CJ0+CR
              IF (CDABS(CR).LT.CDABS(CJ0)*1.0D-15) GO TO 15
10         CONTINUE
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
        ENDIF
 15     RETURN
        END
