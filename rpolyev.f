      SUBROUTINE RPOLYEV(NN,SR,SI,P,QR,QI,PVR,PVI)                   
C EVALUATES A REAL POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
C PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
<<<<<<< HEAD
      implicit none
=======
>>>>>>> 6df2dda88b5c6a538b0f2354ef01e1bc5a2597ff
      INTEGER NN,I
      DOUBLE PRECISION P(NN),QR(NN),QI(NN),
     *    SR,SI,PVR,PVI,T
      QR(1) = P(1)
      QI(1) = 0
      PVR = QR(1)
      PVI = QI(1)
      DO 10 I = 2,NN
          T = PVR*SR-PVI*SI+P(I)
          PVI = PVR*SI+PVI*SR
          PVR = T
          QR(I) = PVR
          QI(I) = PVI
   10 CONTINUE
      RETURN
      END
