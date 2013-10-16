      SUBROUTINE RPOLYEV(NN,SR,SI,P,QR,QI,PVR,PVI)                   
C EVALUATES A REAL POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
C PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
      implicit none
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
