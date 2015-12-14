SUBROUTINE gmmatc( a,rowa,cola,mta, b, rowb,colb,ntb, c )
!*****
!     GMMATC - G E N E R A L  M A T R I X  M U L T I P L Y
!                                 A N D
!                           T R A N S P O S E
!            S I N G L E  P R E C I S I O N  V E R S I O N
!     COMPLEX VERSION
 
!     PERFORMS                                     WHEN
!               A            *  B            =  C     MTA=0  NTB= 0
!               A            *  B TRANSPOSE  =  C          0       1
!               A TRANSPOSE  *  B            =  C          1       0
!               A TRANSPOSE  *  B TRANSPOSE  =  C          1       1
!*****
!     A -  IS A MATRIX (ROWA) ROWS BY (COLA) COLUMNS
!     B -  IS A MATRIX (ROWB) ROWS BY (COLB) COLUMNS
!     A,B AND C ARE STORED BY ROWS (EXAMPLE)
!              MATRIX                   STORED
!         A=   1    2              A=   1
!              3    4                   2
!              5    6                   3
!                                       4
!                                       5
!                                       6
!*****
 
 
!     IF MTA .LT. 0, C IS NOT ZEROED OUT.  HENCE THE ROUTINE, IN THIS
!     CASE, COMPUTES  A * B  +  D  =  C  WHERE THE MATRIX  D  HAS BEEN
!     STORED ROW-WISE AT  C  BY THE CALLING PROGRAM.  IF MTA = -1,  A
!     IS TRANSPOSED.  IF MTA = -2,  A  IS NOT TRANSPOSED.  NTB IS
!     DEFINED AS ABOVE AND IS INDEPENDENT OF MTA.
 
 
 
 COMPLEX, INTENT(IN)                      :: a(1)
 INTEGER, INTENT(IN)                      :: rowa
 INTEGER, INTENT(IN)                      :: cola
 INTEGER, INTENT(IN)                      :: mta
 COMPLEX, INTENT(IN)                      :: b(1)
 INTEGER, INTENT(IN)                      :: rowb
 INTEGER, INTENT(IN)                      :: colb
 INTEGER, INTENT(IN)                      :: ntb
 COMPLEX, INTENT(OUT)                     :: c(1)
 
 
 
 
 INTEGER :: iparm(2)
 
 
 
 
 nta = IABS(mta)
 IF (mta == (-2)) nta = 0
 IF ( nta /= 0 ) GO TO 10
 
! A IS NOT TRANSPOSED
 
 nrowa = rowa
 ncola = cola
 incrik = 1
 ikn = cola
 incik1 = cola
 GO TO 20
 
! A IS TRANSPOSED
 
 10 nrowa = cola
 ncola = rowa
 incrik = cola
 ikn = ( rowa-1 )*cola + 1
 incik1 = 1
 20 IF( ntb /= 0 ) GO TO 30
 
! B IS NOT TRANSPOSED
 
 nrowb = rowb
 ncolb = colb
 incrkj = colb
 inckj1 = 1
 GO TO 40
 
! B IS TRANSPOSED
 
 30 nrowb = colb
 ncolb = rowb
 incrkj = 1
 inckj1 = colb
 
! CHECK CONSISTANT DIMENSIONS AND ZERO C IF NO D MATRIX
 
 40 IF( ncola /= nrowb ) GO TO 80
 IF( mta < 0 ) GO TO 50
 nterms = nrowa*ncolb
 DO  i=1,nterms
   c(i) = 0
 END DO
 
! PERFORM MATRIX MULTIPLICATION
 
 50 ij1 = 1
 ijn = ncolb
 ik1 = 1
 DO  i=1,nrowa
   kj1 = 1
   DO  ij =ij1,ijn
     kj = kj1
     DO  ik=ik1,ikn,incrik
       c(ij) = c(ij) + a(ik)*b(kj)
       kj = kj + incrkj
     END DO
     kj1 = kj1 + inckj1
   END DO
   ij1 = ijn + 1
   ijn = ijn + ncolb
   ik1 = ik1 + incik1
   ikn = ikn + incik1
 END DO
 RETURN
 80 iparm(1) = nta
 iparm(2) = ntb
 CALL mesage (-30,21,iparm(1))
 RETURN
END SUBROUTINE gmmatc
