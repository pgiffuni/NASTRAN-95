SUBROUTINE gmmatd (a,irowa,icola,mta, b,irowb,icolb,ntb, c)
!*****
!     GMMATD - G E N E R A L  M A T R I X  M U L T I P L Y
!                                 A N D
!                           T R A N S P O S E
!            D O U B L E  P R E C I S I O N  V E R S I O N
 
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
!*****
 
 
!     IF MTA .LT. 0, C IS NOT ZEROED OUT.  HENCE THE ROUTINE, IN THIS
!     CASE, COMPUTES  A * B  +  D  =  C  WHERE THE MATRIX  D  HAS BEEN
!     STORED ROW-WISE AT  C  BY THE CALLING PROGRAM.  IF MTA = -1,  A
!     IS TRANSPOSED.  IF MTA = -2,  A  IS NOT TRANSPOSED.  NTB IS
!     DEFINED AS ABOVE AND IS INDEPENDENT OF MTA.
 
 
 
 DOUBLE PRECISION, INTENT(IN)             :: a(1)
 INTEGER, INTENT(IN)                      :: irowa
 INTEGER, INTENT(IN)                      :: icola
 INTEGER, INTENT(IN OUT)                  :: mta
 DOUBLE PRECISION, INTENT(IN)             :: b(1)
 INTEGER, INTENT(IN)                      :: irowb
 INTEGER, INTENT(IN)                      :: icolb
 INTEGER, INTENT(IN)                      :: ntb
 DOUBLE PRECISION, INTENT(OUT)            :: c(1)
 INTEGER :: rowa,cola,  rowb,colb
 
 
 
 
 
 
 
 DIMENSION     iparm(2)
 
 
 
 rowa = irowa
 cola = icola
 rowb = irowb
 colb = icolb
 nta = IABS(mta)
 IF (mta == (-2)) nta = 0
 IF (nta == 0  .AND.  ntb == 0) IF (cola - rowb) 80,5,80
 IF (nta == 1  .AND.  ntb == 0) IF (rowa - rowb) 80,5,80
 IF (nta == 0  .AND.  ntb == 1) IF (cola - colb) 80,5,80
 IF (nta == 1  .AND.  ntb == 1) IF (rowa - colb) 80,5,80
 5 IF (nta == 1) GO TO 10
 ilim= rowa
 klim= cola
 inci= cola
 incka= 1
 GO TO 20
 10 ilim= cola
 klim= rowa
 inci= 1
 incka= cola
 20 IF(ntb == 1) GO TO 30
 jlim= colb
 incj= 1
 inckb= colb
 GO TO 40
 30 jlim= rowb
 incj= colb
 inckb= 1
 40 IF (mta < 0) GO TO 47
 lim = ilim * jlim
 DO  i = 1,lim
   c(i) = 0.0D0
 END DO
 47 ij = 0
 i = 0
 50 i = i + 1
 IFIX=i*inci-cola
 j = 0
 60 j = j + 1
 ij=ij+1
 ia=IFIX
 jb=j*incj-colb
 k = 0
 70 k = k + 1
 ia=ia+incka
 jb=jb+inckb
 c(ij)=c(ij)+ a(ia) * b(jb)
 IF (k < klim) GO TO 70
 IF (j < jlim) GO TO 60
 IF (i < ilim) GO TO 50
 RETURN
 80 iparm(1) = nta
 iparm(2) = ntb
 CALL mesage (-30,21,iparm)
 RETURN
END SUBROUTINE gmmatd
