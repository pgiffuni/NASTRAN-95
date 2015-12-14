SUBROUTINE termsd (nnode,gpth,epnorm,egpdt,iorder,mmn,bterms)
     
!     DOUBLE PRECISION ROUTINE TO CALCULATE B-MATRIX TERMS
!     FOR ELEMENTS  QUAD4, QUAD8 AND TRIA6.
 
!     THE INPUT FLAG LETS THE SUBROUTINE SWITCH BETWEEN QUAD4,
!     QUAD8 AND TRIA6 VERSIONS
 
!     ELEMENT TYPE FLAG (LTYPFL) = 1  FOR QUAD4,
!                                = 2  FOR TRIA6 (NOT AVAILABLE),
!                                = 3  FOR QUAD8 (NOT AVAILABLE).
 
!     THE OUTPUT CONSISTS OF THE DETERMINANT OF THE JACOBIAN
!     (DETJ), SHAPE FUNCTIONS AND THEIR DERIVATIVES. THE OUTPUT
!     PARAMETER, BADJAC, IS AN INTERNAL LOGICAL FLAG TO THE CALLING
!     ROUTINE INDICATING THAT THE JACOBIAN IS NOT CORRECT.
!     PART OF THE INPUT IS PASSED TO THIS SUBROUTINE THROUGH THE
!     INTERNAL COMMON BLOCK  /COMJAC/.
 
 
 INTEGER, INTENT(IN)                      :: nnode
 DOUBLE PRECISION, INTENT(IN)             :: gpth(1)
 REAL, INTENT(IN)                         :: epnorm(4,1)
 REAL, INTENT(IN)                         :: egpdt(4,1)
 INTEGER, INTENT(IN)                      :: iorder(1)
 INTEGER, INTENT(IN OUT)                  :: mmn(1)
 DOUBLE PRECISION, INTENT(OUT)            :: bterms(1)
 LOGICAL :: badjac
 INTEGER :: ltypfl, INDEX(3,3)
 
 DOUBLE PRECISION :: xi,eta,zeta,detj,shp(8),jacob(3,3),dshpx(8),  &
     dshpe(8),dshp(16),tshp(8),tdshp(16),  &
     dum,temp,eps,tie(9),tj(3,3),vn(3),cjac, th,gridc(3,8)
 COMMON /comjac/  xi,eta,zeta,detj,badjac,ltypfl
 COMMON /cjacob/  cjac(19)
 EQUIVALENCE      (dshpx(1),dshp(1)), (dshpe(1),dshp(9) )
 EQUIVALENCE      (vn(1)   ,cjac(8)), (tie(1)  ,cjac(11))
 EQUIVALENCE      (th      ,cjac(1))
 
 eps = 1.0D-15
 badjac = .false.
 
 SELECT CASE ( ltypfl )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 20
 END SELECT
 
!     QUAD4 VERSION
 
 10 ngp = 4
 CALL q4shpd (xi,eta,shp,dshp)
 GO TO 40
 
!     QUAD8 VERSION
 
 20 ngp = 8
 GO TO 40
 
!     TRIA6 VERSION
 
 30 ngp = 6
 
 40 DO  i = 1,ngp
   tshp (i  ) = shp(i)
   tdshp(i  ) = dshp(i)
   tdshp(i+8) = dshp(i+ngp)
 END DO
 DO  i = 1,ngp
   io = iorder(i)
   shp (i  ) = tshp(io)
   dshp(i  ) = tdshp(io)
   dshp(i+8) = tdshp(io+8)
 END DO
 
 th = 0.0D0
 DO  i = 1,nnode
   th = th + gpth(i)*shp(i)
   DO  j = 1,3
     j1 = j + 1
     gridc(j,i) = egpdt(j1,i) + zeta*gpth(i)*epnorm(j1,i)*0.5D0
   END DO
 END DO
 
 DO  i = 1,2
   ii = (i-1)*8
   DO  j = 1,3
     tj(i,j) = 0.0D0
     DO  k = 1,nnode
       tj(i,j) = tj(i,j) + dshp(k+ii)*gridc(j,k)
     END DO
   END DO
 END DO
 
 DO  i = 1,3
   tj(3,i) = 0.0D0
   DO  j = 1,nnode
     tj(3,i) = tj(3,i) + 0.5D0*gpth(j)*shp(j)*epnorm(i+1,j)
   END DO
 END DO
 
 DO  i = 1,3
   DO  j = 1,3
     IF (DABS(tj(i,j)) < eps) tj(i,j) = 0.0D0
   END DO
 END DO
 
!     SET UP THE TRANSFORMATION FROM THIS INTEGRATION POINT C.S.
!     TO THE ELEMENT C.S.  TIE
 
 vn(1) = tj(1,2)*tj(2,3) - tj(2,2)*tj(1,3)
 vn(2) = tj(2,1)*tj(1,3) - tj(1,1)*tj(2,3)
 vn(3) = tj(1,1)*tj(2,2) - tj(2,1)*tj(1,2)
 
 temp = DSQRT(vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3))
 
 tie(7) = vn(1)/temp
 tie(8) = vn(2)/temp
 tie(9) = vn(3)/temp
 
 temp = DSQRT(tie(8)*tie(8) + tie(9)*tie(9))
 
 tie(1) = tie(9)/temp
 tie(2) = 0.0D0
 tie(3) =-tie(7)/temp
 
 tie(4) = tie(8)*tie(3)
 tie(5) = temp
 tie(6) =-tie(1)*tie(8)
 
 CALL inverd (3,tj,3,dum,0,detj,ising,INDEX)
 
 
!     NOTE - THE INVERSE OF JACOBIAN HAS BEEN STORED IN TJ
!            UPON RETURN FROM INVERD.
 
 IF (ising == 1 .AND. detj > 0.0D0) GO TO 110
 badjac = .true.
 GO TO 150
 
 110 CONTINUE
 
 DO  i = 1,3
   ii = (i-1)*3
   DO  j = 1,3
     jacob(i,j) = 0.0D0
     DO  k = 1,3
       ik = ii + k
       jacob(i,j) = jacob(i,j) + tie(ik)*tj(k,j)
     END DO
   END DO
 END DO
 
!     MULTIPLY THE INVERSE OF THE JACOBIAN BY THE TRANSPOSE
!     OF THE ARRAY CONTAINING DERIVATIVES OF THE SHAPE FUNCTIONS
!     TO GET THE TERMS USED IN THE ASSEMBLY OF THE B MATRIX.
!     NOTE THAT THE LAST ROW CONTAINS THE SHAPE FUNCTION VALUES.
 
 node3 = nnode*3
 DO  i = 1,nnode
   bterms(node3+i) = shp(i)*jacob(3,3)
 END DO
 
 DO  i = 1,3
   ii = (i-1)*nnode
   DO  j = 1,nnode
     ij = ii + j
     bterms(ij) = 0.0D0
     DO  k = 1,2
       ik = (k-1)*8
       bterms(ij) = bterms(ij) + jacob(i,k)*dshp(ik+j)
     END DO
   END DO
 END DO
 150 RETURN
END SUBROUTINE termsd
