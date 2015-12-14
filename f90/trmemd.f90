SUBROUTINE trmemd
     
!     THIS SUBROUTINE CALCULATES THE STIFFNESS AND MASS MATRICES FOR
!     THE  TRIANGULAR MEMBRANE ELEMENT.  CALCULATIONS ARE PERFORMED
!     PRIMARILY BY SUBROUTINES EKTRMS AND EMASTQ.
!     DOUBLE PRECISION VERSION
 
!     ECPT FOR THE TRMEM ELEMENT
!***********************************************************************
! INDEX   DESCRIPTION                                       TYPE
! *****   ***********                                       ****
!   1     ELEMENT ID                                         I
!   2-4   GRID POINTS A,B,AND C                              I
!   5     THETA = ANGLE OF MATERIAL                          R
!   6     MATERIAL ID                                        I
!   7     T                                                  R
!   8     NON-STRUCTURAL MASS                                R
!   9     COORDINATE SYSTEM ID 1                             I
! 10-12   X1,Y1,Z1                                           R
!  13     COORDINATE SYSTEM ID 2                             I
! 14-16   X2,Y2,Z2                                           R
!  17     COORDINATE SYSTEM ID 3                             I
! 18-20   X3,Y3,Z3                                           R
!  21     ELEMENT TEMPERATURE                                R
!***********************************************************************
 DOUBLE PRECISION :: k,kout,m(9),mout(9),ksave  &
     ,                 a,prod9,temp9,xsub,bfact,e
 LOGICAL :: nogo,heat
 INTEGER :: elid,estid, dict(10), ipart(3), necpt(50), ngrid(3)
 
 COMMON /system /  ksystm (60)
 COMMON /emgprm / dm(15),ismb(3),iprec,nogo,heat,icmbar
 COMMON /emgdic /  qq(3), elid, estid
 COMMON /emgest /  ecpt(50)
 COMMON /emgtrx /  a(225),prod9(9),temp9(9),xsub(3),bfact,  &
     e(18), k(324), kout(324),ksave(81)
 
 EQUIVALENCE   (ecpt(1),necpt(1),ielid), (dict5,dict(5))
 EQUIVALENCE   (k(1),m(1)),(kout(1),mout(1)),(ksystm(2),ioutpt)
 EQUIVALENCE   (ksystm(56), iheat), (ecpt(2), ngrid(1))
 
 DATA  ipart / 1,2, 3/
 
 
 
 ip = iprec
 dict(1) = estid
 
!     CREATE AN ARRAY POINTING TO GRID POINTS IN INCREASING ORDER
 
 100 DO  i=1,2
   ip1 = i+1
   ii =  ipart(i)
   DO   j=ip1,3
     jj = ipart(j)
     IF (ngrid(ii) <= ngrid(jj)) CYCLE
     ipart(i) =jj
     ipart(j) =ii
     ii = jj
     GO TO 100
   END DO
 END DO
 
!     IF STIFFNESS MATRIX IS REQUESTED CALL EKTRMS. OTHERWISE GO TO
!     MASS MATRIX CALCULATION SECTION
 
 IF (ismb(1) == 0 ) GO TO  300
 
 CALL ektrmd (0)
 
 IF (nogo) RETURN
 
!     RE-ORDER  THE STIFFNESS MATRIX BY INCREASING SIL VALUE
 
 IF (heat) GO TO 200
 DO  i=1,3
   ii = ipart(i)
   DO  j=1,3
     jj = ipart(j)
     DO  ka=1,3
       DO  l=1,3
         isave = (ii-1)*27 + (jj-1) *9 + (ka-1)*3  + l
         iout = (i-1)*27 + (j-1)*3  +  (ka-1)*9  + l
         k(iout)=ksave(isave)
       END DO
     END DO
   END DO
 END DO
!    OUTPUT THE MATRIX
 dict(2)=1
 dict(3)=9
 dict(4)=7
 
 CALL emgout(k,k,81,1,dict,1,ip)
 GO TO 300
 
!     OUTPUT HEAT MATRIX HERE
 
 200 DO  i=1,3
   DO  j=1,3
     iout = (i-1)* 3+ j
     ik  =  (ipart(i)-1)* 3 + ipart(j)
     k(iout)=ksave(ik)
   END DO
 END DO
!     OUTPUT   HEAT  K
 dict(2) = 1
 dict(3) = 3
 dict(4) = 1
 
 CALL emgout (k,k,9,1,dict,1,ip)
 
!     PERFORM MASS MATRIX CALCULATIONS HERE
 
 300 IF (ismb(2)  == 0) RETURN
 
!     CONVENTIONAL MASS MATRIX
 
 CALL emadtq (4,m)
!     REORDER THE MASS MATRIX
 IF (heat) GO TO 350
 DO  i=1,3
   ii = (i-1)*3
   ij = ipart(i)
   jj = (ij-1)*3
   DO   j=1,3
     iout = ii + j
     ik = jj + j
     mout(iout) =  m(ik)
   END DO
 END DO
 
 dict(2) =2
 dict(3) = 9
 dict(4) = 7
 
 CALL emgout (mout, mout, 9,1,dict,2,ip)
 RETURN
 
!     HEAT FORMULATION
 
 350 DO  i=1,3
   j=ipart(i)
   mout(i)=m(j)
 END DO
 dict(2)=2
 dict(3)=3
 dict(4)=1
 
 CALL emgout(mout,mout,3,1,dict,2,ip)
 RETURN
 
END SUBROUTINE trmemd
