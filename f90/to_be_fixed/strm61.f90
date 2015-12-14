SUBROUTINE strm61
     
 
!   PHASE I OF  STRESS DATA RECOVERY FOR TRIANGULAR MEMBRANE ELEMENT TRI
 
!   OUTPUTS FROM THIS PHASE FOR USE IN PHASE II ARE THE FOLLOWING
 
!   1) ELEMENT ID             WORDS    1    STORAGE IN PH1OUT   1
!   2) SIX SILS               WORDS    6                        2-7
!   3) THICKNESS T1           WORDS    1                        8
!   4) THICKNESS T2           WORDS    1                        9
!   5) THICKNESS T3           WORDS    1                       10
!   6) REFERENCE TEMP T0      WORDS    1                       11
!   7) S SUB I MATRICES       WORDS    216                     12-227
!   8) THERMAL VECTOR G ALF   WORDS    3                      228-230
 
!    EST ENTRIES SAME AS IN SUBROUTINE KTRM6S
 
 
 REAL :: nsm,ivect,jvect,kvect
 
 DIMENSION iest(45),ind(6,3),ee1(6),nph1ou(990),xc(6),yc(6),zc(6)  &
     ,   q(6,6),qq(36),ivect(3),jvect(3),kvect(3),e( 6 ),eph1(6)  &
     ,   NAME(2),ics(6),nl(6),trans(9),balotr(9),emod(9),tm(3,12) ,   tmm(36)
 
 COMMON /sdr2x5/ est(100),ph1out(250)
 COMMON /matin / matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/ em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy
 
 EQUIVALENCE (nph1ou(1),ph1out(1)),(iest(1),est(1))
 EQUIVALENCE (tm(1,1),tmm(1))
 
 DATA NAME /4HSTRM,4H61  / , BLANK  /4H    /
 DATA degra /0.0174532925/
 
 idele=iest(1)
 DO  i=1,6
   nl(i)=iest(i+1)
 END DO
 thetam=est(8)
 matid1=iest(9)
 tmem1 =est(10)
 tmem3 =est(11)
 tmem5 =est(12)
 
!   IF TMEM3 OR TMEM5 IS 0.0 OR BLANK , IT WILL BE SET EQUAL TO TMEM1
 
 IF (tmem3 == 0.0. OR .tmem3 == BLANK)  tmem3 = tmem1
 
 IF (tmem5 == 0.0. OR .tmem5 == BLANK)  tmem5 = tmem1
 
 nsm = est(13)
 
 j=0
 DO  i=14,34,4
   j=j+1
   ics(j)=iest(i)
   xc(j) = est(i+1)
   yc(j) = est(i+2)
   zc(j)=est(i+3)
 END DO
 eltemp=(est(38)+est(39)+est(40)+est(41)+est(42)+est(43))/6.0
 theta1=thetam*degra
 sinth=SIN(theta1)
 costh=COS(theta1)
 IF (ABS(sinth) <= 1.0E-06) sinth=0.0
 
 
!   EVALUATE MATERIAL PROPERTIES
 
 matflg = 2
 matid  = matid1
 CALL mat (idele)
 TO=tref
 
!   CALCULATIONS FOR THE TRIANGLE
 
 CALL trif (xc,yc,zc,ivect,jvect,kvect,dista,distb,distc,iest(1), NAME)
 
!   TRANSFORMATION MATRIX BETWEEN ELEMENT AND BASIC CO-ORDINATES
 
 e(1)=ivect(1)
 e(2)=jvect(1)
 e(3)=ivect(2)
 e(4)=jvect(2)
 e(5)=ivect(3)
 e(6)=jvect(3)
 
!   CALCULATIONS FOR  Q MATRIX AND ITS INVERSE
 
 DO  i=1,6
   DO  j=1,6
     q(i,j)=0.0
   END DO
 END DO
 DO  i=1,6
   q(i,1)=1.0
   q(i,2)=xc(i)
   q(i,3)=yc(i)
   q(i,4)=xc(i)*xc(i)
   q(i,5)=xc(i)*yc(i)
   q(i,6)=yc(i)*yc(i)
 END DO
 
!     FIND INVERSE OF Q  MATRIX
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL invers (6,q,6,qq(1),0,determ,ising,ind)
 
!   ISING EQUAL TO 2 IMPLIES THAT Q MATRIX IS SINGULAR
 
 DO  i=1,6
   DO  j=1,6
     ij=(i-1)*6+j
     qq(ij)=q(i,j)
   END DO
 END DO
 DO  i=1,9
   balotr(i)=0.0
 END DO
 
 DO  i=1,7
   ph1out(i)=est(i)
 END DO
 ph1out(8)=est(10)
 ph1out(9)=est(11)
 ph1out(10)=est(12)
 ph1out(11)=TO
 emod(1)=em(1)
 emod(2)=em(2)
 emod(3)=em(3)
 emod(4)=em(2)
 emod(5)=em(4)
 emod(6)=em(5)
 emod(7)=em(3)
 emod(8)=em(5)
 emod(9)=em(6)
 
!   STRESSES AND STRAINS ARE EVALUATED AT FOUR POINTS ,VIZ., THE THREE
!   CORNER GRID POINTS AND THE CENTROID
 
 DO  jj=1,4
   j=2*(jj-1)+1
   IF (j == 7) GO TO 103
   x=xc(j)
   y=yc(j)
   GO TO 104
   103 x=(xc(1)+xc(3)+xc(5))/3.0
   y=(yc(1)+yc(3)+yc(5))/3.0
   104 CONTINUE
   DO  i=1,36
     tmm(i)=0.0E0
   END DO
   
!   TM MATRIX IS THE PRODUCT OF B AND QINVERSE MATRICES
   
   DO  j=1,6
     j1=(j-1)*2+1
     j2=j1+1
     tm(1,j1)=q(2,j)+2.0*x*q(4,j)+y*q(5,j)
     tm(2,j2)=q(3,j)+x*q(5,j)+2.0*y*q(6,j)
     tm(3,j1)=tm(2,j2)
     tm(3,j2)=tm(1,j1)
   END DO
   DO  ii=1,6
     IF (ics(ii) == 0) GO TO 130
     CALL transs (iest(4*ii+10),trans)
     CALL gmmats (e,3,2,+1,trans,3,3,0,ee1)
     GO TO 133
     130 CONTINUE
     DO  i=1,3
       DO  j=1,2
         i1=(i-1)*2+j
         j1=(j-1)*3+i
         ee1(j1)=e(i1)
       END DO
     END DO
     133 CONTINUE
     ij1=(jj-1)*54+(ii-1)*9+12
     mz=(ii-1)*6+1
     CALL gmmats (emod,3,3,0,tmm(mz),2,3,+1,eph1)
     CALL gmmats (eph1,3,2,0,ee1,2,3,0,ph1out(ij1))
   END DO
 END DO
 CALL gmmats (emod,3,3,0,alf,3,1,0,ph1out(228))
 RETURN
END SUBROUTINE strm61
