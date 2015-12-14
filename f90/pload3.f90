SUBROUTINE pload3
     
!     COMPUTES THE CONTRIBUTION TO THE LOAD VECTOR DUE TO PRESSURES
!     APPLIED TO THE FACES OF ISOPARAMETRIC SOLID ELEMENTS
 
 INTEGER :: gp(32)     ,seq(32)    ,face       ,sgncol(6)  ,  &
     col        ,TYPE       ,cid(32)    ,n(3)       , ibgpd(4)
 
 DOUBLE PRECISION :: shp(32)    ,dshp(3,32) ,jinv(3,3)  ,detj  &
     ,     s(3,2)     ,absisa     ,pfact      ,f(3,32)
 
 REAL :: bxyz(3,32) ,bgpd(4)    ,p(6)       ,rf(3,32)
 
 COMMON /loadx /         lcore      ,slt        ,bgpdt      , OLD
 COMMON /zzzzzz/         core(1)
 
 EQUIVALENCE (bgpd(1),ibgpd(1),dshp(1,1))
 EQUIVALENCE (seq(1),shp(1))
 EQUIVALENCE (n(1),ni)  ,(n(2),nj)  ,(n(3),nk)
 EQUIVALENCE (f(1,1),rf(1,1))
 
 DATA absisa/0.577350269189626D0/
 DATA sgncol/-3,-2,1,2,-1,3/
 
!     READ PRESSURES AND GRID POINT ID S FROM THE SLT, DETERMINE
!     ELEMENT TYPE AND NUMBER OF GRID POINTS AND GET BASIC COORDINATES.
 
 CALL READ(*500,*500,slt,p,6,0,i)
 CALL READ(*500,*500,slt,gp,32,0,i)
 TYPE=1
 ngp=8
 IF (gp(9) == 0) GO TO 10
 TYPE=2
 ngp=20
 IF (gp(21) == 0) GO TO 10
 TYPE=3
 ngp=32
 10 CALL permut(gp,seq,ngp,OLD)
 DO  i=1,ngp
   j=seq(i)
   CALL fndpnt(bgpd,gp(j))
   cid(j)=ibgpd(1)
   DO  k=1,3
     f(k,i)=0.0
     bxyz(k,j)=bgpd(k+1)
   END DO
 END DO
 
!     LOOP OVER SIX ELEMENT FACES
 
 DO  face=1,6
   IF (p(face) == 0.0) CYCLE
   j=1
   i=ISIGN(j,sgncol(face))
   sgn=FLOAT(i)
   col=IABS(sgncol(face))
   DO  i=1,3
     IF (i /= col) GO TO 40
     s(i,1)=sgn
     n(i)=1
     CYCLE
     40 s(i,1)=-absisa
     s(i,2)= absisa
     n(i)=2
   END DO
   
!     INTEGRATION LOOPS
   
   DO  i=1,ni
     DO  j=1,nj
       DO  k=1,nk
         
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE.
         
         CALL ihexsd(TYPE,shp,dshp,jinv,detj,0,s(1,i),s(2,j),s(3,k),bxyz)
         IF (detj == 0.0) CALL mesage(-61,0,0)
         pfact=detj*DBLE(sgn*p(face))
         
!     LOOP OVER GRID POINTS
         
         DO  l=1,ngp
           IF (shp(l) == 0.0) CYCLE
           DO  m=1,3
             f(m,l)=pfact*jinv(m,col)*shp(l)+f(m,l)
           END DO
         END DO
       END DO
     END DO
   END DO
 END DO
 j=3*ngp
 DO  i=1,j
   rf(i,1)=f(i,1)
 END DO
 
!     TRANSFORM VECTOR TO GLOBAL AND ADD TO GLOBAL LOAD VECTOR.
 
 DO  i=1,ngp
   IF (cid(i) == 0) GO TO 310
   CALL basglb(rf(1,i),rf(1,i),bxyz(1,i),cid(i))
   310 CALL fndsil(gp(i))
   DO  j=1,3
     k=gp(i)+j-1
     core(k)=core(k)+rf(j,i)
   END DO
 END DO
 RETURN
 500 CALL mesage(-61,0,0)
 RETURN
END SUBROUTINE pload3
