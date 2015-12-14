SUBROUTINE flfree (frrec,afe,nafe,kge,nkge)
     
!     CALCULATES THE AREA FACTOR MATRIX AND GRAVITATIONAL STIFFNESS
!     MATRIX FOR A SINGLE FLUID ELEMENT ON THE FREE SURFACE
 
 
 INTEGER, INTENT(IN)                      :: frrec(7)
 DOUBLE PRECISION, INTENT(OUT)            :: afe(16)
 INTEGER, INTENT(OUT)                     :: nafe
 DOUBLE PRECISION, INTENT(OUT)            :: kge(16)
 INTEGER, INTENT(OUT)                     :: nkge
 LOGICAL :: error    ,grav     ,ltilt
 INTEGER :: gf1      ,gf2      ,gf3      ,iz(1)    , grid(3,4)
 DOUBLE PRECISION :: r12(3)   ,r13(3)   ,a        ,rt(3)    ,  afact    ,rhoxg
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ z(1)
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,ngbpdt   ,  &
     isil     ,nsil     ,igrav    ,ngrav
 COMMON /matin / matid    ,inflag
 COMMON /matout/ dum(3)   ,rho
 COMMON /system/ sysbuf   ,nout
 COMMON /BLANK / nograv   ,nofree   ,tilt(2)
 EQUIVALENCE     (z(1),iz(1))
 
!     GRID POINTS DEFINING FOUR OVERLAPING TRIANGLES IN A QUAD
 
 DATA    grid  / 1        ,2        ,3        ,  &
     2        ,3        ,4        , 3        ,4        ,1        ,  &
     4        ,1        ,2        /
 DATA    ltilt / .false.  /
 
 
!     CALCULATE SIZE OF ELEMENT MATRICES
 
 ngridf = 4
 IF (frrec(6) < 0) ngridf = 3
 nafe = ngridf*ngridf*2
 nkge = 0
 
!     OBTAIN MATERIAL PROPERTY AND GRAVITY DATA IF A GRAV ID IS GIVEN
 
 grav = .false.
 IF (frrec(7) == 0) GO TO 6
 inflag = 11
 matid  = frrec(2)
 CALL mat (frrec(1))
 
 IF (ngrav == 0) GO TO 70
 lgrav = igrav + ngrav - 1
 DO  i = igrav,lgrav,6
   IF (iz(i) == frrec(7)) GO TO 4
 END DO
 
 GO TO 70
 
 4 g = SQRT(z(i+3)**2 + z(i+4)**2 + z(i+5)**2)
 
!     USING THE FIRST GRAV VECTOR DETERMING THE FREE SURFACE PLOTTING
!     ANGLE
 
 IF (ltilt) GO TO 5
 tilt(1) = z(i+5)/g
 tilt(2) = z(i+3)/g
 ltilt   = .false.
 
 5 g     = g*z(i+2)
 rhoxg = DBLE(rho)*DBLE(g)
 nkge  = nafe
 nograv= 1
 grav  = .true.
 
!     DETERMINE NUMBER OF OVERLAPING TRIANGLES TO BE UESED
 
!     1 IF TRIANGLAR FLUID FACE
!     4 IF QUADRATIC FLUID FACE
 
 6 itria = 4
 IF (ngridf /= 4) itria = 1
 
!     ZERO OUT GRAVITATIONAL STIFFNESS AND AREA FACTOR MATRIX
 
 DO  i = 1,16
   kge(i) = 0.0D0
   afe(i) = 0.0D0
 END DO
 
!     LOOP OVER TRIANGLES
 
!     FIRST LOCATE GRID POINT COORDINATES FOR CORNERS FO THIS TRIANGLE
 
 DO  it = 1,itria
   
   i   = grid(1,it)
   gf1 = ibgpdt + (frrec(i+2)-1)*4
   i   = grid(2,it)
   gf2 = ibgpdt + (frrec(i+2)-1)*4
   i   = grid(3,it)
   gf3 = ibgpdt + (frrec(i+2)-1)*4
   
!     CALCUATE AREA OF TRIAGLE
!     DIVIDE AREA BY TWO IF OVERLAPPING TRIAGLES USED
   
   DO  i = 1,3
     r12(i) = z(gf2+i) - z(gf1+i)
     r13(i) = z(gf3+i) - z(gf1+i)
   END DO
   
   CALL dcross (r12,r13,rt)
   
   a = DSQRT(rt(1)*rt(1) + rt(2)*rt(2) + rt(3)*rt(3))/2.0D0
   IF (itria == 4) a = a/2.0D0
   
!     INSERT AREA AND STIFFNESS CONTRIBUTIONS INTO FULL SIZE
!     ELEMTENT MATRICES
   
   DO  i = 1,3
     icol = grid(i,it)
     iloc = ngridf*(icol-1)
     DO  j = 1,3
       irow = grid(j,it)
       IF (irow == icol) afact = a/6.0D0
       IF (irow /= icol) afact = a/12.0D0
       afe(iloc+irow) = afe(iloc+irow) + afact
       IF (grav) kge(iloc+irow) = kge(iloc+irow) + rhoxg*afact
     END DO
   END DO
   
 END DO
 
 RETURN
 
!     ERROR CONDITIONS
 
 70 WRITE  (nout,80) ufm,frrec(1),frrec(7)
 80 FORMAT (a23,' 8012, FLUID ELEMENT',i9,  &
     ' ON A CFFREE CARD REFERENCES UNDEFINED GRAVITY ID',i9)
 error = .true.
 RETURN
END SUBROUTINE flfree
