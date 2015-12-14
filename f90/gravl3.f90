SUBROUTINE gravl3 (nvect,gvect,sr1,iharm)
     
!     BUILD GRAVITY LOADS FOR AXISYMMETRIC SHELL
 
!     DEFINITION OF VARIABLES
 
!     NVECT    NUMBER OF GRAVITY LOADS
!     GVECT    ARRAY OF G VECTORS
!     SR1      FILE TO PUT ACCELERATION VECTOR ON
!     IHARM    SINE OR COSINE SET FLAG -- 1 = SINE SET
!     LUSET    LENGTH OF G SET
!     MCB      MATRIX CONTROL BLOCK FOR SR1
!     M        NUMBER OF RINGS
!     N        NUMBER OF HARMONICS
!     IL       POINTER IN GVECT ARRAY
 
 
 INTEGER, INTENT(IN)                      :: nvect
 REAL, INTENT(IN)                         :: gvect(1)
 INTEGER, INTENT(IN OUT)                  :: sr1
 INTEGER, INTENT(IN OUT)                  :: iharm
 EXTERNAL        rshift,andf
 INTEGER :: andf,rshift,sysbuf,mcb(7)
 
 DIMENSION       isystm(175)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /BLANK / luset
 COMMON /system/ sysbuf,ix(25),mn
 COMMON /zzzzzz/ z(1)
 COMMON /zblpkx/ b(4),ii
 EQUIVALENCE (sysbuf, isystm(1))
 
!     INITIALIZE STUFF
 
 ibuf = korsz(z) - sysbuf + 1
 CALL gopen (sr1,z(ibuf),1)
 CALL makmcb (mcb,sr1,luset,2,1)
 il    = 1
 n = mn
 m = isystm(161)
 
!     BUILD NVECT GRAVITY VECTORS
 
 DO  iloop = 1,nvect
   CALL bldpk (1,1,mcb(1),0,0)
   
!     COMPUTE VALUES
   
   sinth = 0.0
   sinph = 0.0
   cosph = 1.0
   g = SQRT(gvect(il)*gvect(il)+gvect(il+1)*gvect(il+1)+gvect(il+2)*  &
       gvect(il+2))
   costh = gvect(il+2)/g
   IF (gvect(il) == 0.0 .AND. gvect(il+1) == 0.0) GO TO 30
   gxy   = SQRT(gvect(il)*gvect(il)+ gvect(il+1)*gvect(il+1))
   sinth = gxy/g
   sinph = gvect(il+1)/gxy
   cosph = gvect(il  )/gxy
   30 CONTINUE
   SELECT CASE ( iharm )
     CASE (    1)
       GO TO 40
     CASE (    2)
       GO TO 50
   END SELECT
   
!     SINE SET
   
   40 b(1) = g*sinth*sinph
   ii = luset - m*(n-1)*6 + 1
   DO  i = 1,m
     CALL zblpki
     ii = ii +1
     CALL zblpki
     ii = ii +5
   END DO
   GO TO 110
   
!     COSINE SET
   
   50 b(1)=  g*costh
   ii  = luset - m*n*6 + 3
   
!     LOAD ZERO HARMONIC
   
   DO  i = 1,m
     CALL zblpki
     ii = ii + 6
   END DO
   
!     LOAD 2-D HARMONIC
   
   ii = ii - 2
   b(1) = g*sinth*cosph
   DO  i = 1,m
     CALL zblpki
     ii = ii +1
     b(1) = -b(1)
     CALL zblpki
     b(1) = -b(1)
     ii = ii +5
   END DO
   
!     END OF COLUMN
   
   110 CALL bldpkn (mcb(1),0,mcb(1))
   il = il + 3
 END DO
 CALL CLOSE (mcb(1),1)
 CALL wrttrl (mcb)
 RETURN
 
END SUBROUTINE gravl3
