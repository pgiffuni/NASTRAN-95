SUBROUTINE trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,id,elem)
     
!     CALCULATEIONS FOR THE TRIANGLE USED IN TRIM6,TRPLT1,TRSHL - THE HI
!     LEVEL PLATE ELEMENTS.  COMPUTATIONS IN SINGLE PRECISION ONLY
 
!     IVECT, JVECT, AND KVECT ARE UNIT VECTORS OF THE TRIANGLE
!     B IS THE DISTANCE OF THE GRID POINT 1
!     A IS THE DISTANCE OF THE GRID POINT 3
!     C IS THE DISTANCE OF THE GRID POINT 5
 
 
 REAL, INTENT(IN OUT)                     :: xc(6)
 REAL, INTENT(IN OUT)                     :: yc(6)
 REAL, INTENT(IN)                         :: zc(6)
 REAL, INTENT(OUT)                        :: ivect(3)
 REAL, INTENT(OUT)                        :: jvect(3)
 REAL, INTENT(OUT)                        :: kvect(3)
 REAL, INTENT(OUT)                        :: a
 REAL, INTENT(OUT)                        :: b
 REAL, INTENT(OUT)                        :: c
 INTEGER, INTENT(IN OUT)                  :: id
 REAL, INTENT(IN OUT)                     :: elem(2)
 LOGICAL :: nogo
 
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,nout,nogo
 
!     EVALUATE DIRECTIONAL COSINES
 
 x1 = xc(3) - xc(1)
 y1 = yc(3) - yc(1)
 z1 = zc(3) - zc(1)
 x2 = xc(5) - xc(1)
 y2 = yc(5) - yc(1)
 z2 = zc(5) - zc(1)
 temp = x1*x1 + y1*y1 + z1*z1
 IF (temp <= 1.0E-10) GO TO 40
 temp = SQRT(temp)
 
!     I-VECTOR
 
 ivect(1) = x1/temp
 ivect(2) = y1/temp
 ivect(3) = z1/temp
 SAVE = temp
 
!     NON-NORMALIZED K-VECTOR
 
 kvect(1) = ivect(2)*z2 - y2*ivect(3)
 kvect(2) = ivect(3)*x2 - z2*ivect(1)
 kvect(3) = ivect(1)*y2 - x2*ivect(2)
 temp = SQRT(kvect(1)**2 + kvect(2)**2 + kvect(3)**2)
 IF (temp <= 1.0E-10) GO TO 50
 
!     NORMALIZE K-VECTOR
!     DISTANCE C OF THE TRAINGLE IS TEMP
 
 kvect(1) = kvect(1)/temp
 kvect(2) = kvect(2)/temp
 kvect(3) = kvect(3)/temp
 c = temp
 
!     J-VECTOR = K X I VECTORS
 
 jvect(1) = kvect(2)*ivect(3) - ivect(2)*kvect(3)
 jvect(2) = kvect(3)*ivect(1) - ivect(3)*kvect(1)
 jvect(3) = kvect(1)*ivect(2) - ivect(1)*kvect(2)
 temp = SQRT(jvect(1)**2 + jvect(2)**2 + jvect(3)**2)
 IF (temp <= 1.0E-10) GO TO 60
 
!     NORMALIZE J-VECTOR TO MAKE SURE
 
 jvect(1) = jvect(1)/temp
 jvect(2) = jvect(2)/temp
 jvect(3) = jvect(3)/temp
 
!     DISTANCE B OF THE TRIANGLE IS OBTAINED BY DOTTING (X2,Y2,Z2) WITH
!     THE IVECT UNIT VECTOR
 
 b = x2*ivect(1) + y2*ivect(2) + z2*ivect(3)
 
!     THE LOCAL X AND Y COORINATES OF THE SIX GRID PTS. ARE AS FOLLOWS
 
 yc(1) = 0.0
 yc(2) = 0.0
 yc(3) = 0.0
 yc(4) = c*0.5
 yc(5) = c
 yc(6) = yc(4)
 
!     THE TRIANGLE SHOULD BELONG TO
 
!     KASE1 (ACUTE ANGLES AT GRID POINTS 1 AND 3),
!     KASE2 (OBTUSE ANGLE AT GRID POINT 3), OR
!     KASE3 (OBTUSE ANGLE AT GRID POINT 1)
 
!     KASE  = 1
!     IF (B .GT. SAVE) KASE = 2
!     IF (B .LT.  0.0) KASE = 3
 temp  = -b
!     IF (KASE .EQ. 3) TEMP = ABS(B)
!     IF (B .LT.  0.0) TEMP = ABS(B)
 xc(1) = temp
 xc(2) = temp + SAVE*0.5
 xc(3) = temp + SAVE
 xc(4) = xc(3)*0.5
 xc(5) = 0.0
 xc(6) = xc(1)*0.5
 
!     RE-SET DISTANCE A AND B
 
 b = ABS(b)
 a = ABS(xc(3))
 RETURN
 
!     GEOMETRY ERRORS
 
 40   WRITE (nout,140) ufm,elem,id
 GO TO 80
 50   WRITE (nout,150) ufm,elem,id
 GO TO 80
 60   WRITE (nout,160) ufm,elem,id
 80   nogo = .true.
 
 140  FORMAT (a23,' 2404, GRID POINTS 1 AND 3 OF ',a4,a2,  &
     ' WITH ELEMENT ID =',i9,' HAVE SAME COORDINATES.')
 150  FORMAT (a23,' 2405, GRID POINTS 1, 3, AND 5 OF ',a4,a2,' WITH ',  &
     'ELEMENT ID =',i9,' APPEAR TO BE ON A STRAIGHT LINE.')
 160  FORMAT (a23,' 2406, GRID POINTS 1 AND 5 OF ',a4,a2,  &
     ' WITH ELEMENT ID =',i9,' HAVE SAME COORDINATES.')
 RETURN
END SUBROUTINE trif
