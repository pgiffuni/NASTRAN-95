SUBROUTINE tetra (temps,pg,iopt)
     
!     ELEMENT THERMAL LOAD GENERATOR FOR THE TETRAHEDRON SOLID ELEMENT
 
!     LOOKING DOWN ON THIS ELEMENT, GRIDS 1,2,3 ARE THE BASE AND MUST BE
!     LABELED COUNTERCLOCKWISE. GRID 4 MUST BE ABOVE THE PLANE FORMED BY
!     GRIDS 1,2,3 AND CLOSEST TO THIS OBSERVER.
 
!     ECPT FOR THE TETRAHEDRON SOLID ELEMENT
 
!     ECPT( 1) = ELEMENT ID
!     ECPT( 4) = SIL GRID POINT 3
!     ECPT( 5) = SIL GRID POINT 4
!     ECPT( 2) = MATERIAL ID (MAT1 MATERIAL TYPE)
!     ECPT( 3) = SIL GRID POINT 1
!     ECPT( 4) = SIL GRID POINT 2
!     ECPT( 5) = SIL GRID POINT 3
!     ECPT( 6) = SIL GRID POINT 4
!     ECPT( 7) = COORD SYS ID GRID PT 1
!     ECPT( 8) = X1
!     ECPT( 9) = Y1
!     ECPT(10) = Z1
!     ECPT(11) = COORD SYS ID GRID PT 2
!     ECPT(12) = X2
!     ECPT(13) = Y2
!     ECPT(14) = Z2
!     ECPT(15) = COORD SYS ID GRID PT 3
!     ECPT(16) = X3
!     ECPT(17) = Y3
!     ECPT(18) = Z3
!     ECPT(19) = COORD SYS ID GRID PT 4
!     ECPT(20) = X4
!     ECPT(21) = Y4
!     ECPT(22) = Z4
!     ECPT(23) = ELEMENT TEMPERATURE
 
 
 REAL, INTENT(IN)                         :: temps(4)
 REAL, INTENT(OUT)                        :: pg(6)
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER :: necpt(2) ,out
 REAL :: p(6)    ,c(72)    ,g(36)    ,  &
     h(16)    ,ctg(18) ,nu      ,alfa(6)  ,temp(12)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf   ,out
 COMMON /trimex/ ecpt(23)
 COMMON /matin / matid    ,inflag  ,eltemp
 COMMON /matout/ e        ,gg      ,nu      ,rho       ,alpha   ,  &
     tsub0    ,gsube   ,sigt    ,sigc      ,sigs
 EQUIVALENCE     (necpt(1),ecpt(1))
 
!     FILL THE 4 X 4 H MATRIX.
 
 h( 1) = 1.0
 h( 2) = ecpt( 8)
 h( 3) = ecpt( 9)
 h( 4) = ecpt(10)
 h( 5) = 1.0
 h( 6) = ecpt(12)
 h( 7) = ecpt(13)
 h( 8) = ecpt(14)
 h( 9) = 1.0
 h(10) = ecpt(16)
 h(11) = ecpt(17)
 h(12) = ecpt(18)
 h(13) = 1.0
 h(14) = ecpt(20)
 h(15) = ecpt(21)
 h(16) = ecpt(22)
 
!     INVERT H AND GET THE DETERMINANT
 
 ising = 0
 
 CALL invers (4,h(1),4,dum,0,hdeter,ising,temp(1))
 
!     IF THE MATRIX IS SINGULAR TETRAHEDRON IS BAD
 
 hdeter = ABS(hdeter)
 IF (ising /= 2) GO TO 200
 WRITE  (out,150) ufm,necpt(1)
 150 FORMAT (a23,' 4002, MODULE SSG1 DETECTS BAD OR REVERSE GEOMETRY ',  &
     'FOR ELEMENT ID =',i9)
 GO TO 900
 
!     GET THE MATERIAL DATA AND FILL THE 6X6 G MATERIAL STRESS-STRAIN
!     MATRIX.
 
 200 inflag = 1
 matid  = necpt(2)
 eltemp = ecpt(23)
 CALL mat (necpt(1))
 DO  i = 1,36
   g(i)  = 0.0
 END DO
 temp1 = (1.0+nu)*(1.0-2.0*nu)
 IF (temp1 /= 0.0) GO TO 240
 WRITE  (out,230) ufm,matid,ecpt(1)
 230 FORMAT (a23,' 4003, AN ILLEGAL VALUE OF -NU- HAS BEEN SPECIFIED ',  &
     'UNDER MATERIAL ID =',i9,' FOR ELEMENT ID =',i9)
 GO TO 900
 240 g( 1) = e*(1.0-nu)/temp1
 g( 8) = g(1)
 g(15) = g(1)
 g( 2) = e*nu/temp1
 g( 3) = g(2)
 g( 7) = g(2)
 g( 9) = g(2)
 g(13) = g(2)
 g(14) = g(2)
 g(22) = gg
 g(29) = gg
 g(36) = gg
 
!     FILL 4 C-MATRICES. (6X3) EACH.
 
 DO  i = 1,72
   c(i) = 0.0
 END DO
 DO  i = 1,4
   j = 18*i - 18
   c(j+ 1) = h(i+ 4)
   c(j+ 5) = h(i+ 8)
   c(j+ 9) = h(i+12)
   c(j+11) = h(i+12)
   c(j+12) = h(i+ 8)
   c(j+13) = h(i+12)
   c(j+15) = h(i+ 4)
   c(j+16) = h(i+ 8)
   c(j+17) = h(i+ 4)
 END DO
 
!     DIVIDE DETERMINANT BY 6.0, AND BY AN ADDITIONAL 2.0 IF A SUB-TETRA
!     FOR THE HEXA-10 ELEMENT.
 
 IF (iopt == 0) THEN
   GO TO   601
 ELSE
   GO TO   602
 END IF
 601 hdeter = hdeter/6.0
 GO TO 610
 602 hdeter = hdeter/12.0
 
!     INTRODUCE TBAR AND ALPHA
 
 610 hdeter = hdeter*(0.25*(temps(1)+temps(2)+temps(3)+temps(4))-tsub0) *alpha
 
!     FILL ALPHA VECTOR
 
 alfa(1) = hdeter
 alfa(2) = hdeter
 alfa(3) = hdeter
 alfa(4) = 0.0
 alfa(5) = 0.0
 alfa(6) = 0.0
 
!     LOOP FOR THE FOUR GRID POINTS
 
 DO  i = 1,4
   CALL gmmats (c(18*i-17),6,3,1, g(1),6,6,0, ctg(1))
   CALL gmmats (ctg(1),3,6,0, alfa(1),6,1,0, p(1))
   
!     TRANSFORM TO GLOBAL
   
   p(4) = 0.0
   p(5) = 0.0
   p(6) = 0.0
   k    = 4*i + 3
   IF (necpt(k) /= 0) CALL basglb (p(1),p(1),necpt(k+1),necpt(k))
   
!     INSERT LOAD VECTOR FOR GRID POINT
   
   l = necpt(i+2) - 1
   DO  j = 1,3
     l = l + 1
     pg(l) = pg(l) + p(j)
   END DO
 END DO
 RETURN
 
 900 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE tetra
