SUBROUTINE qvol
     
!     CALCULATES THERMAL LOADS DUE TO QVOL CARDS
 
 INTEGER :: ip(3),nsil(8),map(4,14),slt,reason,TYPE,bg,OLD, order(8)
 REAL :: r(4,8),data4(4,9),p(8),d12(3),d13(3),d14(3), card(12)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /condas/ consts(5)
 COMMON /loadx / lc,slt,bg,OLD,nx(12), ifm
 COMMON /zzzzzz/ core(1)
 COMMON /system/ sysbuf,iout
 EQUIVALENCE     (consts(2),twopi),(npts,card(1)),(id,card(2)),  &
     (nsil(1),card(3)),(coef,card(11)),(TYPE,card(12)),  &
     (r(1,1),data4(2,1)),(i1,ip(1)),(i2,ip(2)), (i3,ip(3))
 DATA    map   / 1  ,2  ,3  ,4  , 1  ,2  ,3  ,6  ,  &
     1  ,2  ,6  ,5  , 1  ,4  ,5  ,6  ,  &
     1  ,2  ,3  ,6  , 1  ,3  ,4  ,8  ,  &
     1  ,3  ,8  ,6  , 1  ,5  ,6  ,8  ,  &
     3  ,6  ,7  ,8  , 2  ,3  ,4  ,7  ,  &
     1  ,2  ,4  ,5  , 2  ,4  ,5  ,7  ,  &
     2  ,5  ,6  ,7  , 4  ,5  ,7  ,8  /
 
!     READ AND PROCESS ONE ELEMENT OF ONE QVOL CARD PER CALL
!     THE LOAD COEFFICIENTS ARE GENERATED AND INSERTED HERE
 
!     THE INPUT DATA ON FILE SLT IS
 
!     FIELD       DATA
!       1         NO. OF POINTS
!       2         EL. ID.
!      3-10       1 TO 8 SILS
!                        *  A*Q  FOR  TYPE=1 (RODS,ETC)
!       11        COEF = *  T*Q  FOR  TYPE=2 (TRIANGLES ETC)
!                        *    Q  FOR  TYPE=3 (BELL) OR 4 (SOLID)
!       12        TYPE
 
 CALL fread (slt,card,12,0)
 reason = 1
 IF (npts <= 1) GO TO 240
 CALL permut (nsil(1),order(1),npts,OLD)
 reason = 2
 DO  i = 1,npts
   l = order(i)
   CALL fndpnt (data4(1,l),nsil(l))
   n = nsil(l)
   CALL fndsil (n)
   IF (n /= nsil(l)) GO TO 240
   p(i) = 0.0
 END DO
 reason = 3
 IF (TYPE < 1 .OR. TYPE > 4) GO TO 240
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 40
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 120
 END SELECT
 
!     RODS, CONRODS, TUBES, BARS
 
 20 el = 0.0
 DO  i = 1,3
   el = el + (r(i,1) - r(i,2))**2
 END DO
 p(1) = coef*SQRT(el)*0.5
 p(2) = p(1)
 GO TO 200
 
!     MEMBRANES, PLATES, AND AXISYMMETRIC SOLIDS
 
 40 IF (npts == 3) GO TO 50
 IF (npts == 4) GO TO 60
 reason = 4
 GO TO 240
 50 nel  = 1
 fact =  coef/6.0
 GO TO 70
 60 nel  = 4
 fact =  coef/12.0
 70 DO  iel = 1,nel
   DO  i = 1,3
     ip(i) = i + iel - 1
     IF (ip(i) > 4) ip(i) = ip(i) - 4
   END DO
   DO  i = 1,3
     d12(i) = r(i,i2) - r(i,i1)
     d13(i) = r(i,i3) - r(i,i1)
   END DO
   CALL saxb (d12(1),d13(1),d12(1))
   el =  fact*SQRT(d12(1)**2 + d12(2)**2 + d12(3)**2)
   IF (TYPE == 2) GO TO 100
   
!     SPECIAL FACTOR FOR  AXISYMMETRIC ELEMENTS
   
   el    = el*twopi*(r(1,i1) + r(1,i2) + r(1,i3))/3.0
   100 p(i1) = p(i1) + el
   p(i2) = p(i2) + el
   p(i3) = p(i3) + el
 END DO
 GO TO 200
 
!     SOLID ELEMENTS
 
 120 IF (npts == 4) GO TO 130
 IF (npts == 6) GO TO 140
 IF (npts == 8) GO TO 150
 reason = 5
 GO TO 240
 130 nel  = 1
 fact = coef/24.0
 imap = 1
 GO TO 160
 140 nel  = 3
 imap = 2
 fact = coef/24.0
 GO TO 160
 150 imap = 5
 nel  = 10
 fact = coef/48.0
 160 DO  iel = 1,nel
   im = imap + iel - 1
   i1 =  map(1,im)
   i2 =  map(2,im)
   i3 =  map(3,im)
   i4 =  map(4,im)
   DO   i = 1,3
     
     d12(i) = r(i,i2) - r(i,i1)
     d13(i) = r(i,i3) - r(i,i1)
     d14(i) = r(i,i4) - r(i,i1)
   END DO
   
   CALL saxb (d12(1),d13(1),d12(1))
   el = fact*ABS(d12(1)*d14(1) + d12(2)*d14(2) + d12(3)*d14(3))
   DO  i = 1,4
     l = map(i,im)
     p(l) = p(l) + el
   END DO
 END DO
 
!     INSERT THE LOADS
 
 200 DO  i = 1,npts
   isil = nsil(i)
   core(isil) = core(isil) + p(i)
 END DO
 RETURN
 
!     ERROR MESSAGE
 
 240 WRITE  (iout,250) sfm,id,reason
 250 FORMAT (a25,' 3093, ELEMENT =',i9,'.   REASON =',i7)
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE qvol
