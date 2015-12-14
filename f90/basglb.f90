SUBROUTINE basglb (vin1,vout1,pont,icstm)
     
    !     THIS ROUTINE CONTAINS FOUR ENTRY POINTS
 
    !     1- BASGLB TRANSFORMS A VECTOR FROM BASIC TO GLOBAL
    !     2- GLBBAS TRANSFORMS A VECTOR FROM GLOBAL TO BASIC
    !     3- FDCSTM FINDS THE LOGICAL RECORD ON THE CSTM FOR A PARTICULAR ID
    !     4- GBTRAN FINDS A PARTICULAR GLOBAL TO BASIC TRANSFORMATION AND
    !        RETURNS IT AS A 3 X 3 STORED BY ROWS.
 
 
 
    REAL, INTENT(IN)                         :: vin1(3)
    REAL, INTENT(IN OUT)                     :: vout1(3)
    REAL, INTENT(IN)                         :: pont(3)
    INTEGER, INTENT(IN)                      :: icstm
    LOGICAL :: tonly
    INTEGER :: cstm,tysys,check
    REAL :: t(9)
    DIMENSION       vin(3), ti(3,3),tl(3,3),  pont1(3),tz(3,3),iparm(2)
    COMMON /xcstm / tz
    COMMON /loadx / lc(4),cstm,lc1(10),idum(3),icm
    COMMON /tranx / nsys,tysys,ro(3),TO(3,3)
    COMMON /system/ ibuf,nout
    DATA    iparm/ 4HBASG,2HLB /
 
    !     NSYS  IS SYSTEM NUMBER
    !     TYSYS IS SYSTEM TYPE
    !     RO IS LOCATION OF ORIGIN
    !     TO IS ROTATION MATRIX
 
    tonly = .false.
    check = 123456789
    ASSIGN 90 TO iexit
    GO TO 10
 
 
    ENTRY gbtran (icstm,pont,t)
    !     ===========================
 
    IF (icstm == 0) GO TO 300
    IF (tysys >= 2 .AND. check /= 123456789) WRITE (nout,5)
5   FORMAT ('0*** SYSTEM POTENTIAL ERROR, GBTRAN WAS CALLED WITHOUT',  &
        ' FIRST CALLING BASGLB')
    check = 123456789
    tonly = .true.
    GO TO 235
 
 
    ENTRY fdcstm (icstm)
    !     ====================
 
    tonly = .false.
    ASSIGN 50 TO iexit
 
!     FDCSTM WILL FIND REQUESTED SYSTEM (ICSTM)
 
10 CONTINUE
   IF (icstm == 0) GO TO 81
   IF (icm   /= 0) GO TO 80
   IF (icstm-nsys == 0) THEN
       GO TO    40
   END IF
20 CALL READ (*60,*80,cstm,nsys,14,0,flag)
   IF (icstm-nsys == 0) THEN
       GO TO    30
   ELSE
       GO TO    20
   END IF
30 CALL bckrec (cstm)
40 GO TO iexit, (90,240,50)
50 RETURN
 
60 n1 = -2
   iparm1 = cstm
70 CALL mesage (n1,iparm1,iparm)
 
   !     UNABLE TO FIND REQUESTED COORDINATE SYSTEM
 
80 n1      =-30
   iparm1  = 25
   iparm(1)= icstm
   GO TO 70
 
   !     REQUEST FOR BASIC COORDINATE SYSTEM
 
81 CONTINUE
   tysys = 1
   nsys  = 0
   ro(1) = 0.0
   ro(2) = 0.0
   ro(3) = 0.0
   DO  i = 1,3
       DO  j = 1,3
           TO(j,i) = 0.0
       END DO
   END DO
   TO(1,1) = 1.0
   TO(2,2) = 1.0
   TO(3,3) = 1.0
   GO TO 40
 
   !     CONVERTS BASIC TO GLOBAL
 
90 ioth = 0
 
   !     RECTANGULAR
 
   100 DO  i = 1,3
       DO  j = 1,3
           tz(i,j) = TO(j,i)
           ti(i,j) = TO(j,i)
       END DO
       vin(i)  = vin1(i)
   END DO
   IF (tysys-2 < 0.0) THEN
       GO TO   130
   ELSE
       GO TO   140
   END IF
130 CALL mpyl (ti(1,1),vin(1),3,3,1,vout1(1))
   GO TO 50
 
   !     CYLINDRICAL
 
   140 DO  i = 1,3
       pont1(i) = pont(i) - ro(i)
   END DO
   CALL mpyl (ti(1,1),pont1(1),3,3,1,vin(1))
   DO  i = 1,3
       DO  j = 1,3
           tl(i,j) = 0.0
       END DO
   END DO
   r = SQRT(vin(1)*vin(1) + vin(2)*vin(2))
   IF (r  ==  0.0) GO TO 210
   IF (tysys > 2) GO TO 230
   tl(3,3) = 1.0
   tl(1,1) = vin(1)/r
   tl(2,2) = tl(1,1)
   tl(2,1) = vin(2)/r
   tl(1,2) =-tl(2,1)
170 CALL mpyl (tl(1,1),ti(1,1),3,3,3,tz(1,1))
180 IF (tonly) GO TO 201
   IF ( ioth == 0) THEN
       GO TO   190
   ELSE
       GO TO   270
   END IF
   190 DO  i = 1,3
       vin(i) = vin1(i)
   END DO
   CALL mpyl (tz(1,1),vin(1),3,3,1,vout1(1))
   GO TO 50
 
   !     RETURN THE TRANSFORMATION ONLY
 
201 t(1) = tz(1,1)
   t(2) = tz(1,2)
   t(3) = tz(1,3)
   t(4) = tz(2,1)
   t(5) = tz(2,2)
   t(6) = tz(2,3)
   t(7) = tz(3,1)
   t(8) = tz(3,2)
   t(9) = tz(3,3)
   GO TO 50
 
   !     ORIENTATION ARBITARY   TL = I   I.E. TZ = TI
 
   210 DO  i = 1,3
       DO  j = 1,3
           tz(i,j) = ti(i,j)
       END DO
   END DO
   GO TO 180
 
   !     SPHERICAL
 
230 xl = SQRT(vin(1)*vin(1) + vin(2)*vin(2) + vin(3)*vin(3))
   xr = vin(1)/r
   yr = vin(2)/r
   zl = vin(3)/xl
 
   !     BUILD TL TRANSPOSE
 
   tl(1,1) = vin(1)/xl
   tl(1,2) = xr*zl
   tl(1,3) =-yr
   tl(2,1) = vin(2)/xl
   tl(2,2) = yr*zl
   tl(2,3) = xr
   tl(3,1) = zl
   tl(3,2) =-r/xl
   GO TO 170
 
 
   ENTRY glbbas (vin1,vout1,pont,icstm)
   !     ====================================
 
   tonly = .false.
235 ASSIGN 240 TO iexit
   ioth = 1
   GO TO 10
 
   !     CONVERTS FROM GLOBAL TO BASIC
 
240 IF (tysys-2 < 0.0) THEN
       GO TO   250
   ELSE
       GO TO   100
   END IF
250 IF ( tonly ) GO TO 261
   DO  i = 1,3
       vin(i) = vin1(i)
   END DO
   CALL mpyl (TO(1,1),vin(1),3,3,1,vout1(1))
   GO TO 50
 
   !     RETURN THE TRANSFORMATION ONLY.
 
261 t(1) = TO(1,1)
   t(2) = TO(2,1)
   t(3) = TO(3,1)
   t(4) = TO(1,2)
   t(5) = TO(2,2)
   t(6) = TO(3,2)
   t(7) = TO(1,3)
   t(8) = TO(2,3)
   t(9) = TO(3,3)
   GO TO 50
 
   !     COMPUTE TL TRANSPOSE
 
   !     TRANSPOSE ROTATION PRODUCT
 
   270 DO  i = 1,3
       vin(i) = vin1(i)
       DO  j = 1,3
           ti(i,j) = tz(j,i)
       END DO
   END DO
   CALL mpyl (ti(1,1),vin(1),3,3,1,vout1(1))
   GO TO 50
 
   !     COORDINATE SYSTEM 0
 
   300 DO  i = 2,8
       t(i) = 0.
   END DO
   t(1) = 1.
   t(5) = 1.
   t(9) = 1.
   GO TO 50
   END SUBROUTINE basglb
