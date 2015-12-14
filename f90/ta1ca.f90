SUBROUTINE ta1ca(koz)
!*****
! THIS ROUTINE, CALLED BY SUBROUTINE TA1C, COMPUTES THE S MATRIX OF A
! GENERAL ELEMENT FROM INFORMATION IN THE CSTM AND BGPDT DATA BLOCKS.
! SEE FMMS-57 FOR EQUATIONS.
!*****
 
 INTEGER, INTENT(IN OUT)                  :: koz
 DOUBLE PRECISION :: v(3)               ,t(9)  &
     ,                  e(18)              ,d(42)  &
     ,                  s(6)               ,det  &
     ,                  b(6)               ,INDEX(18)  &
     ,                  dd(30)             ,dl(25) ,                  du(25)
 
 
 
 INTEGER :: cstm               ,bgpdt  &
     ,                  gei                ,FILE  &
     ,                  clsrw              ,eor  &
     ,                  bufr1              ,bufr2 ,                  bufr3
 
 
 
 DIMENSION NAME(2)            ,ssp(6)  &
     ,                  z(1) ,                  lrow(5)            ,icol(6)
 
 
 
 COMMON   /ta1com/ dum3(3)            ,bgpdt  &
     ,                  dum2(2)            ,cstm  &
     ,                  dum22(2)           ,gei
 
! OPEN CORE
 
 COMMON   /zzzzzz/ iz(1)
 
 
 
 COMMON   /names / dummy1             ,inrw  &
     ,                  dummy2             ,outrw ,                  clsrw
 
 
 
 COMMON   /tac1ax/ bufr1              ,bufr2  &
     ,                  bufr3              ,iui  &
     ,                  nui                ,iud  &
     ,                  nud                ,izzz  &
     ,                  nogo               ,idgenl
 
 
 
 EQUIVALENCE (z(1),iz(1))
 
 
 
 DATA     eor,neor /1,0/
 DATA     NAME(1)/4HTA1C/ , NAME(2)/4HA   /
 
! INITIALIZE
 
 ncstm = 0
 icstm = izzz
 left  = bufr3 - icstm
 
! ATTEMPT TO OPEN THE CSTM
 
 FILE = cstm
 CALL OPEN(*20,cstm,z(bufr3),inrw)
 CALL fwdrec(*9020,cstm)
 CALL READ(*9020,*10,cstm,z(icstm+1),left,eor,ncstm)
 CALL mesage (-8,0,NAME(1))
 10 CALL CLOSE (cstm,clsrw)
 
! PRETRD SETS UP SUBSEQUENT CALLS TO TRANSD
 
 CALL pretrd (z(icstm+1),ncstm)
 left = left - ncstm
 
! READ THE BGPDT INTO CORE
 
 20 ibgpdt = icstm + ncstm
 FILE = bgpdt
 CALL OPEN(*9010,bgpdt,z(bufr3),inrw)
 CALL fwdrec(*9020,bgpdt)
 CALL READ(*9020,*30,bgpdt,z(ibgpdt+1),left,eor,nbgpdt)
 CALL mesage (-8,0,NAME(1))
 30 CALL CLOSE (bgpdt,clsrw)
 
! ZERO OUT THE E MATRIX
 
 DO  i = 1,18
   e(i) = 0.0D0
 END DO
 e(1)  = 1.0D0
 e(8)  = 1.0D0
 e(15) = 1.0D0
 ind = 0
 50 ind = ind + 1
!*****
! IF IND = 1, THE D MATRIX IS FORMED IN THE DO 200 LOOP.
! IF IND = 2, THE S MATRIX IS FORMED AND OUTPUT A ROW AT A TIME IN THE
! DO LOOP.
!*****
 IF (ind - 2 < 0) THEN
   GO TO    60
 ELSE IF (ind - 2 == 0) THEN
   GO TO    70
 ELSE
   GO TO   300
 END IF
 60 CONTINUE
 
!     IF STIFFNESS IS INPUT,CALCULATE LIM
 
 IF (koz == 1) GO TO 65
 lim = 6
 lima = 6
 ibeg = iud
 GO TO 80
 65 lim = (nud - iud) / 4 + 1
 lima = lim
 ibeg = iud
 GO TO 80
 70 lim = (nui - iui) / 4  +  1
 ibeg = iui
 irow = 37
 80 j = ibeg - 2
 i = 1
 85 CONTINUE
 IF (ind == 1) irow = 6*i - 5
 j = j + 4
 jj = iz(j+1)
 k = ibgpdt + 4*(iz(j) - 1)
 
! COMPUTE THE V VECTOR
 
 v(1) = 0.0D0
 v(2) = 0.0D0
 v(3) = 0.0D0
 kk = jj
 IF (jj > 3) kk = jj - 3
 IF (iz(k+1) == 0) GO TO 120
 CALL transd (iz(k+1),t)
 SELECT CASE ( kk )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 110
 END SELECT
 90 v(1) = t(1)
 v(2) = t(4)
 v(3) = t(7)
 GO TO 130
 100 v(1) = t(2)
 v(2) = t(5)
 v(3) = t(8)
 GO TO 130
 110 v(1) = t(3)
 v(2) = t(6)
 v(3) = t(9)
 GO TO 130
 120 v(kk) = 1.0D0
 
! FORM THE E MATRIX IF THE DEGREE OF FREEDOM IS A TRANSLATION.
 
 130 IF (jj > 3) GO TO 150
 e( 5) =  z(k+4)
 e( 6) = -z(k+3)
 e(10) = -z(k+4)
 e(12) =  z(k+2)
 e(16) =  z(k+3)
 e(17) = -z(k+2)
 IF (iz(k+1) == 0) GO TO 140
 CALL gmmatd (v,3,1,1, e,3,6,0, d(irow) )
 GO TO 180
 140 ierow = 6*jj - 5
 d(irow  ) = e(ierow  )
 d(irow+1) = e(ierow+1)
 d(irow+2) = e(ierow+2)
 d(irow+3) = e(ierow+3)
 d(irow+4) = e(ierow+4)
 d(irow+5) = e(ierow+5)
 GO TO 180
 
! THE DEGREE OF FREEDOM IS A ROTATION.
 
 150 ll = irow
 DO  l = 1,6
   d(ll) = 0.0D0
   ll = ll + 1
 END DO
 IF (iz(k+1) == 0) GO TO 170
 d(irow+3) = v(1)
 d(irow+4) = v(2)
 d(irow+5) = v(3)
 GO TO 180
 170 ll = irow + jj - 1
 d(ll) = 1.0D0
 
! IF IND = 2 FORM A ROW OF THE S MATRIX AND WRITE IT OUT.
 
 180 IF (ind == 1) GO TO 200
 
!     IF STIFFNESS MATRIX INPUT AND LESS THAN 6 RIGID BODY DEGREES OF
!     FREEDOM, BRANCH
 
 IF (koz == 1.AND.lima < 6) GO TO 410
 CALL gmmatd (d(37),6,1,1, d(1),6,6,0, s(1) )
 DO  l = 1,6
   ssp(l) = s(l)
 END DO
 CALL WRITE (gei,ssp,6,neor)
 195 CONTINUE
 200 CONTINUE
 i = i + 1
 IF (i <= lim) GO TO 85
 IF (ind /= 1) GO TO 300
 
!     IF STIFFNESS MATRIX WAS INPUT AND LESS THAN 6 RIGID BODY DEGREES
!     OF FREEDOM, BRANCH
 
 IF (koz == 1.AND.lim < 6) GO TO 310
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL inverd (6,d(1),6,b(1),0,det,ising,INDEX(1))
 IF (ising == 1) GO TO 50
 nogo = 1
 CALL mesage (30,82,idgenl)
 GO TO 300
 310 CONTINUE
 
!     SWITCH FROM ROW STORED TO COLUMN STORED
 
 DO  i=1,lim
   DO  j=1,6
     indxz = i+(j-1)*lim
     dd(indxz) = d(6*i+j-6)
   END DO
 END DO
 
!     DETERMINE RANK OF DD AND EXPRESS MATRIX OF MAXIMAL RANK AS A
!     PRODUCT OF TRIANGULAR FACTORS
 
 CALL dmfgr(dd,lim,6,1.05E-05,irank,lrow,icol)
 IF (irank == lim) GO TO 325
 nogo = 1
 CALL mesage (30,152,idgenl)
 GO TO 300
 
!     EXTRACT LOWER AND UPPER TRIANGULAR FACTORS, FORM PRODUCT,RESTORE
!     ROWS TO THEIR POSITION BEFORE FACTORIZATION AND INVERT. THEN
!     EXPAND MATRIX TO BE OF DIMENSION  6 BY IRANK
 
 325 IF (irank == 1) GO TO 365
 DO  i=1,25
   dl(i) = 0.0D0
   du(i) = 0.0D0
 END DO
 DO  i=1,irank
   DO  j=1,irank
     IF (i > j) GO TO 340
     IF (i == j) GO TO 330
     indxz = i+(j-1)*irank
     dl(indxz) = dd(i*irank+j-irank)
     CYCLE
     330 indxz = i+(j-1)*irank
     dl(indxz) = 1.0D0
     340 indxz = i+(j-1)*irank
     du(indxz) = dd(i*irank+j-irank)
   END DO
 END DO
 CALL gmmatd (dl(1),lim,lim,0,du(1),lim,lim,0,dd)
 DO  i=1,lim
   k = lrow(i)
   DO  j=1,lim
     indxz = j+(k-1)*lim
     d(indxz) = dd(i*lim+j-lim)
   END DO
 END DO
!     AGAIN NO NEED TO COMPUTE DETERMINANT
 ising = -1
 CALL inverd (lim,d(1),lim,b(1),0,det,ising,INDEX(1))
 IF (ising == 1) GO TO 370
 nogo = 1
 CALL mesage (30,153,idgenl)
 GO TO 300
 365 d(1) = 1.0D0/dd(1)
 370 k = lim * lim + 1
 j = lim * 6
 DO   i = k,j
   d(i) = 0.0D0
 END DO
 GO TO 50
 410 CONTINUE
 
!     REARRANGE COLUMNS TO AGREE WITH ORDER OF DD AFTER MATRIX FACTOR-
!     IZATION
 
 DO  l = 1,6
   lk = icol(l)
   b(l) = d(36+lk)
 END DO
 
!     MULTIPLY DI BY THE EXPANDED INVERSE OF DD
 
 CALL gmmatd (b(1),1,6,0,d(1),6,lima,0,s(1))
 
!     WRITE OUT THIS ROW OF THE S MATRIX
 
 DO  l =1,lima
   ssp(l) = s(l)
 END DO
 CALL WRITE (gei,ssp,lima,neor)
 GO TO 195
 300 RETURN
 
!     ERROR MESSAGES
 
 9010 CALL mesage (-1,FILE,NAME(1))
 9020 CALL mesage (-1,FILE,NAME(1))
 RETURN
END SUBROUTINE ta1ca
