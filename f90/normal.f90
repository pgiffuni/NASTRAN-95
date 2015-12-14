SUBROUTINE normal
     
!     THIS IS THE DRIVER FOR THE NORM MODULE.
 
!     NORM        INMAT/OUTMAT/S,N,NCOL/S,N,NROW/S,N,XNORM/V,Y,IOPT $
 
!     DEPENDING ON THE VALUE OF IOPT, THIS MODULE PERFORMS THE
!     FOLLOWING FUNCTIONS --
 
!     IOPT = 'MAX'
!                 NORM GENERATES A MATRIX.  EACH COLUMN OF THIS OUTPUT
!                 MATRIX REPRESENTS A COLUMN OF THE INPUT MATRIX
!                 NORMALIZED BY ITS LARGEST ROW ELEMENT. (DEFAULT)
 
!     IOPT = 'SRSS'
!                 NORM GENERATES A COLUMN VECTOR.  EACH ELEMENT OF THIS
!                 VECTOR REPRESENTS THE SQUARE ROOT OF THE SUM OF THE
!                 SQUARES (SRSS) OF THE CORRESPONDING ROW OF THE INPUT
!                 MATRIX.
 
 
 
!     INPUT DATA BLOCK --
 
!     INMAT     - ANY MATRIX
 
!     OUTPUT DATA BLOCK --
 
!     OUTMAT    - OUTPUT MATRIX GENERATED AS DESCRIBED BELOW
 
!     PARAMETERS --
 
!     NCOL      - NO. OF COLUMNS OF THE INPUT MATRIX (OUTPUT/INTEGER)
 
!     NROW      - NO. OF ROWS OF THE INPUT MATRIX (OUTPUT/INTEGER)
 
!     XNORM     - MAX. NORMALIZING OR SRSS VALUE, DEPENDING UPON THE
!                 IOPT VALUE SPECIFIED (OUTPUT/REAL)
!     IOPT      - OPTION INDICATING WHETHER EACH COLUMN OF THE INPUT
!                 MATRIX IS TO BE NORMALIZED BY THE MAXIMUM ROW ELEMENT
!                 IN THAT COLUMN OR WHETHER THE SRSS VALUE FOR EACH ROW
!                 OF THE INPUT MATRIX IS TO BE COMPUTED (INPUT/BCD)
 
!     THIS MODULE DEVELOPED BY P. R. PAMIDI OF RPK CORPORATION,
!     MARCH 1988
 
 DIMENSION        mcb(7), z(1)   , isubnm(2)
 DOUBLE PRECISION :: dxmax , zd(1)  , dzero
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  ncol   , nrow   , xxmax , iopt(2)
 COMMON /packx /  ipkot1 , ipkot2 , ip1   , ip2   , incrp
 COMMON /system/  isysbf , nout
 COMMON /TYPE  /  iprc(2), nwds(4), irc(4)
 COMMON /unpakx/  iunout , iu1    , iu2   , incru
 COMMON /zzzzzz/  iz(1)
 EQUIVALENCE      (iz(1),z(1),zd(1)), (iopt1,iopt(1))
 DATA    matin ,  matout / 101, 201/
 DATA    isubnm,           MAX     , isrss , iblnk  , dzero  /  &
     4HNORM,  4HAL   , 4HMAX   , 4HSRSS, 4H     , 0.0D+0 /
 
 IF (iopt(2) == iblnk .AND. (iopt1 == MAX .OR. iopt1 == isrss)) GO TO 20
 WRITE  (nout,10) ufm,iopt
 10 FORMAT (a23,', ILLEGAL BCD VALUE (', 2A4,') FOR THE 4TH PARAMATER'  &
     ,      ' IN MODULE NORM')
 CALL mesage (-61,0,0)
 20 incru  = 1
 incrp  = 1
 icore  = korsz(iz)
 ibuf1  = icore - isysbf + 1
 ibuf2  = ibuf1 - isysbf
 icore  = ibuf2 - 1
 CALL gopen (matin ,iz(ibuf1),0)
 CALL gopen (matout,iz(ibuf2),1)
 mcb(1) = matin
 CALL rdtrl (mcb)
 ncol   = mcb(2)
 nrow   = mcb(3)
 nrow2  = 2*nrow
 itype  = mcb(5)
 iprec  = itype
 IF (iprec > 2) iprec = iprec - 2
 iunout = itype
 ipkot1 = itype
 ipkot2 = itype
 nrowp  = iprec*nrow
 nwords = nwds(itype)
 mwords = nrow*nwords
 kwords = mwords
 IF (iopt1 /= MAX) kwords = kwords + nrowp
 icrreq = kwords - icore
 IF (icrreq > 0) CALL mesage (-8,icrreq,isubnm)
 ivec   = mwords
 ivec1  = ivec + 1
 ivec2  = ivec + nrowp
 IF (iopt1 == MAX) GO TO 40
 mcb(5) = iprec
 ipkot1 = iprec
 ipkot2 = iprec
 DO  i= ivec1,ivec2
   z(i)   = 0.0
 END DO
 40 mcb(1) = matout
 mcb(2) = 0
 mcb(6) = 0
 mcb(7) = 0
 iu1    = 1
 iu2    = nrow
 
 xxmax  = 0.0
 DO  i = 1,ncol
   xx   = 0.0
   CALL unpack (*50,matin,z)
   ip1  = iu1
   ip2  = iu2
   xmax =-1.0
   GO TO 70
   50 ip1  = 1
   ip2  = 1
   xmax = 0.0
   DO  j = 1,nwords
     z(j) = 0.0
   END DO
   
   70 IF (iopt1 == isrss) GO TO 600
   IF (xmax  ==   0.0) GO TO 510
   
!     OPTION IS MAX
   
   SELECT CASE ( itype )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO 200
     CASE (    3)
       GO TO 300
     CASE (    4)
       GO TO 400
   END SELECT
   
   100 xmax = 0.0
   DO  j = 1,nrow
     x = ABS(z(j))
     IF (x > xmax) xmax = x
   END DO
   IF (xmax == 0.0) GO TO 510
   xx = xmax
   DO  j = 1,nrow
     z(j) = z(j)/xmax
   END DO
   GO TO 500
   
   200 dxmax = dzero
   DO  j = 1,nrow
     dx = DABS(zd(j))
     IF (dx > dxmax) dxmax = dx
   END DO
   IF (dxmax == dzero) GO TO 510
   xx = dxmax
   DO  j = 1,nrow
     zd(j) = zd(j)/dxmax
   END DO
   GO TO 500
   
   300 xmax = 0.0
   DO  j = 1,nrow2,2
     x = SQRT(z(j)*z(j) + z(j+1)**2)
     IF (x > xmax) xmax = x
   END DO
   IF (xmax == 0.0) GO TO 510
   xx = xmax
   DO  j = 1,nrow2,2
     z(j  ) = z(j  )/xmax
     z(j+1) = z(j+1)/xmax
   END DO
   GO TO 500
   
   400 dxmax = dzero
   DO  j = 1,nrow2,2
     dx = DSQRT(zd(j)*zd(j) + zd(j+1)**2)
     IF (dx > dxmax) dxmax = dx
   END DO
   IF (dxmax == dzero) GO TO 510
   xx = dxmax
   DO  j = 1,nrow2,2
     zd(j  ) = zd(j  )/dxmax
     zd(j+1) = zd(j+1)/dxmax
   END DO
   
   500 IF (xx > xxmax) xxmax = xx
   510 CALL pack (z,matout,mcb)
   CYCLE
   
!     OPTION IS SRSS
   
   600 IF (xmax == 0.0) CYCLE
   SELECT CASE ( itype )
     CASE (    1)
       GO TO 610
     CASE (    2)
       GO TO 630
     CASE (    3)
       GO TO 650
     CASE (    4)
       GO TO 670
   END SELECT
   
   610 DO  j = 1,nrow
     k = ivec + j
     z(k) = z(k) + z(j)*z(j)
   END DO
   CYCLE
   
   630 DO  j = 1,nrow
     k = ivec + j
     zd(k) = zd(k) + zd(j)*zd(j)
   END DO
   CYCLE
   
   650 k = ivec
   DO  j = 1,nrow2,2
     k = k + 1
     z(k) = z(k) + z(j)*z(j) + z(j+1)**2
   END DO
   CYCLE
   
   670 k = ivec
   DO  j = 1,nrow2,2
     k = k + 1
     zd(k) = zd(k) + zd(j)*zd(j) + zd(j+1)**2
   END DO
   
 END DO
 CALL CLOSE (matin, 1)
 IF (iopt1 == MAX) GO TO 760
 
 ip1 = iu1
 ip2 = iu2
 SELECT CASE ( iprec )
   CASE (    1)
     GO TO 710
   CASE (    2)
     GO TO 730
 END SELECT
 
 710 DO  i = ivec1,ivec2
   z(i) = SQRT(z(i))
   IF (z(i) > xxmax) xxmax = z(i)
 END DO
 GO TO 750
 
 730 DO  i = ivec1,ivec2
   zd(i) = DSQRT(zd(i))
   IF (zd(i) > xxmax) xxmax = zd(i)
 END DO
 
 750 CALL pack (z(ivec1),matout,mcb)
 
 760 CALL CLOSE (matout,1)
 CALL wrttrl (mcb)
 RETURN
END SUBROUTINE normal
