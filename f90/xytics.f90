SUBROUTINE xytics (iout,out,ndevis,r1,r2,iskip,LOG,iflag)
     
!     THIS SUBROUTINE PERFORMS ONLY TIC COMPUTATIONS FOR XYDUMP.
 
 
 INTEGER, INTENT(OUT)                     :: iout(8)
 REAL, INTENT(OUT)                        :: out(8)
 INTEGER, INTENT(IN)                      :: ndevis
 REAL, INTENT(IN)                         :: r1
 REAL, INTENT(OUT)                        :: r2
 INTEGER, INTENT(IN)                      :: iskip
 INTEGER, INTENT(IN)                      :: LOG
 INTEGER, INTENT(IN OUT)                  :: iflag
 
 REAL :: length
 
 IF (LOG /= 0) GO TO 70
 IF (r1 == r2) r2 = r1 + 1.0
 div = ndevis
 IF (div <= 0.0) div = 5.0
 length = r2 - r1
 IF (iflag  /= 0  ) div = length
 IF (length <= 0.0) GO TO 50
 finc = 1.0001*length/div
 
!     CONVERT FINC TO SCIENTIFIC AND ROUND OFF TO 1 DIGIT (1 TO 10)
 
 ipower = 0
 IF (finc <  1.0) GO TO 20
 10 IF (finc < 10.0) GO TO 30
 ipower = ipower + 1
 finc   = finc/10.0
 GO TO 10
 20 ipower = ipower - 1
 finc = finc*10.0
 IF (finc < 1.0) GO TO 20
 30 iinc = 10
 IF (finc < 7.5) iinc = 5
 IF (finc < 3.5) iinc = 2
 IF (finc < 1.5) iinc = 1
 
!     ACTUAL INCREMENT
 
 finc = FLOAT(iinc)*10.0**ipower
 
!     COMPUTE FIRST DIVISION POINT
 
 nfirst = r1*10.0**(-ipower) + SIGN(0.555,r1)
 
!     GUARANTEE THAT TICKS WILL STEP THROUGH ZERO
 
 ntemp  = nfirst/iinc
 nfirst = ntemp*iinc
 first  = FLOAT(nfirst)*10.0**(ipower)
 
!     GET LOWEST VALUE OF FRAME
 
 IF (first <= r1) GO TO 35
 
!     CHECK ABAINST EPSILON DIFFERENCE.  SENSITIVE TO TRUNCATION
 
 length = finc*1.0E-4
 IF (first-r1 < length) first = first - length
 IF (first-r1 >= length) first = first - finc
 nfirst = first*10.0**(-ipower) + SIGN(0.5,r1)
 35 itics  = (r2-first)/finc + 1.5
 temp   = FLOAT(itics-1)*finc + first
 endv   = temp
 IF (endv >= r2) GO TO 37
 length = finc*2.0E-4
 IF (endv+length >= r2)  endv = endv + length
 IF (endv+length < r2)  endv = endv + finc
 itics = (endv-first)/finc + 0.5
 temp  = FLOAT(itics-1)*finc + first
 37 CONTINUE
 IF (endv-temp < finc/4.0) itics = itics - 1
 
!     FIND MAXIMUM NUMBER OF DIGITS
 
 last   = nfirst + iinc*itics
 last   = MAX0(IABS(last),IABS(nfirst))
 maxdig = 1
 40 IF (last < 10) GO TO 60
 maxdig = maxdig + 1
 last   = last/10
 GO TO 40
 
!     LENGTH = 0
 
 50 itics   = 0
 60 out(1)  = first
 out(2)  = finc
 out(3)  = endv
 iout(4) = maxdig
 iout(5) = ipower
 iout(6) = itics
 iout(7) = iskip
 RETURN
 
!     LOG SCALE - INITIAL LABELING CALCULATED
 
 70 first  = r1
 itics  = LOG
 endv   = r2
 finc   = 10.0
 maxdig = 1
 ipower = 0
 GO TO 60
 
END SUBROUTINE xytics
