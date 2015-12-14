SUBROUTINE diagon
     
!     DMAP FUNCTIONAL MODULE
 
!     DIAGONAL  A / B / V,Y,OPT=COLUMN / V,Y,POWER $
 
!     INPUT  - A IS ANY MATRIX, EXCEPT RECTANGULAR AND ROW VECTOR
!            - OPT IS OUTPUT MATRIX TYPE V,Y FLAG
!            - POWER IS A VALUE TO WHICH THE REAL PART OF EACH ELEMENT
!              ON THE DIAGONAL OF A IS RAISED. (DEFAULT OF POWER IS 1.0)
!     OUTPUT - B IS A REAL SYMMETRIC MATRIX (OPT='SQUARE'), OR A COLUMN
!              VECTOR CONTAINING THE DIAGONAL OF A (OPT='COLUMN'), OR
!              A DIAGONAL MATRIX (OPT='DIAGONAL'
 
!     WRITTEN BY R. MITCHELL, CODE 324, GSFC, DECEMBER 7,1972
 
!     LAST MODIFIED BY G.CHAN/UNISYS   11/1991
!     TO MAKE SUERE  0.0**0 = 1.0, NOT 0.0
 
 INTEGER :: sysbuf,col,sq,ia(7),ib(7),NAME(2),opt(2)
 DOUBLE PRECISION :: d(2),dval(2),dcore(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm,swm
 COMMON /unpakx/  itypeu,iu,ju,incru
 COMMON /zntpkx/  a(4),ii,last
 COMMON /zblpkx/  val(4),jrow
!ZZ   COMMON /ZZDIAG/  CORE(1)
 COMMON /zzzzzz/  core(20000)
 COMMON /system/  ksystm(60)
 COMMON /BLANK /  param(3)
 EQUIVALENCE      (ksystm(1),sysbuf), (ksystm(2),nout),  &
     (ia(2),incol), (ia(3),inrow), (ia(4),iform),  &
     (ia(5),itype), (a(1),d(1)), (val(1),dval(1)),  &
     (param(1),opt(1)), (param(3),power), (core(1), dcore(1))
 DATA    col,sq/  4HCOLU,4HSQUA  /,  in1,iout / 101,201 /
 DATA    NAME  /  4HDIAG,4HONAL  /
 
 
!     CHECK FOR VALID PARAMETER.
 
 IF (opt(1) == sq .OR. opt(1) == col .OR. opt(1) == NAME(1)) GO TO 10
 WRITE  (nout,5) swm,opt
 5 FORMAT (a27,' 3300, INVALID PARAMETER ',2A4,  &
     ' SUPPLIED TO MODULE DIAGONAL, COLUMN SUBSTITUTED')
 opt(1) = col
 
!     GET INFO ON INPUT MATRIX
 
 10 ia(1) = in1
 CALL rdtrl (ia)
 
!     CHECK FOR PURGED INPUT.
 
 IF (ia(1) < 0) GO TO 210
 
!     CHECK FOR PROPER FORM OF MATRIX
 
 SELECT CASE ( iform )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 220
   CASE (    3)
     GO TO 20
   CASE (    4)
     GO TO 20
   CASE (    5)
     GO TO 20
   CASE (    6)
     GO TO 20
   CASE (    7)
     GO TO 220
   CASE (    8)
     GO TO 200
 END SELECT
 
!     SET OUTPUT CONTROL BLOCK TO MATCH INPUT AND REQUESTS.
 
 20 ib4 = 6
 IF (opt(1) == col) ib4 = 2
 IF (opt(1) /= NAME(1)) GO TO 25
 ib4 = 3
 opt(1) = col
 25 ib5 = 1
 IF (itype == 2 .OR. itype == 4) ib5 = 2
 CALL makmcb (ib,iout,inrow,ib4,ib5)
 
!     CHECK FOR SPECIAL CASES OF POWER PARAMETER.
 
!     CHECK FOR 1.0 = NO ARITHMETIC REQUIRED.
 
 IF (ABS(power-1.0)-1.0E-6 > 0.0) THEN
   GO TO    40
 END IF
 30 ipow = 1
 GO TO 100
 
!     CHECK FOR 0.5 = SQUARE ROOT
 
 40 IF (ABS(power-0.5)-1.0E-6 > 0.0) THEN
   GO TO    50
 END IF
 45 ipow = 2
 GO TO 100
 
!     CHECK FOR 2.0 = SQUARE
 
 50 IF (ABS(power-2.0)-1.0E-6 > 0.0) THEN
   GO TO    60
 END IF
 55 ipow = 3
 GO TO 100
 
!     CHECK FOR 0.0 = IDENTITY MATRIX
 
 60 IF (power == 0.0) THEN
   GO TO    65
 ELSE
   GO TO    70
 END IF
 65 ipow = 4
 GO TO 100
 
!     GENERAL CASE
 
 70 ipow = 5
 
!     DO OPEN CORE BOOKKEEPING
 
!     OBTAIN LENGTH OF OPEN CORE
 
 100 lcore = korsz(core)
 
!     NEED ROOM FOR 2 GINO BUFFERS
 
 IF (lcore < 2*sysbuf) GO TO 230
 
!     IF INPUT MATRIX IS A DIAGONAL MATRIX, NEED ADDITIONAL
!     ROOM FOR A FULL COLUMN
 
 IF (iform == 3 .AND. lcore < (2*sysbuf + ib(5)*inrow + 1)) GO TO 230
 ibuf = lcore - sysbuf + 1
 
!     OPEN INPUT FILE AND SKIP HEADER
 
 CALL gopen (ia,core(ibuf),0)
 
!     OPEN OUTPUT FILE AND WRITE HEADER
 
 nprec = ib(5)
 ibuf  = ibuf - sysbuf
 CALL gopen (ib,core(ibuf),1)
 
!     PRIME PACK ROUTINE IF COLUMN OUTPUT
 
 IF (opt(1) == col) CALL bldpk (nprec,nprec,iout,0,0)
 
!     READ INPUT MATRIX AND SEARCH COLUMNS FOR DIAGONAL ELEMENTS.
 
 DO  nowcol = 1,incol
   
!     CHECK IF THE INPUT MATRIX IS A DIAGONAL MATRIX (IFORM = 3)
   
   IF (iform /= 3) GO TO 118
   
!     UNPACK THE FULL COLUMN OF THE INPUT DIAGONAL MATRIX
   
   itypeu = nprec
   iu = 1
   ju = inrow
   incru = 1
   CALL unpack (*105,ia,core)
   GO TO 110
   105 jju = nprec*ju
   DO  i = 1,jju
     core(i) = 0.0
   END DO
   IF (ipow /= 4) GO TO 110
   IF (nprec ==  1) core (nowcol) = 1.0
   IF (nprec ==  2) dcore(nowcol) = 1.0D0
   110 ii = 0
   115 ii = ii + 1
   a(1) = core(ii)
   IF (nprec == 2) d(1) = dcore(ii)
   IF (opt(1) == sq) CALL bldpk (nprec,nprec,iout,0,0)
   GO TO 140
   
!     START A NEW COLUMN IF SYMMETRIC OUTPUT MATRIX.
   
   118 IF (opt(1) == sq) CALL bldpk (nprec,nprec,iout,0,0)
!WKBI 9/93
   inull = 0
   
!     START READING A COLUMN
   
!     NOTE THAT NULL INPUT COLUMN RESULTS IN NULL OUTPUT ELEMENT ONLY
!     IF POWER IS NOT ZERO.
   
   CALL intpk (*120,ia,0,itype,0)
   GO TO 130
   120 IF (ipow /= 4) GO TO 175
!WKBI 9/93
   inull = 1
   val(2)  = 0.0
   dval(2) = 0.0D0
   IF (nprec == 1)  val(1) = 1.0
   IF (nprec == 2) dval(1) = 1.0D0
   ii = nowcol
   GO TO 170
   
!     GET AN ELEMENT
   
   130 CALL zntpki
   
!     CHECK FOR DESIRED ELEMENT (ROW = COLUMN)
   
   IF (ii-nowcol < 0) THEN
     GO TO   132
   ELSE IF (ii-nowcol == 0) THEN
     GO TO   140
   ELSE
     GO TO   135
   END IF
   
!     CHECK FOR LAST NON-ZERO ELEMENT IN COLUMN.
   
   132 IF (last > 0) THEN
     GO TO   135
   ELSE
     GO TO   130
   END IF
   
!     SET ELEMENT VALUE TO 0. IF NOT IN COLUMN
   
   135 val(1)  = 0.
   dval(1) = 0.0D0
   GO TO 170
   
!     PROCESS RETURNED VALUE.
   
!     CHECK FOR PRECISION REQUIRED
   
   140 SELECT CASE ( nprec )
     CASE (    1)
       GO TO 150
     CASE (    2)
       GO TO 160
   END SELECT
   
!     SINGLE PRECISION PROCESSING OF REAL PART OF DIAGONAL ELEMENT
   
!     PERFORM REQUESTED OPERATION
   
   150 SELECT CASE ( ipow )
     CASE (    1)
       GO TO 152
     CASE (    2)
       GO TO 154
     CASE (    3)
       GO TO 156
     CASE (    4)
       GO TO 158
     CASE (    5)
       GO TO 159
   END SELECT
   152 val(1) = a(1)
   GO TO 170
   154 val(1) = SQRT(a(1))
   GO TO 170
   156 val(1) = a(1)*a(1)
   GO TO 170
   158 val(1) = 1.0
   GO TO 170
   159 val(1) = a(1)**power
   GO TO 170
   
!     DOUBLE PRECISION PROCESSING OF REAL PART OF DIAGONAL ELEMENT
   
!     PERFORM REQUESTED OPERATION
   
   160 SELECT CASE ( ipow )
     CASE (    1)
       GO TO 162
     CASE (    2)
       GO TO 164
     CASE (    3)
       GO TO 166
     CASE (    4)
       GO TO 168
     CASE (    5)
       GO TO 169
   END SELECT
   162 dval(1) = d(1)
   GO TO 170
   164 dval(1) = DSQRT(d(1))
   GO TO 170
   166 dval(1) = d(1)*d(1)
   GO TO 170
   168 dval(1) = 1.0D0
   GO TO 170
   169 dval(1) = d(1)**power
   
!     PACK COMPUTED VALUE INTO OUTPUT MATRIX
   
   170 jrow = nowcol
   IF (iform == 3) jrow = ii
   CALL zblpki
   
!     TEST FOR SPECIAL CASE OF DIAGONAL INPUT MATRIX (1 COLUMN).
   
   IF (iform == 3) GO TO 175
   
!     SKIP REST OF INPUT COLUMN IF NOT ON LAST ELEMENT.
   
!WKBI 9/93
   IF ( inull == 1 ) GO TO 171
   IF (last == 0) CALL skprec (in1,1)
!WKBI 9/93
   171 CONTINUE
   
!     TEST FOR SQUARE MATRIX CASE
!     FINISHED WITH COLUMN IF SQUARE MATRIX
   
   175 IF (opt(1) == sq) CALL bldpkn (iout,0,ib)
   
!     FINISHED WITH ONE OUTPUT ELEMENT.
   
   IF (iform == 3 .AND. ii < inrow) GO TO 115
 END DO
 
!     FINISH PACKING VECTOR IF COLUMN OUTPUT OPTION.
 
 IF (opt(1) == col) CALL bldpkn (iout,0,ib)
 
!     WRITE TRAILER IN FIAT.
 
 CALL wrttrl (ib)
 
!     FINISHED WITH ALL OF MATRIX, CLOSE UNITS
 
 CALL CLOSE(in1,1)
 CALL CLOSE(ib ,1)
 200 RETURN
 
!     ERROR MESSAGES
 
 210 RETURN
 
!     WRONG TYPE OF INPUT MATRIX = MSG.3016
 
 220 numm = -16
 GO TO 300
 
!     NOT ENOUGH CORE (MESSAGE 3008)
 
 230 numm = -8
 GO TO 300
 
 300 CALL mesage (numm,num,NAME)
 RETURN
 
END SUBROUTINE diagon
