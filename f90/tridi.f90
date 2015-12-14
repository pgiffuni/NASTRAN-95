SUBROUTINE tridi (d,o,c,a,b,aa)
     
!     MODIFIED GIVENS REAL SYMMETRIC TRIDIAGONALIZATION
!     THIS ROUTINE IS CALLED ONLY BY VALVEC
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: d(1)
 DOUBLE PRECISION, INTENT(OUT)            :: o(1)
 DOUBLE PRECISION, INTENT(OUT)            :: c(1)
 REAL, INTENT(OUT)                        :: a(2)
 DOUBLE PRECISION, INTENT(OUT)            :: b(1)
 DOUBLE PRECISION, INTENT(IN)             :: aa(1)
 INTEGER :: savemr,ENTRY,rstrt,row,xentry,filcor,rot,row1,  &
     row2,rowp1,rowp2,sysbuf,mcb(7),count
 
 DIMENSION        vvcom(150)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /givn  /  title(1),mo,md,mr1,m1,m2,m3,m4,savemr,t10,ENTRY,  &
     t12(5),rstrt,row,t19,xentry
 COMMON /system/  sysbuf,nout,idummy(52),iprec
 COMMON /packx /  it1,it2,ii,jj,incr
 COMMON /unpakx/  it3,iii,jjj,incr1
 EQUIVALENCE      (vvcom(1),title(1)), (n,vvcom(101))
 DATA    count ,  MAX,mcb / 0, 10, 7*0/
 
 
!     DEFINITION OF VARIABLES
 
!     D       = LOCATION OF DIAGIONAL
!     O       = LOCATION OF OFF DIAGONAL
!     C       = LOCATION OF COSINES
!     A       = REST OF OPEN CORE
!     B       = O**2
!     SAVEMR
!     RSTRT
!     ROW
!     XENTRY
!     FILCOR
!     ROT
!     ROW1
!     ROW2
!     MO      = RESTART DATA - SINES AND COSINES
!     MD      = INPUT  MATRIX
!     MR1     = RESTART TAPE
!     M1      = SCRATCH TAPE
!     M2
!     M3
!     M4
!     MR2
!     MIDIN
!     COUNT   = NUMBER OF ROWS ROTATED
!     MAX     = NUMBER OF ROWS TO ROTATE BEFORE CHECKPOINTING
 
 
!     INITIALIZATION
 
 nz   = korsz(a)
 ibuf1= nz - sysbuf + 1
 ibuf2= ibuf1 - sysbuf
 nz   = nz - 2*sysbuf
 nzz  = nz/iprec
 nzsq = SQRT(FLOAT((nzz-1)*2))
 im1  = 1
 nm1  = n - 1
 nm2  = n - 2
 m3   = 305
 mss  = mr1
 ms1  = m1
 ms2  = m2
 ms3  = m3
 ms4  = m4
 
!     INITIALIZE  TRANSFORMATION ROUTINES
 
!     SICOX AND ROTAX ARE NOT USED ANY MORE. SEE SINC0S AND ROTATE
 
!     CALL SICOX (D,O,C)
!     CALL ROTAX (O,D,C)
 
 midin= n
 mr   = mr1
 
!     START AT THE BEGINNING
 
 row = 0
 
!     OPEN MD
 
 CALL gopen (md,a(ibuf1),0)
 CALL gopen (mr,a(ibuf2),1)
 
!     SET UP FOR UNPACK
 
 it3  = 2
 iii  = 1
 jjj  = n
 incr1= 1
 CALL unpack (*102,md,d)
 
!     COPY REST OF MD ONTO MR
 
 103 CONTINUE
 it1 = 2
 it2 = 2
 incr= 1
 k   = n - 1
 DO  i = 1,k
   iii = 0
   CALL unpack (*107,md,a)
   ii = iii
   jj = jjj
   106 CALL pack (a,mr,mcb)
   CYCLE
   107 ii = 1
   jj = 1
   a(1) = 0.0
   a(2) = 0.0
   GO TO 106
 END DO
 iii = 1
 jjj = n
 ii  = 1
 jj  = n
 GO TO 104
 102 DO  i = 1,n
   d(i) = 0.0D0
 END DO
 GO TO 103
 
!     END OF MATRIX MD
 
 104 CALL WRITE (mr,row,1,1)
 
!     ATTACH DIAGONALS
 
 CALL pack  (d,mr,mcb)
 CALL CLOSE (md,1)
 CALL CLOSE (mr,1)
 ms = mr
 CALL gopen (ms,a(ibuf1),0)
 
!     TRIDIAGONALIZATION PROCEDURE UNTIL THE MATRIX FITS IN CORE
 
 200 row  = row + 1
 rowp1= row + 1
 rowp2= row + 2
 it3  = 2
 iii  = rowp1
 CALL unpack (*201,ms,o(rowp1))
 GO TO 203
 201 DO  i = rowp1,n
   o(i) = 0.0D0
 END DO
 
!     FIND SINES AND COSINES
 
 203 CALL sinc0s (row,rot, d,o,c)
 CALL gopen (mo,a(ibuf2),im1)
 im1 = 3
 ii  = rowp2
 it1 = 2
 it2 = 2
 CALL pack (d(rowp2),mo,mcb)
 CALL CLOSE (mo,2)
 
!     WILL THE REST OF MATRIX FIT IN CORE
 
 IF ((n-rowp1)*(n-rowp1+1)/2+1 <= nzz) GO TO 225
 
!         (N-ROWP1)*(N-ROWP1  )     < (NZZ-1)*2
!                   (N-ROWP1  )     < SQRT((NZZ-1)*2) (=NZSQ)
!                    N              < NZSQ + ROWP1
!                    N-NZSQ         < ROWP1
!                    N-NZSQ         = NUMBER OF ROTAIONS NEEDED
 
!     NO-- MUST REST OF MATRIX BE ROTATED
 
 IF (rot == 0) GO TO 215
 count = count + 1
 IF (count == MAX) count = 0
 
!     ROTATE THE REST OF THE MATRIX
 
 midout = rowp1 + (n-rowp1+3)/4
 row1   = rowp2
 CALL gopen (ms3,a(ibuf2),1)
 
!     HERE THRU 217 WILL BE VERY TIME COMSUMING. THE ROTATION IS ONE
!     ROW AT A TIME. COMPUTE HOW MANY ROTATIONS NEEDED. IF TOO MANY,
!     ISSUE A USER FATAL MESSAGE AND GET OUT
 
 i = n - nzsq
 IF (i <= 25) GO TO 205
 j = (n*n - nzsq*nzsq)*iprec
 WRITE  (nout,204) ufm,n,n,i,j
 204 FORMAT (a23,' FROM GIVENS EIGENSOLVER - EXCESSIVE CPU TIME IS ',  &
     'NEEDED FOR TRIDIAGONALIZE THE DYNAMIC', /5X,  &
     'MATRIX, WHICH IS',i6,' BY',i6, 15X,1H(,i6,' LOOPS)', /5X,  &
     'RERUN JOB WITH',i8,' ADDITIONAL CORE WORDS, OR USE FEER,',  &
     ' OR OTHER METHOD')
 CALL mesage (-61,0,0)
 
!     FILL CORE WITH AS MUCH OF MATRIX AS POSSIBLE--UP TO ROW -ROW2-
 
 205 row2 = filcor(mss,ms2,iprec,row1,midin,n,a,nz,a(ibuf1))
 
!     ROTATE ROWS ROW1 TO ROW2
 
 CALL rotate (aa,row,row1,row2,o,d,c)
 
!     EMPTY THE ROTATED ROWS ONTO MS3 AND MS4
 
 CALL empcor (ms3,ms4,iprec,iprec,row1,midout,row2,n,a,a(ibuf2))
 row1 = row2 + 1
 IF (row2 < n) GO TO 205
 
!     SWITCH TAPES
 
 ms  = ms1
 ms1 = ms3
 ms3 = ms
 ms  = ms2
 ms2 = ms4
 ms4 = ms
 mss = ms1
 midin = midout
 215 DO  i = rowp1,n
   d(i) =  o(i)
 END DO
 ms = mss
 IF (row > midin) GO TO 217
 IF (rot ==     0) GO TO 200
 218 CALL gopen (ms,a(ibuf1),0)
 GO TO 200
 217 ms = ms2
 GO TO 218
 
!     TRIDIAGONALIZATION PROCEDURE WHEN MATRIX FITS IN CORE
 
 
!     FILL CORE WITH THE REST OF THE MATRIX
 
 225 row2 = filcor(mss,ms2,iprec,rowp2,midin,n,a,nz,a(ibuf1))
 na = 1
 CALL gopen (mo,a(ibuf2),3)
 GO TO 235
 230 row   = row + 1
 rowp1 = row + 1
 rowp2 = row + 2
 232 DO  i = rowp1,n
   o(i) = aa(na)
   na   = na + 1
 END DO
 234 CALL sinc0s (row,rot, d,o,c)
 
!     WRITE SINES ON MO
 
 ii  = rowp2
 it1 = 2
 it2 = 2
 CALL pack (d(rowp2),mo,mcb)
 235 IF (rot == 0) GO TO 236
 row1  = rowp2
 CALL rotate (aa(na),row,row1,row2,o,d,c)
 236 DO  i = rowp1,n
   d(i) = o(i)
 END DO
 IF (row /= nm2) GO TO 230
 
!     ALL DONE.
 
 d(n) = aa(na)
 o(n-1) = o(n)
 o(n  ) = 0.0D0
 CALL CLOSE (mo,3)
 DO  i = 1,n
   c(i) = d(i)
   b(i) = o(i)**2
 END DO
 xentry = -ENTRY
 rstrt  = 0
 savemr = 0
 RETURN
END SUBROUTINE tridi
