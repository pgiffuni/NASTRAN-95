SUBROUTINE table v (*,in,ll,trl,NAME,p4,ibuf,z5)
     
!     TABLE-V IS CALLED ONLY BY INPUT5 TO GENERATE A GINO TABLE
!     DATA BLOCK IN 'OUT' FROM AN INPUT FILE 'IN' - A REVERSE PROCESS
!     OF TABLE-5.
!     THE INPUT FILE WAS FORTRAN WRITTEN, FORMATTED OR UNFORMATTED
 
!     IN     = INPUT FILE, INTEGERS
!     LL     = (200+LL) IS THE OUTPUT FILE, INTEGER
!     TRL    = AN ARRAY OF 7 WORDS FOR TRAILER
!     NAME   = ORIGINAL FILE NAME FROM INPUT FILE, 2 BCD WORDS, PLUS 1
!     P4     = 0, INPUT FILE WAS WRITTEN UNFORMATTED, BINARY, INTEGER
!            = 1, INPUT FILE WAS WRITTEN FORMATTED, ASCII, INTEGER
!     IBUF   = OPEN CORE AND GINO BUFFER POINTER, INTEGER
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: in
 INTEGER, INTENT(IN OUT)                  :: ll
 INTEGER, INTENT(OUT)                     :: trl(7)
 INTEGER, INTENT(OUT)                     :: NAME(3)
 INTEGER, INTENT(IN OUT)                  :: p4
 INTEGER, INTENT(IN)                      :: ibuf
 CHARACTER (LEN=5), INTENT(IN)            :: z5(1)
 LOGICAL :: debug
 INTEGER :: sysbuf, z, out, namex(2),sub(2), END,tble,fuf,fu(2)
 REAL :: rz(1),z4(2)
 DOUBLE PRECISION :: dz
 CHARACTER (LEN=1) :: z1,i1,r1,b1,d1,f1
 CHARACTER (LEN=5) :: z5l,end5
 CHARACTER (LEN=10) :: z10
 CHARACTER (LEN=15) :: z15
 COMMON /system/  sysbuf,nout
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (z1,z5l), (z(1),rz(1)), (dz,z4(1))
 DATA    i1,r1,   b1,d1,f1  / 'I', 'R', '/', 'D', 'X'    /
 DATA    fu,      END,end5  / 2H  ,2HUN, 4H*END, ' *END' /
 DATA    sub,     tble      / 4HTABL,4HEV  ,     4HTBLE  /
 DATA    debug              / .false.                    /
 
 IF (debug) WRITE (nout,10)
 10   FORMAT (///,' *** IN TABLE-V, DEBUG ***')
 kore  = ibuf-1
 kore9 = (kore/9)*9
 out   = 200+ll
 ll    = ll+1
 kount = 0
 
!     OPEN GINO OUTPUT FILE AND WRITE A FILE HEADER
 
 CALL OPEN (*180,out,z(ibuf),1)
 CALL fname (out,namex)
 CALL WRITE (out,namex,2,1)
 IF (debug) WRITE (nout,20) namex
 20   FORMAT (/5X,'GENERATING...',2A4,/)
 NAME(3) = tble
 IF (p4 == 1) GO TO 40
 
!     UNFORMATED READ
 
 30   READ (in,ERR=150,END=130) ln,(z(j),j=1,ln)
 IF (ln > kore) GO TO 170
 IF (ln == 1 .AND. z(1) == END) GO TO 130
 CALL WRITE  (out,z(1),ln,1)
 kount = kount+1
 GO TO 30
 
!     FORMATTED READ
 
 40   READ  (in,50,ERR=150,END=130) ln,(z5(j),j=1,ln)
 50   FORMAT (i10,24A5,/,(26A5))
 IF (ln > kore) GO TO 170
 IF (ln == 1 .AND. z5(1) == end5) GO TO 130
 IF (ln <= -1) GO TO 130
 lb = (ln*5)/4+1
 k  = 0
 l  = 1
 60   IF (l > ln) GO TO 120
 k  = k+1
 z5l= z5(l)
 IF (z1 == i1) GO TO 90
 IF (z1 == r1) GO TO 100
 IF (z1 == b1) GO TO 70
 IF (z1 == f1) GO TO 80
 IF (z1 == d1) GO TO 110
 WRITE  (nout,65) z5l
 65   FORMAT (/,' SYSTEM ERROR/TABLEV @65  Z5L=',a5)
 GO TO 150
 
!     BCD
 
 70   READ (z5l,75) z(lb+k)
 75   FORMAT (1X,a4)
 
!     FILLER
 
 80   l = l+1
 GO TO 60
 
!     INTEGER
 
 85   FORMAT (3A5)
 90   WRITE  (z10,85) z5(l),z5(l+1)
 READ   (z10,95) z(lb+k)
 95   FORMAT (1X,i9)
 l = l+2
 GO TO 60
 
!     REAL, SINGLE PRECISION
 
 100  WRITE  (z15, 85) z5(l),z5(l+1),z5(l+2)
 READ   (z15,105) rz(lb+k)
 105  FORMAT (1X,e14.7)
 l = l+3
 GO TO 60
 
!     REAL, DOUBLE PRECISION
 
 110  WRITE (z15, 85) z5(l),z5(l+1),z5(l+2)
 READ  (z15,115) dz
 115  FORMAT (1X,d14.7)
 rz(lb+k  ) = z4(1)
 rz(lb+k+1) = z4(2)
 k = k+1
 l = l+3
 GO TO 60
 
 120  IF (k <= 0) GO TO 40
 CALL WRITE (out,z(lb+1),k,1)
 kount = kount+1
 GO TO 40
 
!     ALL DONE.
!     CLOSE OUTPUT GINO FILE AND WRITE TRAILER
 
 130  CALL CLOSE (out,1)
 IF (debug) WRITE (nout,135) trl(2),kount
 135  FORMAT (/,' DEBUG ECHO - OLD AND NEW COLUMN COUNTS =',2I5)
 trl(1) = out
 trl(2) = kount
 CALL wrttrl (trl)
 fuf = fu(1)
 IF (p4 == 0) fuf = fu(2)
 WRITE  (nout,140) fuf,namex
 140  FORMAT (/5X,'DATA TRANSFERED SUCCESSFULLY FROM ',a2,'FORMATTED ',  &
     'TAPE TO GINO OUTPUT FILE ',2A4)
 GO TO 200
 
!     ERROR
 
 150  CALL CLOSE (out,1)
 WRITE  (nout,160) namex
 160  FORMAT (//5X,'ERROR IN READING INPUT TAPE IN TABLEV. NO ',2A4,  &
     /5X,'FILE GENERATED')
 GO TO 200
 170  CALL mesage (8,0,sub)
 GO TO 200
 180  CALL mesage (1,out,sub)
 
 200  RETURN 1
END SUBROUTINE table v
