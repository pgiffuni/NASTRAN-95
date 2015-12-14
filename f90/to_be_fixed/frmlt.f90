SUBROUTINE frmlt (ifile,z,y,zm)
     
!     FEER MATRIX TRANSPOSE MULTIPLY  (SINGLE PREC)
!              T
!     Y = IFILE * Z        WHERE Z IS A VECTOR ALREADY IN CORE
!                          IFILE IS A GINO MATIRX FILE
 
!     LAST REVISED  11/91, BY C.CHAN/UNISYS
!     ADDITION OF A NEW TRANSPOSE MULTIPLY METHOD WHICH IS MORE
!     EFFECIENT, AND IS ALREADY GOOD FOR VECTORIZATION
 
!DB   LOGICAL          DEBUG
 
 INTEGER, INTENT(IN)                      :: ifile(7)
 REAL, INTENT(IN)                         :: z(1)
 REAL, INTENT(OUT)                        :: y(1)
 REAL, INTENT(IN)                         :: zm(1)
 REAL :: dp     ,sum
 DIMENSION  nam(2)
 COMMON  /unpakx/ ityp     ,ip       ,np       ,incr
 COMMON  /system/ ibuf     ,nout
 COMMON  /feerxx/ dum18(18),nzm
 COMMON  /zzzzzz/ iz(1)
 EQUIVALENCE      (dp,idp)
 DATA     nam   / 4HFRML   ,4HT    /
!DB   DATA     DEBUG , ITER     ,MAX    / .FALSE.   ,0       ,3     /
 
!DB   IF (.NOT.DEBUG) GO TO 20
!     ITER = ITER + 1
!     IF (ITER .GT. MAX) DEBUG = .FALSE.
!     IF (DEBUG) WRITE (NOUT,10) NZM,IFILE(5)
!  10 FORMAT ('  .... IN FRMLT DEBUG.   NZM,IFILE(5) =',2I8)
!  20 CONTINUE
 n    = ifile(2)
 ifl  = ifile(1)
 IF (ifile(7) < 0) ifl = -ifile(7)
 CALL REWIND (ifl)
 CALL skprec (ifl,1)
 IF (ifile(7) < 0) GO TO 50
 ityp = ifile(5)
 
!     NASTRAN ORIGIANL METHOD
 
 incr = 1
 DO  i = 1,n
   y(i) = 0.0
   ip   = 0
   CALL unpack (*40,ifl,zm(1))
   sum  = 0.0
   ii   = 0
   DO  j = ip,np
     ii   = ii + 1
     sum  = sum + zm(ii)*z(j)
   END DO
   y(i) = sum
 END DO
 GO TO 200
 
!     NEW METHOD, READ ONLY AND NO UNPACK
 
!     UNLIKE FRMLTA, IFL WAS UNPACKED FORWARD BY UNPSCR
 
 50 nrec = 0
!     NWDS = IFILE(5)
!DB   N20  = N - 20
!     IF (DEBUG) WRITE (NOUT,60) IFILE(5),NZM
!  60 FORMAT ('   /@60   NWDS,NZM =',2I8)
 ll2  = 0
 next = 1
 DO  i = 1,n
   IF (next < ll2) GO TO 100
   nrec = nrec + 1
!DB   IF (DEBUG) WRITE (NOUT,70) NREC,I
!  70 FORMAT ('  ...READING RECORD',I5,'.   I =',I7)
   CALL READ (*150,*80,ifl,zm,nzm,1,ll)
   CALL mesage (-8,0,nam)
!  50 LL2  = LL/NWDS
   80 ll2  = ll
!DB   IF (DEBUG) WRITE (NOUT,90) LL,NREC,LL2
!  90 FORMAT (1X,I10,'WORDS READ FROM RECORD',I5,'.   LL2 =',I10)
   next = 1
   100 dp   = zm(next)
   ii   = idp
   dp   = zm(next+1)
   jj   = idp
!DB   IF (DEBUG .AND. (I.LT.20 .OR. I.GT.N20)) WRITE (NOUT,110) I,II,JJ,
!     1                                                         NEXT
! 110 FORMAT ('   @110  I,II,JJ,NEXT =',4I8)
   IF (jj == ii) GO TO 130
   sum  = 0.0
   ll   = next + 1
   DO  j = ii,jj
     ll   = ll + 1
     sum  = sum + zm(ll)*z(j)
   END DO
   y(i) = sum
   GO TO 140
   130 y(i) = zm(next+2)*z(ii)
   140 next = next + jj - ii + 3
 END DO
 GO TO 200
 
 150 j = ifile(4)/10
 WRITE  (nout,160) nrec,i,n,j
 160 FORMAT ('*** TRY TO READ RECORD',i5,'.  I,N,IFILE(4) =',2I7,i5)
 CALL mesage (-2,ifl,nam)
 
 200 RETURN
END SUBROUTINE frmlt
