SUBROUTINE frmltd (ifile,dz,dy,zm)
     
!     FEER MATRIX TRANSPOSE MULTIPLY  (DOUBLE PREC)
!               T
!     DY = IFILE * DZ        WHERE DZ IS A VECTOR ALREADY IN CORE
!                            IFILE IS A GINO MATIRX FILE
 
!     LAST REVISED  11/91, BY C.CHAN/UNISYS
!     ADDITION OF A NEW TRANSPOSE MULTIPLY METHOD WHICH IS MORE
!     EFFECIENT, AND IS ALREADY GOOD FOR VECTORIZATION
 
!DB   LOGICAL          DEBUG
 
 INTEGER, INTENT(IN)                      :: ifile(7)
 DOUBLE PRECISION, INTENT(IN)             :: dz(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dy(1)
 DOUBLE PRECISION, INTENT(IN)             :: zm(1)
 DOUBLE PRECISION :: dp     ,dsum
 DIMENSION  idp(2)   ,nam(2)
 COMMON  /unpakx/ ityp     ,ip       ,np       ,incr
 COMMON  /system/ ibuf     ,nout
 COMMON  /feerxx/ dum18(18),nzm
 COMMON  /zzzzzz/ iz(1)
 EQUIVALENCE      (dp,idp(1))
 DATA     nam   / 4HFRML   ,4HTD    /
!DB   DATA     DEBUG , ITER     ,MAX     / .FALSE.  ,0      ,4       /
 
!DB   IF (.NOT.DEBUG) GO TO 20
!     ITER = ITER + 1
!     IF (ITER .GT. MAX) DEBUG = .FALSE.
!     IF (DEBUG) WRITE (NOUT,10) NZM,IFILE(5)
!  10 FORMAT ('  .... IN FRMLTD DEBUG.   NZM,IFILE(5) =',2I8)
!  20 CONTINUE
 n    = ifile(2)
 ifl  = ifile(1)
 IF (ifile(7) < 0) ifl = -ifile(7)
 CALL REWIND (ifl)
 CALL skprec (ifl,1)
 IF (ifile(7) < 0) GO TO 50
 ityp = ifile(5)
 
!     NASTRAN ORIGIANL METHOD
 
 incr  = 1
 DO  i = 1,n
   dy(i) = 0.0D+0
   ip    = 0
   CALL unpack (*40,ifl,zm(1))
   dsum  = 0.0D+0
   ii    = 0
   DO  j = ip,np
     ii    = ii + 1
     dsum  = dsum + zm(ii)*dz(j)
   END DO
   dy(i) = dsum
   40 CONTINUE
 END DO
 GO TO 200
 
!     NEW METHOD, READ ONLY AND NO UNPACK
 
!     UNLIKE FRMLTX, IFL WAS UNPACKED FORWARD BY UNPSCR
 
 50 nrec = 0
 nwds = ifile(5)
!DB   N20  = N - 20
!     IF (DEBUG) WRITE (NOUT,60) NWDS,NZM
!  60 FORMAT ('  /@60   NWDS,NZM =',2I8)
 ll2  = 0
 next = 1
 DO  i = 1,n
   IF (next < ll2) GO TO 100
   nrec = nrec + 1
!DB   IF (DEBUG) WRITE (NOUT,70) NREC,I
!  70 FORMAT ('  ...READING RECORD',I5,'.  I =',I7)
   CALL READ (*150,*80,ifl,zm,nzm,1,ll)
   CALL mesage (-8,0,nam)
   80 ll2  = ll/nwds
!DB   IF (DEBUG) WRITE (NOUT,90) LL,NREC,LL2
!  90 FORMAT (1X,I10,' WORDS READ FROM RECORD NO.',I5,'   LL2 =',I10)
   next = 1
   100 dp   = zm(next)
   ii   = idp(1)
   jj   = idp(2)
!DB   IF (DEBUG .AND. (I.LT.20 .OR. I.GT.N20)) WRITE (NOUT,110) I,II,JJ,
!     1                                                         NEXT
! 110 FORMAT ('   @110  I,II,JJ,NEXT =',4I8)
   IF (ii == jj) GO TO 130
   dsum = 0.0D+0
   ll   = next
   DO  j = ii,jj
     ll   = ll + 1
     dsum = dsum + zm(ll)*dz(j)
   END DO
   dy(i)= dsum
   GO TO 140
   130 dy(i)= zm(next+1)*dz(ii)
   140 next = next + jj - ii + 2
 END DO
 GO TO 200
 
 150 j = ifile(4)/10
 WRITE  (nout,160) nrec,i,n,j
 160 FORMAT ('0*** TRY TO READ RECORD',i5,'.   I,N,IFILE(4) =',2I7,i5)
 CALL mesage (-2,ifl,nam)
 
 200 RETURN
END SUBROUTINE frmltd
