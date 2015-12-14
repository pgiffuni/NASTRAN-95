SUBROUTINE frmlta (ifile,z,y,zm)
     
!     LOWER TRIANGULAR TRANSPOSE WITH OFF-DIAGONAL SWITCH
!     SINGLE PRECISION VERSION
 
!     LAST REVISED  11/91, BY G.CHAN/UNISYS
!     ADDITIONAL OF A NEW METHODS WHICH IS MORE EFFICIENT, AND IS
!     ALREADY GOOD FOR VECTORIZATION
 
 
 INTEGER, INTENT(IN)                      :: ifile(7)
 REAL, INTENT(IN)                         :: z(1)
 REAL, INTENT(OUT)                        :: y(1)
 REAL, INTENT(OUT)                        :: zm(1)
 
 DIMENSION  nam(2)
 COMMON  /unpakx/ ityp    ,ip      ,np      ,incr
 COMMON  /feerxx/ dm18(18),nzm
 COMMON  /zzzzzz/ iz(1)
 COMMON  /system/ ibuf    ,nout
 EQUIVALENCE      (dp,idp)
 DATA     nam   / 4HFRML  ,4HTA  /
 
 n    = ifile(2)
 ifl  = ifile(1)
 IF (ifile(7) < 0) ifl = -ifile(7)
 CALL REWIND (ifl)
 IF (ifile(7) < 0) GO TO 30
 CALL skprec (ifl,1)
 ityp = ifile(5)
 
!     NASTRAN ORIGINAL METHOD
 
 incr = 1
 DO  i = 1,n
   y(i) = 0.0
   ip   = 0
   CALL unpack (*30,ifl,zm(1))
   IF (ip == i) zm(1) = -zm(1)
   sum  = 0.0
   ii   = 0
   DO  j = ip,np
     ii   = ii + 1
     sum  = sum - zm(ii)*z(j)
   END DO
   y(i) = sum
 END DO
 GO TO 150
 
!     NEW METHOD
 
!     UNLIKE FRMLT, IFL WAS UNPACKED BACKWARD FIRST, THEN FORWARD BY
!     UNPSCR/FEER3. SO WE SKIP BACKWARD PASS BEFORE READING DATA
 
 30 nrec = ifile(4)/10
 CALL skprec (ifl,nrec+1)
 nrec = 0
 ll2  = 0
 ntms = 1
 DO  i = 1,n
   IF (ntms < ll2) GO TO 50
   nrec = nrec + 1
   CALL READ (*100,*40,ifl,zm,nzm,1,ll)
   CALL mesage (-8,0,nam)
   40 ll2  = ll
   ntms = 1
   50 dp   = zm(ntms)
   ii   = idp
   IF (ii /= i) GO TO 120
   dp   = zm(ntms+1)
   jj   = idp
   zm(ntms+2) = -zm(ntms+2)
   sum  = 0.0
   ll   = ntms + 1
   DO  j = ii,jj
     ll   = ll + 1
     sum  = sum - zm(ll)*z(j)
   END DO
   y(i) = sum
   ntms = ntms + jj - ii + 3
 END DO
 GO TO 150
 
 100 j = ifile(4)/10
 WRITE  (nout,110) nrec,i,n,j
 110 FORMAT ('0*** TRY TO READ RECORD',i5,'.  I,N,IFILE(4) =',2I7,i5)
 CALL mesage (-2,ilf,nam)
 120 WRITE  (nout,130) ii,i
 130 FORMAT ('0*** II AND I MISMATCH =',2I8)
 CALL mesage (-37,0,nam)
 
 150 RETURN
END SUBROUTINE frmlta
