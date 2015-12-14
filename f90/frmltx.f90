SUBROUTINE frmltx (ifile,dz,dy,zm)
     
!     LOWER TRIANGULAR TRANSPOSE WITH OFF-DIAGONAL SWITCH
!     DOUBLE PRECISION VERSION
 
!     LAST REVISED  11/91, BY G.CHAN/UNISYS
!     ADDITIONAL OF A NEW METHOD WHICH IS MORE EFFICIENT, AND IS
!     ALREADY GOOD FOR VECTORIZATION
 
 
 INTEGER, INTENT(IN)                      :: ifile(7)
 DOUBLE PRECISION, INTENT(IN)             :: dz(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dy(1)
 DOUBLE PRECISION, INTENT(OUT)            :: zm(1)
 DOUBLE PRECISION :: dp     ,dsum
 DIMENSION  idp(2)  ,nam(2)
 COMMON  /unpakx/ ityp    ,ip      ,np      ,incr
 COMMON  /feerxx/ dm18(18),nzm
 COMMON  /zzzzzz/ iz(1)
 COMMON  /system/ ibuf    ,nout
 EQUIVALENCE      (dp,idp(1))
 DATA     nam   / 4HFRML  ,4HTX    /
 
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
   dy(i)= 0.0D+0
   ip   = 0
   CALL unpack (*30,ifl,zm(1))
   IF (ip == i) zm(1) = -zm(1)
   dsum = 0.d0
   ii   = 0
   DO  j = ip,np
     ii   = ii + 1
     dsum = dsum - zm(ii)*dz(j)
   END DO
   dy(i)= dsum
 END DO
 GO TO 150
 
!     NEW METHOD
 
!     UNLIKE FRMLTD, IFL WAS UNPACKED BACKWARD FIRST, THEN FORWARD BY
!     UNPSCR/FEER3. SO WE SKIP BACKWARD PASS BEFORE READING DATA
 
 30 nrec = ifile(4)/10
 CALL skprec (ifl,nrec+1)
 nwds = ifile(5)
 nrec = 0
 ll2  = 0
 ntms = 1
 DO  i = 1,n
   IF (ntms < ll2) GO TO 50
   nrec = nrec + 1
   CALL READ (*100,*40,ifl,zm,nzm,1,ll)
   CALL mesage (-8,0,nam)
   40 ll2  = ll/nwds
   ntms = 1
   50 dp   = zm(ntms)
   ii   = idp(1)
   jj   = idp(2)
   IF (ii /= i) GO TO 120
   zm(ntms+1) = -zm(ntms+1)
   dsum = 0.0D+0
   ll   = ntms
   DO  j = ii,jj
     ll   = ll + 1
     dsum = dsum - zm(ll)*dz(j)
   END DO
   dy(i)= dsum
   ntms = ntms + jj - ii + 2
 END DO
 GO TO 150
 
 100 j = ifile(4)/10
 WRITE  (nout,110) nrec,i,n,j
 110 FORMAT ('0*** TRY TO READ RECORD',i5,'.  I,N,IFILE(4) =',2I7,i5)
 CALL mesage (-2,ifl,nam)
 120 WRITE  (nout,130) ii,i
 130 FORMAT ('0*** II AND I MISMATCH =',2I8)
 CALL mesage (-37,0,nam)
 
 150 RETURN
END SUBROUTINE frmltx
