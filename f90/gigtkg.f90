SUBROUTINE gigtkg
     
 EXTERNAL        WRITE
 INTEGER :: gsize,scr2,scr3,trl(7),iz(1),sysbuf,out,buf1,  &
     buf2,nam(2),sdtab(6,5),ctype
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,out
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /gicom / spline,dum(8),ksize,gsize,scr1,scr2,scr3
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE    (z(1),iz(1))
 DATA    nam   / 4HGIGT,4HKG     /
 DATA    sdtab / 9, 9, 0, 9, 1, 9, 9, 0, 1, 9, 2, 3,  &
     9, 9, 0, 9, 9, 9, 9, 9, 0, 9, 1, 2,  &
     9, 9, 0, 9, 1, 2/
 
 ncore  = korsz(z) - 2*sysbuf
 buf1   = ncore
 buf2   = buf1 + sysbuf
 iti    = 1
 ito    = 1
 ii     = 1
 incr   = 1
 trl(1) = scr2
 trl(2) = 0
 trl(3) = gsize
 trl(4) = 2
 trl(5) = 1
 trl(6) = 0
 trl(7) = 0
 
!     BUILD A G BY K MATRIX PUT OUT SPLINE3 COLUMNS WHEN NECESSARY
 
 CALL gopen (scr2,z(buf1),1)
 CALL gopen (scr3,z(buf2),0)
 iss  = gsize + 1
 ncore= ncore - iss
 kcol = 0
 DO  i = 1,ksize
   IF (kcol < i) GO TO 20
   10 IF (kcol == i) GO TO 50
   nn   = 1
   z(1) = 0.0
   CALL pack (z,scr2,trl)
   CYCLE
   20 CALL READ (*30,*40,scr3,z(iss),ncore,0,nwr)
   GO TO 90
   30 kcol = ksize +1
   GO TO 10
   40 kst  = iz(iss+2)
   ctype= iz(iss+nwr-9)
   icm  = iz(iss+3)
   k    = sdtab(icm,ctype)
   IF(k == 9) GO TO 100
   kcol = kst + k
   GO TO 10
   
!     BUILD COLUMN FOR SPLINE CARD
   
   50 DO  j = 1,gsize
     z(j) = 0.0
   END DO
   nn   = gsize
   jj   = iss+4
   jjj  = iss+nwr-19
   DO  j = jj,jjj,3
     k    = iz(j) + iz(j+1) -1
     z(k) = z(j+2)
   END DO
   CALL pack (z,scr2,trl)
 END DO
 CALL CLOSE (scr2,1)
 CALL CLOSE (scr3,1)
 CALL wrttrl (trl)
 GO TO 120
 
!     ERROR MESSAGES
 
 90 CALL mesage (-8,ncore,nam)
 100 WRITE  (out,110) ufm,iz(iss),ctype,icm
 110 FORMAT (a23,' 2263, SPLINE3',i9,' FOR CAERO',i1,  &
     ' HAS ILLEGAL COMPONENT',i6)
 CALL mesage (-37,0,nam)
 120 RETURN
END SUBROUTINE gigtkg
