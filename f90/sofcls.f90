SUBROUTINE sofcls
     
!     WRITES OUT AT THE TERMINATION OF A MODULE ALL THE IN CORE BUFFERS
!     AND COMMON BLOCKS.
 
 LOGICAL :: ditup,mdiup,nxtup,opnsof,nxtrst
 INTEGER :: buf,a,b,filnam,filsiz,psswrd,dit,ditpbn,ditlbn, mdi,mdipbn,mdilbn
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / a(37)
 COMMON /sys   / b(6)
 COMMON /itemdt/ nitem,item(7,1)
 COMMON /system/ nbuff
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10),STATUS,psswrd(2), first,opnsof
 EQUIVALENCE     (dit   ,a(1) ),(ditpbn,a(2) ),(ditlbn,a(3) ),  &
     (mdi   ,a(15)),(mdipbn,a(16)),(mdilbn,a(17)),  &
     (nxt   ,a(19)),(nxtpbn,a(20)),(nxtlbn,a(21)),  &
     (ditup ,a(34)),(mdiup ,a(35)),(nxtup ,a(36)), (nxtrst,a(37))
 DATA    iwrt  / 2 /
 
 IF (.NOT.opnsof) RETURN
 IF (ditpbn == 0) GO TO 20
 IF (.NOT.ditup) GO TO 20
 CALL sofio (iwrt,ditpbn,buf(dit-2))
 ditup = .false.
 GO TO 40
 20 IF (nxtpbn == 0) GO TO 40
 IF (.NOT.nxtup) GO TO 40
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup = .false.
 40 IF (mdipbn == 0) GO TO 60
 IF (.NOT.mdiup) GO TO 60
 CALL sofio (iwrt,mdipbn,buf(mdi-2))
 mdiup = .false.
 
!     WRITE OUT COMMON BLOCKS.
 
 60 last = nbuff - 4
 DO  i = 1,last
   buf(dit+i) = 0
 END DO
 buf(dit+1) = psswrd(1)
 buf(dit+2) = psswrd(2)
 buf(dit+4) = nfiles
 DO  i = 1,nfiles
   buf(dit+ 4+i) = filnam(i)
   buf(dit+14+i) = filsiz(i)
   buf(dit+33+i) = a(22+i)
 END DO
 DO  i = 1,4
   buf(dit+24+i) = b(i)
 END DO
 buf(dit+29) = a(4)
 buf(dit+30) = a(5)
 buf(dit+31) = a(6)
 buf(dit+32) = a(18)
 buf(dit+33) = a(22)
 buf(dit+44) = a(33)
 nxtrst      = .false.
 buf(dit+45) = a(37)
 buf(dit+46) = b(5)
 buf(dit+47) = b(6)
 
 buf(dit+100) = nitem
 k = 100
 DO  i = 1,nitem
   DO  j = 1,7
     buf(dit+k+j) = item(j,i)
   END DO
   k = k + 7
 END DO
 ibl = 1
 DO  i = 1,nfiles
   buf(dit+3) = i
   CALL sofio (iwrt,ibl,buf(dit-2))
   ibl = ibl + filsiz(i)
 END DO
 CALL sofio (7, 0, 0)
 opnsof = .false.
 RETURN
END SUBROUTINE sofcls
