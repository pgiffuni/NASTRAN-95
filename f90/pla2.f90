SUBROUTINE pla2
!*****
! THIS ROUTINE IS THE SECOND FUNCTIONAL MODULE UNIQUE TO THE PIECE-WISE
! LINEAR ANALYSIS (PLA) RIGID FORMAT FOR THE DISPLACEMENT APPROACH.
 
! DMAP CALL...
 
! PLA2   DELTAUGV,DELTAPG,DELTAQG/UGV1,PGV1,QG1/V,N,PLACOUNT/ $
 
! CONCERNING DELTAUGV AND UGV1, THE ROUTINE WORKS AS FOLLOWS...
! DELTAUGV IS THE CURRENT INCREMENTAL DISPLACEMENT VECTOR IN THE PLA
! RIGID FORMAT LOOP AND UGV1 IS AN APPENDED FILE OF DISPLACEMENT VECTORS
! IF PLACOUNT .EQ. 1, THAT IS, THIS IS THE FIRST TIME PLA2 HAS BEEN
! CALLED IN THE PLA LOOP, THEN DELTAUGV IS COPIED ONTO UGV1.  IF
! PLACOUNT .GT. 1, THE PREVIOUS, OR LAST, DISPLACEMENT VECTOR IS READ
! INTO CORE FROM THE UGV1 DATA BLOCK, AND THE UGV1 FILE IS CLOSED WITH-
! OUT REWIND, THEN OPENED WITHOUT REWIND TO WRITE.  THE DELTAUGV VECTOR
! IS READ AN ELEMENT AT A TIME USING SUBROUTINE ZNTPKI AND ADDED TO
! THE VECTOR IN CORE.  THEN THE NEW DISPLACEMENT VECTOR IS WRITTEN ONTO
! THE UGV1 FILE.
 
! THEN THE PLA DMAP LOOP COUNTER PLACOUNT IS INCREMENTED.
 
! DELTAPG IS THE CURRENT INCREMENTAL LOAD VECTOR AND PGV1 IS THE
! CORRESPONDING MATRIX OF RUNNING SUM LOAD VECTORS.  PROCESSING IS
! SIMILAR TO THE ABOVE.  NOTE THAT NEITHER DATA BLOCK, LIKE THE TWO
! DISCUSSED ABOVE, CAN BE PURGED.
 
! DELTAQG IS THE CURRENT INCREMENTAL VECTOR OF SINGLE POINT CONSTRAINT
! FORCES AND QG1 IS THE APPENDED FILE OF VECTORS OF SPCF.  THESE TWO
! DATA BLOCKS ARE PROCESSED IDENTICALLY TO DELTAUGV AND UGV1 EXCECT
! THAT NO FATAL ERROR EXISTS IF ONE OR THE OTHER HAS BEEN PURGED.
!*****
 
 INTEGER :: bufsz  &
     ,                  buffr1             ,buffr2  &
     ,                  ofile              ,placnt  &
     ,                  eor                ,clsrw  &
     ,                  clsnrw             ,outrw  &
     ,                  eol                ,typea  &
     ,                  typeb              ,outnrw
 INTEGER :: inblk(15),oublk(15)
 
 DIMENSION NAME(2)            ,dummy(2)  &
     ,                  mcb(7)
 COMMON /BLANK/placnt
 COMMON   /system/  bufsz
 COMMON   /zzzzzz /  z(1)
 COMMON   /zntpkx/ a(4)               ,INDEX  &
     ,                  eol                ,idummy
 COMMON   /packx / typea              ,typeb  &
     ,                  ipack              ,jpack ,                  incpk
 COMMON   /unpakx/ iunpkb             ,iunpk  &
     ,                  junpk              ,incupk
 
 DATA               NAME /4HPLA2,4H    /
 DATA               inrw,outrw,outnrw,clsrw,clsnrw,eor/0,1,3,1,2,1/
 
! INITIALIZE
 
 izmax = korsz (z)
 buffr1 = izmax - bufsz
 buffr2 = buffr1 - bufsz
 left = buffr2 - 1
 iloop = 1
 ifile = 101
 ofile = 201
 
! OPEN INPUT FILE TO READ AND OUTPUT FILE TO WRITE (IF PLACNT .EQ. 1)
! OR TO READ (IF PLACNT .GT. 1)
 
 10 jfile = ifile
 inblk(1) = ifile
 oublk(1) = ofile
 DO  i = 2,7
   mcb(i) = 0
 END DO
 mcb(1) = ofile
 IF (placnt == 1) mcb(1) = ifile
 CALL rdtrl (mcb)
 CALL OPEN(*100,ifile,z(buffr1),inrw)
 CALL fwdrec(*9020,ifile)
 iopt = inrw
 IF (placnt == 1) iopt = outrw
 CALL OPEN(*110,ofile,z(buffr2),iopt)
 IF (placnt /= 1) GO TO 30
 
! THIS IS THE FIRST TIME THROUGH THE PLA LOOP.  COPY THE VECTOR ON THE
! INPUT FILE ONTO THE OUTPUT FILE.
 
 CALL fname (ofile,dummy)
 CALL WRITE (ofile,dummy,2,eor)
 CALL cpystr(inblk,oublk,0,0)
 CALL CLOSE (ofile,clsrw)
 CALL CLOSE(ifile,clsrw)
 GO TO 70
 
! THIS IS NOT THE FIRST PASS
 
 30 jfile = ofile
 CALL fwdrec(*9020,ofile)
 nrecs = placnt - 2
 IF (nrecs <= 0) GO TO 50
 DO  i = 1,nrecs
   CALL fwdrec(*9020,ofile)
 END DO
 50 mcb(1) = ofile
 CALL rdtrl (mcb)
 mcb(6) = 0
 mcb(7) = 0
 IF (left < mcb(3)) CALL mesage (-8,0,NAME)
 iunpkb = 1
 iunpk  = 1
 junpk  = mcb(3)
 incupk = 1
 CALL unpack(*9030,ofile,z)
 CALL CLOSE (ofile,clsnrw)
 CALL OPEN(*9010,ofile,z(buffr2),outnrw)
 
! READ THE INCREMENTAL VECTOR.  INTPK INITIALIZES AND ZNTPKI RETURNS
! ONE ELEMENT AT A TIME
 
 ktype = 1
 CALL intpk(*9040,ifile,0,ktype,0)
 
! READ AND ADD LOOP.
 
 60 CALL zntpki
 z(INDEX) = z(INDEX) + a(1)
 IF (eol == 0) GO TO 60
 
! ADDITION HAS BEEN COMPLETED
 
 CALL CLOSE (ifile,clsrw)
 
! WRITE VECTOR ON OUTPUT FILE IN PACKED FORMAT.
 
 typea = 1
 typeb = 1
 ipack = 1
 jpack = mcb(3)
 incpk = 1
 CALL pack(z,ofile,mcb)
 CALL CLOSE (ofile,clsrw)
 
! WRITE TRAILER
 
 70 mcb(1) = ofile
 CALL wrttrl (mcb)
 90 iloop = iloop + 1
 IF (iloop > 3) GO TO 200
 ifile = ifile + 1
 ofile = ofile + 1
 GO TO 10
 100 IF (iloop == 1  .OR.  iloop == 2) CALL mesage (-30,127,ifile)
 GO TO 90
 110 IF (iloop == 1  .OR.  iloop == 2) CALL mesage (-30,128,ofile)
 GO TO 90
 
! INCREMENT THE PLA DMAP LOOP COUNTER
 
 200 placnt = placnt + 1
 RETURN
 
! FATAL ERRORS
 
 9010 CALL mesage (-1,jfile,NAME)
 9020 CALL mesage (-2,jfile,NAME)
 9030 CALL mesage (-30,129,iloop+200)
 9040 CALL mesage (-30,130,iloop+100)
 RETURN
END SUBROUTINE pla2
