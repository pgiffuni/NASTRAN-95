SUBROUTINE sofopn (b1,b2,b3)
     
!     READS THE SOF AND SYS COMMON BLOCKS FROM THE DIRECT ACCESS STORAGE
!     DEVICE, AND INITIALIZES THE POINTERS TO THE THREE BUFFERS NEEDED
!     BY THE SOF UTILITY SUBROUTINES
 
 
 INTEGER, INTENT(IN OUT)                  :: b1(1)
 INTEGER, INTENT(IN OUT)                  :: b2(1)
 INTEGER, INTENT(IN OUT)                  :: b3(1)
 LOGICAL :: first,opnsof
 INTEGER :: buf,dit,a,b,corwds,ginobl
 DIMENSION       NAME(2),iptr(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /machin/ mach
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / a(37)
 COMMON /sys   / b(6)
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10),STATUS,psswrd(2), first,opnsof
 COMMON /itemdt/ nitem,item(7,1)
 COMMON /system/ nbuff,nout
 COMMON /ginox / c(161),ginobl
 DATA    NAME  / 4HSOFO,4HPN  /
 DATA    ird   / 1 /
 
 IF (opnsof) GO TO 1000
 
!     CHECK IF THE OPEN CORE BUFFERS ARE LARGE ENOUGH AND DO NOT OVERLAP
 
 iptr(1) = corwds(buf,b1) + 2
 iptr(2) = corwds(buf,b2) + 2
 iptr(3) = corwds(buf,b3) + 2
 isiz    = korsz(buf)
 DO  i = 1,3
   IF (isiz-iptr(i) < nbuff-3) CALL mesage (-8,0,NAME)
 END DO
 DO  i = 1,2
   k = i + 1
   DO  j = k,3
     isiz = iptr(i) - iptr(j)
     IF (isiz <     0) isiz = -isiz
     IF (isiz < nbuff) CALL mesage (-8,0,NAME)
   END DO
 END DO
 a( 1) = iptr(1)
 a( 7) = iptr(2)
 a(15) = iptr(3)
 a(19) = iptr(1)
 
!     SET SOF BUFFER SIZE FROM /GINOX/
!     ON IBM USE /SYSTEM/ BECAUSE /GINOX/ IS IN SUPER LINK
 
 b(1) = ginobl
 IF (mach == 2 .OR. mach >= 5) b(1) = nbuff - 4
!WKBD 3/94      IF (MACH .EQ. 12) B(1) =NBUFF -28
 IF (first) CALL sofint (iptr(1),iptr(2),numb,ibl1)
 
!     READ AND INITIALIZE THE COMMON BLOCKS SYS AND SOF
 
 dit = iptr(1)
 CALL sofio (ird,1,buf(dit-2))
 DO  i = 1,4
   b(i) = buf(dit+24+i)
 END DO
 b(5) = buf(dit+46)
 b(6) = buf(dit+47)
 a(1) = iptr(1)
 a(2) = 0
 a(3) = 0
 a(4) = buf(dit+29)
 a(5) = buf(dit+30)
 a(6) = buf(dit+31)
 a(7) = iptr(2)
 DO  i = 8,14
   a(i) = 0
 END DO
 a(15) = iptr(3)
 a(16) = 0
 a(17) = 0
 a(18) = buf(dit+32)
 a(19) = iptr(1)
 a(20) = 0
 a(21) = 0
 a(22) = buf(dit+33)
 DO  i = 1,nfiles
   a(22+i) = buf(dit+33+i)
 END DO
 a(33) = buf(dit+44)
 a(34) = 0
 a(35) = 0
 a(36) = 0
 a(37) = buf(dit+45)
 
!     INITILIZE COMMON BLOCK ITEMDT
 
 nitem = buf(dit+100)
 k = 100
 DO  i = 1,nitem
   DO  j = 1,7
     item(j,i) = buf(dit+k+j)
   END DO
   k = k + 7
 END DO
 opnsof = .true.
 IF (.NOT. first) RETURN
 first  = .false.
 IF (numb == 0) RETURN
 
!     ADD THE NUMBER NUMB OF BLOCKS TO THE SUPERBLOCK WHOSE SIZE
!     NEEDED TO BE INCREASED
 
 DO  i = 1,numb
   CALL retblk (ibl1+i-1)
 END DO
 b(4) = b(4) - numb
 RETURN
 
!     ERROR MESSAGE
 
 1000 WRITE  (nout,1001) ufm
 1001 FORMAT (a23,' 6222 - ATTEMPT TO CALL SOFOPN MORE THAN ONCE ',  &
     'WITHOUT CALLING SOFCLS.')
 CALL sofcls
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE sofopn
