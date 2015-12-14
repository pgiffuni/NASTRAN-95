SUBROUTINE sofint (ib1,ib2,numb,ibl1)
     
!     CALLED ONCE BY EVERY RUN USING THE SOF UTILITY SUBROUTINES.
!     SHOULD BE CALLED BEFORE ANY OF THEM IS CALLED.  IF THE SOF IS
!     NOT EMPTY, SOME SECURITY CHECKS WILL BE TAKEN CARE OF, AND THE
!     SOF COMMON BLOCKS WILL BE UPDATED AND WRITTEN OUT ON THE FIRST
!     BLOCK OF EACH OF THE SOF FILES.  IF THE SOF IS EMPTY, THE DIT
!     MDI, AND ARRAY NXT WILL BE INITIALIZED AND WRITTEN OUT ON THE
!     THIRD, FOURTH, AND SECOND BLOCKS OF THE FIRST FILE OF THE SOF,
!     AND THE SOF COMMON BLOCKS WILL BE INITIALIZED AND WRITTEN OUT
!     ON THE FIRST BLOCK OF EACH OF THE SOF FILES.
 
!     THE FIRST BLOCK OF EACH OF THE SOF FILES CONTAINS THE FOLLOWING
!     INFORMATION
!       WORD                   WORD                   WORD
!      NUMBER  CONTENTS       NUMBER  CONTENTS       NUMBER  CONTENTS
!      ------  --------       ------  --------       ------  --------
!       1- 2   PASSWORD          26   DIRSIZ            32   MDIBL
!          3   FILE NUMBER       27   SUPSIZ            33   NXTTSZ
!          4   NFILES            28   AVBLKS         34-43   NXTFSZ
!       5-14   FILNAM            29   DITSIZ            44   NXTCUR
!      15-24   FILSIZ            30   DITNSB            45   NXTRST
!         25   BLKSIZ            31   DITBL             46   HIBLK
!                                                       47   IFRST
 
!     STARTING AT LOCATION 100 THE CONTENTS OF THE ITEMDT COMMON BLOCK
!     ARE STORED
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ib1
 INTEGER, INTENT(IN OUT)                  :: ib2
 INTEGER, INTENT(OUT)                     :: numb
 INTEGER, INTENT(OUT)                     :: ibl1
 EXTERNAL        lshift,rshift,orf
 LOGICAL :: first
 INTEGER :: filnam,filsiz,STATUS,FILE,psswrd,orf,hiblk, buf,rshift,NAME(2)
 CHARACTER (LEN=31) :: sim
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm,sim
 COMMON /machin/ mac,ihalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10),STATUS,psswrd(2), first
 COMMON /system/ nbuff,nout,x1(36),nbpc,nbpw,ncpw
 COMMON /sys   / nsbuff,x4(3),hiblk,ifrst
 COMMON /itemdt/ nitem,item(7,1)
 DATA    ird,iwrt    /1, 2  /
 DATA    iempty,NAME /4H    ,4HSOFI,4HNT  /
 
 IF (ncpw <= 4) GO TO 5
 n = nbpw - nbpc*4
 DO  i = 1,10
   filnam(i) = lshift(rshift(filnam(i),n),n)
 END DO
 5 IF (nfiles <= 0) GO TO 1000
 IF (STATUS == 0) GO TO 250
 
!     THE SOF IS NOT EMPTY.  READ THE FIRST BLOCK OF THE FIRST SOF FILE
!     AND VERIFY THE SECURITY VARIABLES.
 
 FILE = filnam(1)
 CALL sofio (ird,1,buf(ib1-2))
 IF ((buf(ib1+1) /= psswrd(1)) .OR. (buf(ib1+2) /= psswrd(2))) GO TO 1050
 IF (buf(ib1+3) /= 1) GO TO 1060
 IF (buf(ib1+25) /= nsbuff) GO TO 1040
 
!     CHECK IF THE SPECIFIED NUMBER OF FILES AND THEIR SIZES IS ADEQUATE
 
 IF (buf(ib1+4) >= nfiles) GO TO 10
 MAX = buf(ib1+4) - 1
 GO TO 20
 10 MAX = nfiles - 1
 20 IF (MAX < 1) GO TO 50
 DO  i = 1,MAX
   IF (buf(ib1+14+i) == filsiz(i)) CYCLE
   FILE = filnam(i)
   GO TO 1070
 END DO
 
!     CHECK IF ALL SOF FILES HAVE THE CORRECT PASSWORD AND SEQUENCE
!     NUMBER
 
 MAX = MAX + 1
 ibl = 1
 DO  i = 2,MAX
   FILE = filnam(i)
   ibl  = ibl + filsiz(i-1)
   CALL sofio (ird,ibl,buf(ib1-2))
   IF ((buf(ib1+1) /= psswrd(1)) .OR. (buf(ib1+2) /= psswrd(2))) GO TO 1050
   IF (buf(ib1+3) /= i) GO TO 1060
 END DO
 CALL sofio (ird,1,buf(ib1-2))
 MAX = MAX - 1
 50 IF (buf(ib1+14+MAX+1) == filsiz(MAX+1)) GO TO 130
 maxnxt = 0
 IF (MAX < 1) GO TO 70
 DO  i = 1,MAX
   maxnxt = maxnxt+buf(ib1+33+i)
 END DO
 70 lastsz = (filsiz(MAX+1)-1)/buf(ib1+27)
 IF (filsiz(MAX+1)-1 == lastsz*buf(ib1+27)) GO TO 80
 lastsz = lastsz + 1
 80 maxnxt = maxnxt + lastsz
 IF (buf(ib1+33) > maxnxt) GO TO 1080
 maxold = maxnxt - lastsz + buf(ib1+33+MAX+1)
 IF (buf(ib1+33) /= maxold) GO TO 130
 IF (buf(ib1+14+MAX+1) > filsiz(MAX+1)) GO TO 1080
 lstsiz = MOD(buf(ib1+14+MAX+1)-2,buf(ib1+27)) + 1
 IF (lstsiz == buf(ib1+27)) GO TO 130
 
!     THE SIZE OF THE LAST SUPERBLOCK THAT WAS USED ON FILE (MAX+1)
!     SHOULD BE INCREASED.
 
 IF (filsiz(MAX+1)-buf(ib1+14+MAX+1) >= buf(ib1+27)-lstsiz) GO TO 90
 numb = filsiz(MAX+1) - buf(ib1+14+MAX+1)
 GO TO 100
 90 numb = buf(ib1+27) - lstsiz
 100 ibl1 = 0
 IF (MAX < 1) GO TO 120
 DO  i = 1,MAX
   ibl1 = ibl1 + filsiz(i)
 END DO
 120 ibl1 = ibl1 + buf(ib1+14+MAX+1) + 1
 GO TO 135
 130 numb = 0
 
!     UPDATE THE VARIABLE WHICH INDICATES THE NUMBER OF FREE BLOCKS ON
!     THE SOF.
 
 135 IF (nfiles-buf(ib1+4) < 0) THEN
   GO TO   140
 ELSE IF (nfiles-buf(ib1+4) == 0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 140 idiff = buf(ib1+14+nfiles) - filsiz(nfiles)
 MIN   = nfiles + 1
 last  = buf(ib1+4)
 DO  i = MIN,last
   idiff = idiff + buf(ib1+14+i)
 END DO
 GO TO 190
 160 idiff = buf(ib1+14+nfiles) - filsiz(nfiles)
 GO TO 190
 170 ihere1 = buf(ib1+4)
 idiff  = buf(ib1+14+ihere1) - filsiz(ihere1)
 MIN    = buf(ib1+4) + 1
 DO  i = MIN,nfiles
   idiff = idiff - filsiz(i)
 END DO
 190 buf(ib1+28) = buf(ib1+28) - idiff
 
!     IF NO ITEM STRUCTURE IS ON THE SOF (THE SOF WAS CREATED BEFORE
!     LEVEL 17.0) THEN USE THE LEVEL 16.0 ITEM STRUCTURE.
 
 IF (buf(ib1+100) > 0 .AND. buf(ib1+100) <= 100) GO TO 198
 WRITE (nout,6235) uwm
 buf(ib1+ 47) = 3
 buf(ib1+100) = 18
 k = 100
 DO  i = 1,18
   DO  j = 1,7
     buf(ib1+k+j) = item(j,i)
   END DO
   k = k + 7
 END DO
 GO TO 200
 
!     CHECK IF THE DIRECTORY SIZE HAS BEEN CHANGED
 
 198 IF (nitem == buf(ib1+100)) GO TO 200
 WRITE (nout,6233) uwm
 
!     UPDATE THE COMMON BLOCKS USED BY THE SOF UTILITY SUBROUTINES.
 
 200 buf(ib1+4) = nfiles
 DO  i = 1,nfiles
   buf(ib1+4 +i) = filnam(i)
   buf(ib1+14+i) = filsiz(i)
   buf(ib1+33+i) = (filsiz(i)-1)/buf(ib1+27)
   IF (filsiz(i)-1 == buf(ib1+33+i)*buf(ib1+27)) CYCLE
   buf(ib1+33+i) = buf(ib1+33+i) + 1
 END DO
 
!     WRITE THE UPDATED ARRAY A ON THE FIRST BLOCK OF EACH OF THE SOF
!     FILES.
 
 ibl = 1
 DO  i = 1,nfiles
   buf(ib1+3) = i
   CALL sofio (iwrt,ibl,buf(ib1-2))
   ibl = ibl + filsiz(i)
 END DO
 GO TO 340
 
!     THE SOF IS EMPTY.  INITIALIZE THE SOF COMMON BLOCKS WHICH ARE
!     STORED IN THE ARRAY A.
!     CHECK IF THE NASTRAN BUFFER SIZE IS LARGE ENOUGH
 
 250 MIN = 100 + 7*nitem + (nbuff-nsbuff)
 IF (nbuff < MIN) GO TO 1090
 last  = nbuff - 4
 hiblk = 0
 ifrst = 3
 DO  i = 1,last
   buf(ib1+ i) = 0
 END DO
 buf(ib1+ 1) = psswrd(1)
 buf(ib1+ 2) = psswrd(2)
 buf(ib1+25) = nsbuff
 buf(ib1+26) = nitem + ifrst - 1
 buf(ib1+27) = 2*(buf(ib1+25)-1)
 buf(ib1+28) = -4
 DO  i = 1,nfiles
   buf(ib1+28) = buf(ib1+28) + filsiz(i)
 END DO
 buf(ib1+29) = 0
 buf(ib1+30) = 0
 buf(ib1+31) = 3
 buf(ib1+32) = 4
 buf(ib1+33) = 1
 buf(ib1+44) = 1
 buf(ib1+45) = 0
 buf(ib1+46) = 4
 buf(ib1+47) = ifrst
 
 buf(ib1+100) = nitem
 k = 100
 DO  i = 1,nitem
   DO  j = 1,7
     buf(ib1+k+j) = item(j,i)
   END DO
   k = k + 7
 END DO
 
!     INITIALIZE THE ARRAY NXT AND WRITE IT ON THE SECOND BLOCK OF THE
!     FIRST SOF FILE.
 
 DO  i = 1,last
   buf(ib2+i) = 0
 END DO
 IF (buf(ib1+27)+1 > filsiz(1)) GO TO 302
 MAX = buf(ib1+25) - 1
 buf(ib2+MAX+1) = lshift(buf(ib1+27)+1,ihalf)
 buf(ib2+1) = buf(ib1+27) + 1
 GO TO 308
 302 IF (MOD(filsiz(1),2) == 1) GO TO 304
 MAX = filsiz(1)/2
 GO TO 306
 304 MAX = (filsiz(1)-1)/2
 buf(ib2+MAX+1) = lshift(filsiz(1),ihalf)
 306 buf(ib2+1) = filsiz(1)
 308 buf(ib2+1) = orf(buf(ib2+1),lshift(5,ihalf))
 buf(ib2+2) = 0
 buf(ib2+3) = 6
 DO  i = 4,MAX
   buf(ib2+i) = 2*i
   buf(ib2+i) = orf(buf(ib2+i),lshift(2*i-1,ihalf))
 END DO
 CALL sofio (iwrt,1,buf(ib2-2))
 CALL sofio (iwrt,2,buf(ib2-2))
 
!     INITIALIZE THE DIT AND WRITE IT ON THE THIRD BLOCK OF THE FIRST
!     SOF FILE.
 
 DO  i = 1,last
   buf(ib2+i) = iempty
 END DO
 CALL sofio (iwrt,3,buf(ib2-2))
 
!     INITIALIZE THE MDI AND WRITE IT ON THE FOURTH BLOCK OF THE FIRST
!     SOF FILE.
 
 DO  i = 1,last
   buf(ib2+i) = 0
 END DO
 CALL sofio (iwrt,4,buf(ib2-2))
 numb = 0
 GO TO 200
 
!     PRINT MESSAGE INDICATING THE STATUS OF THE CURRENT SOF FILES.
 
 340 WRITE (nout,360) sim,nfiles
 DO  i = 1,nfiles
   WRITE (nout,370) i,filsiz(i)
 END DO
 WRITE  (nout,380) buf(ib1+25)
 360 FORMAT (a31,' 6201,',i3,' FILES HAVE BEEN ALLOCATED TO THE SOF ',  &
     'WHERE --')
 370 FORMAT (18H     size of FILE ,i2,3H = ,i10,7H blocks)
 380 FORMAT (32H     AND where a BLOCK contains ,i4,6H words)
 RETURN
 
!     ERROR MESSAGES.
 
 1000 WRITE  (nout,1001) sfm
 1001 FORMAT (a25,' 6202.  THE REQUESTED NO. OF FILES IS NON POSITIVE.')
 CALL mesage (-37,0,NAME(1))
 RETURN
 
 1040 i = (nbuff-nsbuff) + buf(ib1+25)
 WRITE  (nout,1041) ufm,i
 1041 FORMAT (a23,' 6205, SUBROUTINE SOFINT - THE BUFFER SIZE HAS BEEN',  &
     ' MODIFIED.', /30X, 'THE CORRECT NASTRAN PARAMETER IS BUFFSIZE = ',i6)
 GO TO 1082
 
 1050 WRITE  (nout,1051) ufm,FILE
 1051 FORMAT (a23,' 6206, SUBROUTINE SOFINT - WRONG PASSWORD ON SOF ',  &
     'FILE ',a4,1H.)
 GO TO 1082
 
 1060 WRITE  (nout,1061) ufm,FILE
 1061 FORMAT (a23,' 6207, SUBROUTINE SOFINT - THE SOF FILE ',a4,  &
     ' IS OUT OF SEQUENCE.')
 GO TO 1082
 
 1070 WRITE  (nout,1071) ufm,FILE
 1071 FORMAT (a23,' 6208, SUBROUTINE SOFINT - THE SIZE OF THE SOF FILE '  &
     ,       a4,' HAS BEEN MODIFIED.')
 GO TO 1082
 
 1080 WRITE  (nout,1081) ufm,FILE
 1081 FORMAT (a23,' 6209, SUBROUTINE SOFINT - THE NEW SIZE OF FILE ',a4,  &
     ' IS TOO SMALL.')
 1082 CALL mesage (-61,0,0)
 
 1090 WRITE  (nout,1091) ufm,MIN
 1091 FORMAT (a23,' 6234, THE NASTRAN BUFFER SIZE IS TO SMALL FOR THE',  &
     ' SOF FILE.', /30X,'MINIMUM BUFFER SIZE IS ',i10)
 GO TO 1082
 
 6233 FORMAT (a25,' 6233, THE ITEM STRUCTURE HAS BEEN CHANGED FOR THE ',  &
     'SOF.', /32X,'NEW CAPABILITIES USING THESE ITEMS MAY NOT ',  &
     'BE USED WITH THIS SOF.')
 
 6235 FORMAT (a25,'6235, THE OLD SOF CONTAINS NO ITEM STRUCTURE ',  &
     'INFORMATION.', /27X,'THE LEVEL 16.0 ITEM STRUCTURE WILL ', 'BE USED.')
 
END SUBROUTINE sofint
