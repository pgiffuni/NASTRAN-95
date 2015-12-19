SUBROUTINE gfswch (file1,file2)
     
!     THE PURPOSE OF THIS SUBROUTINE IS TO INTERCHANGE THE NAMES OF
!     TWO FILES.  THIS IS ACCOMPLISHED BY THE DIRECT UPDATEING
!     OF THE FIAT AND THE FIST
 
 
 INTEGER, INTENT(IN)                      :: file1
 INTEGER, INTENT(IN OUT)                  :: file2
 EXTERNAL         lshift,rshift,andf,orf,complf
 INTEGER :: modnam(2),NAME(2),psave1,psave2,  &
     andf,orf,rshift,complf,UNIT,unit1,unit2,unt
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg / ufm,uwm,uim,sfm
 COMMON  /xfiat / ifiat(3)
 COMMON  /xfist / ifist(2)
 COMMON  /xpfist/ ipfist
 COMMON  /system/ sysbuf,nout,skip(21),icfiat
 DATA     modnam/ 4HGFSW,4HCH   /
 
 mask   = lshift(1,30) - 1
 mask   = lshift(rshift(mask,16),16)
 mask1  = complf(mask)
 mask2  = 32767
 mask3  = complf(mask2)
 nuniqe = ifiat(1)*icfiat + 3
 mxe    = ifiat(2)*icfiat + 3
 lastwd = ifiat(3)*icfiat + 3
 
!     LOCATE FILE POINTERS IN THE FIST
 
 nwd    = 2*ipfist   + 2
 nacent = 2*ifist(2) + 2
 nfiles = nacent - nwd
 psave1 = 0
 psave2 = 0
 DO  i = 1,nfiles,2
   IF (ifist(nwd+i) /= file1 .AND. ifist(nwd+i) /= file2) CYCLE
   IF (ifist(nwd+i) - file1 == 0) THEN
     GO TO    15
   END IF
   10 IF (ifist(nwd+i) - file2 == 0) THEN
     GO TO    20
   ELSE
     GO TO    25
   END IF
   15 psave1 = ifist(nwd+i+1) + 1
   iloc1  = i+nwd
   CYCLE
   20 psave2 = ifist(nwd+i+1) + 1
   iloc2  = i+nwd
   25 CONTINUE
 END DO
 
!     CHECK THAT FILES ARE IN FIST
 
 IF (psave1 == 0) CALL mesage (-1,file1,modnam)
 IF (psave2 == 0) CALL mesage (-1,file2,modnam)
 
!     SWITCH THE FIST POINTERS
 
 ifloc = ifist(iloc1+1)
 ifist(iloc1+1) = ifist(iloc2+1)
 ifist(iloc2+1) = ifloc
 
!     SWITCH FILE NAMES IN FIAT
 
 NAME(1)= ifiat(psave1+1)
 NAME(2)= ifiat(psave1+2)
 unit1  = andf(mask2,ifiat(psave1))
 unit2  = andf(mask2,ifiat(psave2))
 nwd    = icfiat*ifiat(3)-2
 ltu1   = andf(mask,ifiat(psave1))
 ltu2   = andf(mask,ifiat(psave2))
 ifiat(psave1  ) = orf(andf(ifiat(psave1),mask2 ),ltu2)
 ifiat(psave1+1) = ifiat(psave2+1)
 ifiat(psave1+2) = ifiat(psave2+2)
 ifiat(psave2  ) = orf(andf(ifiat(psave2),mask2),ltu1)
 ifiat(psave2+1) = NAME(1)
 ifiat(psave2+2) = NAME(2)
 
!     SWITCH STACKED DATA BLOCKS
 
 DO  i = 4,nwd,icfiat
   IF (psave1 == i .OR. psave2 == i) CYCLE
   IF (ifiat(i+1) == 0 .AND. ifiat(i+2) == 0) CYCLE
   UNIT = andf(mask2,ifiat(i))
   IF (UNIT /= unit1 .AND. UNIT /= unit2) CYCLE
   IF (UNIT == unit1) unt = unit2
   IF (UNIT == unit2) unt = unit1
   IF (i > nuniqe) GO TO 70
   
!     DATA BLOCK RESIDES IN UNIQUE PART OF FIAT
!     MOVE ENTRY TO BOTTOM
   
   IF (lastwd+icfiat <= mxe) GO TO 40
   WRITE  (nout,30) sfm
   30 FORMAT (a25,' 1021, FIAT OVERFLOW')
   CALL mesage (-37,0,modnam)
   40 ifiat(lastwd+1) = orf(andf(ifiat(i),mask3),unt)
   DO  k = 2,icfiat
     ifiat(lastwd+k) = ifiat(i+k-1)
   END DO
   lastwd   = lastwd   + icfiat
   ifiat(3) = ifiat(3) + 1
   
!     CLEAR OLD ENTRY IN UNIQUE PART
   
   ifiat(i) = andf(ifiat(i),mask2)
   j1 = i + 1
   j2 = i + icfiat - 1
   DO  k = j1,j2
     ifiat(k) = 0
   END DO
   CYCLE
   
!     DATA BLOCK RESIDES IN NON-UNIQUE PORTION OF FIAT
!     SWITCH UNIT NUMBERS
   
   70 ifiat(i) = orf(andf(ifiat(i),mask3),unt)
 END DO
 RETURN
END SUBROUTINE gfswch
