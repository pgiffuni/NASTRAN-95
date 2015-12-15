SUBROUTINE switch
     
!     THE PURPOSE OF THIS MODULE IS TO INTERCHANGE THE NAMES OF THE
!     TWO INPUT FILES.  THIS IS ACCOMPLISHED BY THE DIRECT UPDATING
!     OF THE FIAT
 
 EXTERNAL        lshift,rshift,andf,orf,complf
 INTEGER :: file1,file2,modnam(2),NAME(2),psave1,psave2,  &
     andf,orf,rshift,complf,UNIT,unit1,unit2,unt
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /xfiat / ifiat(3)
 COMMON /xfist / ifist(2)
 COMMON /xpfist/ ipfist
 COMMON /BLANK / iparam
 COMMON /system/ sysbuf,nout,skip(21),icfiat
 DATA    file1 / 101/, file2 / 102/, modnam/ 4HSWIT,4HCH  /
 
 IF (iparam >= 0) RETURN
 mask2 = 32767
 mask3 = complf(mask2)
 mask  = lshift(1,30) - 1
 mask  = lshift(rshift(mask,16),16)
 mask1 = complf(mask)
 nuniqe= ifiat(1)*icfiat + 3
 mxe   = ifiat(2)*icfiat + 3
 lastwd= ifiat(3)*icfiat + 3
 
!     LOCATE FILE POINTERS IN THE FIST
 
 nwd    = 2*ipfist   + 2
 nacent = 2*ifist(2) + 2
 nfiles = nacent - nwd
 psave1 = 0
 psave2 = 0
 DO  i = 1,nfiles,2
   IF (ifist(nwd+i) /= file1 .AND. ifist(nwd+i) /= file2) CYCLE
   IF (ifist(nwd+i)-file1 == 0) THEN
     GO TO     3
   END IF
   2 IF (ifist(nwd+i)-file2 == 0) THEN
     GO TO     4
   ELSE
     GO TO    10
   END IF
   3 psave1 = ifist(nwd+i+1) + 1
   CYCLE
   4 psave2 = ifist(nwd+i+1) + 1
   10 CONTINUE
 END DO
 
!     CHECK THAT FILES ARE IN FIST
 
 IF (psave1 == 0) CALL mesage (-1,file1,modnam)
 IF (psave2 == 0) CALL mesage (-1,file2,modnam)
 
!     SWITCH FILE NAMES IN FIAT
 
 NAME(1) = ifiat(psave1+1)
 NAME(2) = ifiat(psave1+2)
 unit1   = andf(mask2,ifiat(psave1))
 unit2   = andf(mask2,ifiat(psave2))
 nwd     = icfiat*ifiat(3) - 2
 ltu1    = andf(mask,ifiat(psave1))
 ltu2    = andf(mask,ifiat(psave2))
 ifiat(psave1  ) = orf(andf(ifiat(psave1),mask2),ltu2)
 ifiat(psave1+1) = ifiat(psave2+1)
 ifiat(psave1+2) = ifiat(psave2+2)
 ifiat(psave2  ) = orf(andf(ifiat(psave2),mask2),ltu1)
 ifiat(psave2+1) = NAME(1)
 ifiat(psave2+2) = NAME(2)
 
!     SWITCH STACKED DATA BLOCKS
 
 DO  i = 4,nwd,icfiat
   IF (psave1 == i .OR. psave2 == i) CYCLE
   UNIT = andf(mask2,ifiat(i))
   IF (UNIT /= unit1 .AND. UNIT /= unit2) CYCLE
   IF (UNIT == unit1) unt = unit2
   IF (UNIT == unit2) unt = unit1
   IF (i   > nuniqe) GO TO 50
   
!     DATA BLOCK RESIDES IN UNIQUE PART OF FIAT
!     MOVE ENTRY TO BOTTOM
   
   IF (lastwd+icfiat <= mxe) GO TO 30
   WRITE  (nout,20) sfm
   20 FORMAT (a25,' 1021, FIAT OVERFLOW')
   CALL mesage (-37,0,modnam)
   30 ifiat(lastwd+1) = orf(andf(ifiat(i),mask3),unt)
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
   
   50 ifiat(i) = orf(andf(ifiat(i),mask3),unt)
 END DO
 RETURN
END SUBROUTINE switch
