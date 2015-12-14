INTEGER FUNCTION forfil (NAME)
     
!     FORFIL RETURNS THE LOGICAL UNIT TO WHICH NAME IS ASSIGNED.
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 EXTERNAL        andf
 INTEGER :: andf    ,exfiat,fiat  ,fist  ,sysout
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm     ,uwm   ,uim   ,sfm
 COMMON /system/ sysbuf  ,sysout
 COMMON /xxfiat/ exfiat(1)
 COMMON /xfiat / fiat(1)
 COMMON /xfist / mfist   ,nfist ,fist(1)
 
!     SEARCH FIST FOR NAME. ERROR IF NOT FOUND.
 
 nn = 2*nfist - 1
 DO  i = 1,nn,2
   IF (fist(i) == NAME) GO TO 2010
 END DO
 WRITE  (sysout,2002) sfm,NAME,NAME
 2002 FORMAT (a25,' 2179, ERROR DETECTED IN FUNCTION FORFIL',a4,i4,  &
     ' NOT IN FIST.')
 CALL mesage (-61,0,0)
 forfil = 0
 RETURN
 
!     PICK UP UNIT FROM /XXFIAT/ OR /XFIAT/ AND RETURN.
 
 2010 j = fist(i+1)
 IF (j < 0) THEN
   GO TO  2013
 ELSE IF (j == 0) THEN
   GO TO  2014
 ELSE
   GO TO  2015
 END IF
 2013 j = -j
 2014 forfil = andf(exfiat(j+1),32767)
 RETURN
 2015 forfil = andf(fiat(j+1),32767)
 RETURN
END FUNCTION forfil
