SUBROUTINE fname (FILE,NAME)
!*******
!     GIVEN A FILE NO., FNAME WILL RETURN THE BCD DESCRIPTOR
!*******
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: NAME(2)
 INTEGER :: fiat, fist
 COMMON /xfist / fist(2)
 COMMON /xfiat / fiat(1)
 DATA    nblank/ 4H    /
 DATA    non1  , non2  / 4H (NO,4HNE) /
!*******
!     SEARCH THE FIST FOR THE FILE
!*******
 n = fist(2)*2 + 2
 DO  j=3,n,2
   IF (FILE == fist(j)) GO TO 20
 END DO
!*******
!     FILE DOES NOT EXIST, RETURN -(NONE)-
!*******
 NAME(1) = non1
 NAME(2) = non2
 RETURN
 20 k = fist(j+1)
 IF (k > 0) THEN
   GO TO    30
 END IF
 21 CONTINUE
!*******
!     RETURN BCD DESCRIPTOR
!*******
 NAME(1) = FILE
 NAME(2) = nblank
 RETURN
 
 30 ix = fist(j+1) + 2
 NAME(1) = fiat(ix  )
 NAME(2) = fiat(ix+1)
 
 RETURN
END SUBROUTINE fname
