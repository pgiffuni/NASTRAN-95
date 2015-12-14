SUBROUTINE dsnmrd (iunitu)
     INCLUDE  'DSIOF.COM'
 
 INTEGER, INTENT(IN OUT)                  :: iunitu
 COMMON /gndate/ igndat(2)
 READ (iunitu) fcb, numdev, dev, mdsnam, igndat, maxblk, maxdsk
 READ (iunitu) numopn, numcls, numwri, numrea
 
 DO  i = 1, 80
   fcb( 9, i ) = 0
   fcb(10, i ) = 0
   fcb(11, i ) = 0
 END DO
 RETURN
END SUBROUTINE dsnmrd
