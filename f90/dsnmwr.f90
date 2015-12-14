SUBROUTINE dsnmwr (iunitu)
     
!         THIS SUBROUTINE IS CALLED BY ENDSYS
 
 INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN OUT)                  :: iunitu
 COMMON /gndate/ igndat(2)
 WRITE (iunitu) fcb, numdev, dev, mdsnam, igndat,maxblk,maxdsk
 WRITE (iunitu) numopn, numcls, numwri, numrea
 RETURN
 
END SUBROUTINE dsnmwr
