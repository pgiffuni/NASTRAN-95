SUBROUTINE dsrdnb
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 CALL dbmmgr( 5 )
 nblock = fcb( 4, ifilex )
 indclr = indbas + 5
 indcbp = indclr
 lcw    = ibase( indbas+4 )
 iblk   = ibase( indbas+3 )
 IF ( iblk == nblock ) GO TO 10
 CALL dsmsg ( 102 )
 10      RETURN
END SUBROUTINE dsrdnb
