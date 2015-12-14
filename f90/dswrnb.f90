SUBROUTINE dswrnb
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 ibase( indbas+4 ) = indclr - indbas + 1
 CALL dbmmgr( 4 )
 nblock  = fcb( 4, ifilex )
 indclr  = indbas + 5
 ibase( indbas+3 ) = nblock
 indcbp  = indclr
 RETURN
END SUBROUTINE dswrnb
