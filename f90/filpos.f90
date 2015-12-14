SUBROUTINE filpos ( FILE, ipos )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: ipos
 
 NAME  = FILE
 CALL dsgefl
 nblock = IAND( ipos, maskh2 )
 icblk  = fcb( 4, ifilex )
 IF ( icblk == nblock ) GO TO 10
 CALL dbmmgr( 6 )
 10      CONTINUE
 indclr = ipos/mulq2 + indbas - 1
 indcbp = indclr
 CALL dssdcb
 RETURN
END SUBROUTINE filpos
