SUBROUTINE dsspos( FILE, kcblk, kclr, kcbp )
     
! DSSPOS REPOSITIONS THE "FILE" TO BLOCK "KCBLK" WITH THE CURRENT
! LOGICAL RECORD POINTER SET TO "KCLR" AND THE CURRENT BUFFER
! POINTER SET TO "KCBP"
 
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: kcblk
 INTEGER, INTENT(IN)                      :: kclr
 INTEGER, INTENT(IN)                      :: kcbp
 
 NAME  = FILE
 CALL dsgefl
 icblk  = fcb( 4, ifilex )
 IF ( icblk == kcblk ) GO TO 10
 nblock = kcblk
 CALL dbmmgr( 6 )
 10    CONTINUE
 indclr = kclr + indbas - 1
 indcbp = kcbp + indbas - 1
 CALL dssdcb
 RETURN
END SUBROUTINE dsspos
