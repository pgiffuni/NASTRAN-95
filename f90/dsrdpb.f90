SUBROUTINE dsrdpb
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 nblock = nblock - 1
 CALL dbmmgr( 6 )
 indclr = ibase( indbas+4 ) + indbas - 1
 indcbp = indclr
 iblk   = ibase( indbas+3 )
 IF ( iblk == nblock ) GO TO 10
 CALL dsmsg( 102 )
 10      RETURN
END SUBROUTINE dsrdpb
