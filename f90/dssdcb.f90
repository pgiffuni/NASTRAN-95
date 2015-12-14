SUBROUTINE dssdcb
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 iclr  =  indclr - indbas + 1
 fcb( 3, ifilex )  = iclr
 fcb( 4, ifilex )  = nblock
 ibase( indbas+1 ) = indcbp - indbas + 1
 ibase( indbas+2 ) = iclr
 lasnam = NAME
!        WRITE(6,40646)IFILEX,NBLOCK,ICLR,INDBAS,INDCLR
 40646   FORMAT(' DSSDCB,IFILEX,NBLOCK,ICLR,BAS,CLR=',i3,i5,6I7)
 RETURN
END SUBROUTINE dssdcb
