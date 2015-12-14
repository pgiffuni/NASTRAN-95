SUBROUTINE wrtblk ( FILE, iend )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: iend
 
 NAME = FILE
 CALL dsgefl
!        PRINT *,' WRTBLK,NAME,IFILEX,INDBAS=',NAME,IFILEX,INDBAS
!        WRITE(6,40646)(IBASE(INDBAS+K),K=0,7)
 40646   FORMAT(' WRTBLK,BUFFER=',8(1X,z8))
 CALL dbmmgr( 8 )
 indclr = ibase( indbas+4 )
 nblock = ibase( indbas+3 )
 ibase( indbas+1 ) = ibase( indbas+4 )
 ibase( indbas+2 ) = ibase( indbas+4 )
 fcb( 3,ifilex ) = indclr
 fcb( 4,ifilex ) = nblock
!       INNN = FCB(12, IFILEX)
!       PRINT *,' WRTBLK-2,IFILEX,INNN=',IFILEX,INNN
!       WRITE(6,40646)(IBASE(INNN+K),K=0,7)
 IF ( iend == 1 ) GO TO 700
 CALL dbmmgr( 4 )
 700     RETURN
END SUBROUTINE wrtblk
