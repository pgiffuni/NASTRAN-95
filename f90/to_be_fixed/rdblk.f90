SUBROUTINE rdblk ( *, FILE, ifirst, left )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: ifirst
 INTEGER, INTENT(OUT)                     :: left
 
 NAME = FILE
 CALL    dsgefl
!        PRINT *,' RDBLK,NAME,IFILEX,INDBAS=',NAME,IFILEX,INDBAS
!        WRITE(6,40646)(IBASE(INDBAS+K),K=1,8)
 40646   FORMAT(' BUFFER=',8(1X,z8))
 IF ( iprvop /= 0 ) CALL dsmsg( 4 )
 IF ( ifirst /= 0 ) GO TO 10
 CALL dsrdnb
 10      CALL dbmmgr( 9 )
!      WRITE(6,44771)(FCB(K,IFILEX),K=1,15)
 44771 FORMAT(' RDBLK,FCB=',/,2(5I8,/),2I8,4X,2A4,4X,i8)
!      INNN = FCB( 12, IFILEX )
!       PRINT *,' RDBLK-2,IFILEX,INNN=',IFILEX,INNN
!       WRITE(6,40646)(IBASE(INNN+K),K=0,7)
 ibase( indbas+1 ) = ibase( indbas+4 )
 ibase( indbas+2 ) = ibase( indbas+4 )
 left = nbuff + 3 - lcw
 IF ( ibase( indbas+lcw-2) == idsef ) RETURN 1
 RETURN
END SUBROUTINE rdblk
