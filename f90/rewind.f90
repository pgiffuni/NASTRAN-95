SUBROUTINE REWIND ( FILE )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 
 NAME   = FILE
 CALL dsgefl
! CALL DBMMGR FOR REWIND SO TO SET BUFFER ADDRESS CORRECTLY
 CALL dbmmgr ( 3 )
 nblock = fcb( 4, ifilex )
 IF ( iprvop == 0 ) GO TO 30
! IF FILE OPEN FOR WRITE, THEN INITIAL BUFFER AND BLOCK NUMBER
 ibase( indbas+3 ) = 1
 ibase( indbas+4 ) = 6
!WKBNB NCL93007 11/94
! SET THE COUNTER FOR NUMBER OF STRINGS AND TERMS TO ZERO
 fcb( 16, ifilex ) = 0
 fcb( 17, ifilex ) = 0
!WKBNE NCL93007 11/94
 30    indclr = indbas+5
 indcbp = indclr
 CALL dssdcb
 RETURN
END SUBROUTINE REWIND
