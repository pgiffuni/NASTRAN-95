SUBROUTINE dbmrel
!********************************************************************
!  DBMREL  -   RELEASES IN-MEMORY BLOCKS THAT ARE CURRENTLY
!              ALLOCATED TO AN IN-MEMORY FILE
!********************************************************************
 INCLUDE   'DSIOF.COM'
 INCLUDE   'ZZZZZZ.COM'
 COMMON / system / isysbf, iwr
 IF ( fcb( 9, ifilex ) == 0 .OR. fcb( 10, ifilex ) == 0 ) GO TO 701
 IF ( idbfre /= 0 ) GO TO 10
! FREE CHAIN IS EMPTY, THIS CHAIN BECOMES FREE CHAIN
 idbfre = fcb( 9, ifilex )
 GO TO 777
! SET FIRST OF BLOCKS TO BE FREED AT FIRST OF FREE CHAIN AND
! THEN CONNECT LAST OF BLOCKS TO BE FREED WITH FIRST OF EXISTING
! FREE CHAIN
 10    CONTINUE
 IF ( fcb( 9, ifilex ) == fcb( 10, ifilex ) ) GO TO 20
 isave          = idbfre
 idbfre         = fcb(  9, ifilex )
 mem( isave )   = fcb( 10, ifilex )
 INDEX          = fcb( 10, ifilex )
 mem( INDEX+1 ) = isave
 GO TO 777
! FILE HAD ONLY ONLY ONE BLOCK ALLOCATED TO IT
 20    CONTINUE
 isave          = idbfre
 idbfre         = fcb(  9, ifilex )
 mem( isave )   = idbfre
 mem( idbfre+1) = isave
 GO TO 777
 701   WRITE( iwr, 901 )
 901   FORMAT(///,' ERROR IN ATTEMPT TO FREE BLOCKS TO FREE CHAIN',  &
     /,' CONTENTS OF THE DIRECTORY ARE AS FOLLOWS')
 CALL dbmdmp
 CALL mesage ( -61, 0, 0 )
 777   CONTINUE
 fcb(  9, ifilex ) = 0
 fcb( 10, ifilex ) = 0
 fcb( 11, ifilex ) = 0
 RETURN
END SUBROUTINE dbmrel
