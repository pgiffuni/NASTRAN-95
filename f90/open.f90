SUBROUTINE OPEN(*,namfil,buff,op)
!******
 
! OPEN IS AN INTERMEDIARY TO ENTRY POINT QOPEN IN SUBROUTINE GINO.
! THE MAIN TASK OF OPEN IS TO INSURE THAT DATA BLOCKS WHICH WERE
! WRITTEN AND CLOSED OFF THE LOAD POINT HAVE AN END-OF-FILE BEFORE
! THEY ARE READ.
 
!******

 INTEGER, INTENT(IN)                      :: namfil
 INTEGER, INTENT(IN OUT)                  :: buff(1)
 INTEGER, INTENT(IN)                      :: op
 INTEGER :: xop, xname
 COMMON /system/ isystm(157)
 INCLUDE 'DSIOF.COM'
 
 
! TEST FOR CONDITION IN WHICH END-OF-FILE IS TO BE WRITTEN
 
 DATA init / 0 /
 
 IF ( init /= 0 ) GO TO 5
 CALL dsiodd
 init = 1
 5  CONTINUE
 xname = namfil
 ifilex = 0
 CALL geturn( xname )
 IF(ifilex == 0)RETURN 1
 10 IF( op == 1 .OR. op == 3 ) GO TO 80
 IF( nblock+nlr > 7 ) GO TO 12
 11 IF( op == -2 ) RETURN
 GO TO 80
 12 IF( iprvop == 0 ) GO TO 11
 
! DATA BLOCK WAS PREVIOUSLY OPENED TO WRITE AND IS NOW OFF LOAD POINT.
! WRITE AN END-OF-FILE. IF SPECIAL CALL, RETURN
 
 CALL qopen(*88,namfil,buff,3)
 CALL eof( namfil )
 xop = 2
 IF( op == -2 ) xop = 1
 CALL CLOSE( namfil, xop )
 IF( op == -2 ) RETURN
 
! NOW OPEN ACCORDING TO OP. IF NECESSARY, POSITION PRIOR TO EOF
 
 lasnam = 0
 CALL geturn( namfil )
 CALL qopen(*88,namfil,buff,op)
 IF( op == 2 ) CALL bckrec( namfil )
 RETURN
 
! NORMAL OPEN CALL
 
 80 CALL qopen(*88,namfil,buff,op)
!WKBNB NCL93007 11/94
! SET THE COUNT FOR THE TOTAL NUMBER OF STRINGS AND TERMS
! TO ZERO IF FILE IS BEING OPENED FOR WRITE
 IF ( op /= 1 ) GO TO 70
 fcb( 16, ifilex ) = 0
 fcb( 17, ifilex ) = 0
 70 CONTINUE
!WKBNE NCL93007 11/94
 RETURN
 88 RETURN 1
END SUBROUTINE OPEN
