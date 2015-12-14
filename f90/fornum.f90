SUBROUTINE fornum ( FORM, ICHAR, imult )
     
! THIS SUBROUTINE CONVERTS ALL NUMBERS TO INTEGER FORMAT
 
 
 CHARACTER (LEN=1), INTENT(IN OUT)        :: FORM(200)
 INTEGER, INTENT(OUT)                     :: ICHAR
 INTEGER, INTENT(OUT)                     :: imult
 CHARACTER (LEN=1) :: BLANK, NUMBER(2)
 DATA          BLANK /' ' /
 DATA          NUMBER /'0','9'/
 
 imult = 0
 10    IF ( FORM( ICHAR ) /= BLANK      ) GO TO 20
 ICHAR = ICHAR + 1
 GO TO 10
 20    IF ( FORM( ICHAR ) < NUMBER(1) .OR.  &
     FORM( ICHAR ) > NUMBER(2) ) GO TO 700
 READ ( FORM( ICHAR ), 901 ) ii
 901   FORMAT(i1)
 imult = imult*10 + ii
 ICHAR = ICHAR + 1
 GO TO 20
 700   CONTINUE
 RETURN
END SUBROUTINE fornum
