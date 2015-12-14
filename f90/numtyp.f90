FUNCTION numtyp ( ivalue )
     
 
 INTEGER, INTENT(IN)                      :: ivalue
 CHARACTER (LEN=2) :: byte(4)
 CHARACTER (LEN=8) :: word
 
 
 EQUIVALENCE    ( byte, word )
 
!      WRITE(6,40646) IVALUE
 40646 FORMAT(' NUMTYP,IVALUE=',z9)
 IF ( ivalue == 0 ) GO TO 200
 WRITE ( word, 2000 ) ivalue
 IF ( byte(1) == '  ' ) GO TO 210
 IF ( byte(1) == '00' ) GO TO 210
 IF ((byte(1) == '07'.OR. byte(1) == ' 7') .AND. byte(2) == 'ff' .AND.  &
     byte(3) == 'ff' .AND. byte(4) == 'ff' )  GO TO 210
 IF ( byte(1) == '7F' .AND. byte(2) == 'ff' .AND.  &
     byte(3) == 'ff' .AND. byte(4) == 'ff' )  GO TO 210
 
 IF ( byte(1) == 'ff' ) GO TO 210
 DO  i = 1, 4
   IF ( byte(i) < '1F' .OR. byte(i) > '5E' ) GO TO 220
 END DO
 GO TO 230
 
!     VALUE IS ZERO
 
 200   numtyp = 0
 GO TO 700
 
!     VALUE IS INTEGER
 
 210   numtyp = 1
 GO TO 700
 
!     VALUE IS REAL
 
 220   numtyp = 2
 GO TO 700
 
!     VALUE IS ALPHA
 
 230   numtyp = 3
 
 700   CONTINUE
 RETURN
!*****
 2000  FORMAT(z8)
!*****
END FUNCTION numtyp
