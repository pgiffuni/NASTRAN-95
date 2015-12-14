SUBROUTINE vaxsch (nin,nout)
     
!     TO SEARCH UNIT NIN FOR END OF BULK DATA DECK
 
 
 INTEGER, INTENT(IN OUT)                  :: nin
 INTEGER, INTENT(IN OUT)                  :: nout
 CHARACTER (LEN=8) :: e1,e2,e3,chr
DATA        e1,e2,e3 / 'ENDDATA ', 'END DATA', 'ENDATA  ' /


60 READ (nin,70,END=80) chr
70 FORMAT (a8)
IF (chr == e1 .OR. chr == e2 .OR. chr == e3) GO TO 100
GO TO 60

!     ENDDATA CARD NOT FOUND

80 WRITE  (nout,90)
90 FORMAT ('0*** USER FATAL MESSAGE: "ENDDATA" CARD NOT FOUND BY ',  &
    'INPUT MODULE')
CALL vaxend

100 RETURN
END SUBROUTINE vaxsch
