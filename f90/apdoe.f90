SUBROUTINE apdoe(id,z,start,END,found,count)
     
!     APDOE FINDS AND OPEN ENDED CARD FOR ID
!     GIVEN A LIST Z(START ) TO Z(END)
!     FOUND = 0 IF NOT FOUND
!     FOUND = POINTER TO START OF CARD Z(FOUND)
!     COUNT = NUMBER OF DATA ITEMS NOT COUNTING THE ID
 
 
 INTEGER, INTENT(IN OUT)                  :: id
 INTEGER, INTENT(IN OUT)                  :: z(1)
 INTEGER, INTENT(IN)                      :: start
INTEGER, INTENT(IN)                      :: END
INTEGER, INTENT(OUT)                     :: found
INTEGER, INTENT(OUT)                     :: count

LOGICAL :: look

found = 0
look = .true.
count = 0
IF(start == 0) GO TO 50
DO  i = start,END
IF(look) GO TO 20
IF(z(i) == -1) look = .true.
CYCLE
20 IF(z(i) == id) GO TO 30
look = .false.
END DO
GO TO 50
30 found = i
j = i + 2
count = count + 1

!     START COUNT AT + 2 BECAUSE PAERO4 CARD CAN HAVE -1 IN FIELD 2

DO  i=j,END
IF(z(i) == -1) EXIT
count = count + 1
END DO
50 RETURN
END SUBROUTINE apdoe
