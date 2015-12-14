SUBROUTINE smsg (no,p1,p2)
     
!     MESSAGE WRITER FOR SUBSTRUCTURE DIAGNOSTICS, 61XX SERIES
 
 
 INTEGER, INTENT(IN OUT)                  :: no
 INTEGER, INTENT(IN OUT)                  :: p1
 INTEGER, INTENT(IN OUT)                  :: p2(2)
 INTEGER :: p3(2),pos(2),neg(2),png(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,k
 DATA    pos   , neg/4HWARN,4HING ,4HFATA,4HL   /
 DATA    nmsg  / 8  /, nmsg1 / 11 /
 
 l = IABS(no)
 msgno  = l + 6100
 IF (l < 1 .OR. l > nmsg) GO TO 99
 IF (no > 0) THEN
   GO TO    40
 END IF
 30 png(1) = neg(1)
 png(2) = neg(2)
 GO TO 50
 40 png(1) = pos(1)
 png(2) = pos(2)
 GO TO 50
 
 
 ENTRY smsg1 (no,p1,p2,p3)
!     ========================
 
 l = IABS(no)
 IF (l <= nmsg .OR. l > nmsg1) GO TO 99
 
 50 SELECT CASE ( l )
   CASE (    1)
     GO TO 1
   CASE (    2)
     GO TO 2
   CASE (    3)
     GO TO 3
   CASE (    4)
     GO TO 4
   CASE (    5)
     GO TO 5
   CASE (    6)
     GO TO 6
   CASE (    7)
     GO TO 7
   CASE (    8)
     GO TO 8
   CASE (    9)
     GO TO 9
   CASE (   10)
     GO TO 10
   CASE (   11)
     GO TO 11
   CASE (   12)
     GO TO 12
 END SELECT
 1 WRITE (k,201) png,msgno
 WRITE (k,101) p1,p2
 GO TO 100
 2 WRITE (k,201) png,msgno
 WRITE (k,102) p1,p2
 GO TO 100
 3 WRITE (k,201) png,msgno
 WRITE (k,103) p1,p2
 GO TO 100
 4 WRITE (k,201) png,msgno
 WRITE (k,104) p2
 GO TO 100
 5 WRITE (k,201) png,msgno
 WRITE (k,105) p2
 GO TO 100
 6 WRITE (k,200) png,msgno
 WRITE (k,106) p1,p2
 GO TO 100
 7 WRITE (k,200) png,msgno
 WRITE (k,107) p1,p2
 GO TO 100
 8 WRITE (k,200) png,msgno
 WRITE (k,108) p1,p2
 GO TO 100
 9 WRITE (k,109) p3,p1,p2
 GO TO 100
 10 WRITE (k,110) p3,p1,p2
 GO TO 100
 11 WRITE (k,111) p3,p1,p2
 GO TO 100
 12 WRITE (k,112)
 GO TO 100
 99 WRITE (k,120) no,p1,p2
 100 IF (no > 0) RETURN
 IF (l <= nmsg) CALL sofcls
 WRITE (k,130)
 CALL errtrc ('SMSG    ',130)
 RETURN
 
 101 FORMAT (' REQUESTED SOF ITEM DOES NOT EXIST.  ITEM ',a4,  &
     ', SUBSTRUCTURE ',2A4)
 102 FORMAT (' REQUESTED SUBSTRUCTURE DOES NOT EXIST.  ITEM ',a4,  &
     ', SUBSTRUCTURE ',2A4)
 103 FORMAT (' REQUESTED SOF ITEM HAS INVALID NAME.  ITEM ',a4,  &
     ', SUBSTRUCTURE ',2A4)
 104 FORMAT (' ATTEMPT TO CREATE DUPLICATE SUBSTRUCTURE NAME ',2A4)
 105 FORMAT (' ATTEMPT TO RE-USE SUBSTRUCTURE ',2A4,' IN A REDUCE ',  &
     ' OR COMBINE OPERATION.  USE EQUIV SUBSTRUCTURE COMMAND')
106 FORMAT (' UNEXPECTED END OF GROUP ENCOUNTERED WHILE READING ITEM '  &
    ,       a4,', SUBSTRUCTURE ',2A4)
107 FORMAT (' UNEXPECTED END OF ITEM ENCOUNTERED WHILE READING ITEM ',  &
    a4,', SUBSTRUCTURE ',2A4)
108 FORMAT (' INSUFFICIENT SPACE ON SOF FOR ITEM ',a4,', SUBSTRUCTURE'  &
    ,       1X,2A4)
109 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
    ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
110 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
    ' OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURGED.')
111 FORMAT (a23,' 6215, MODULE ',2A4,' - ITEM ',a4,  &
    ' OF SUBSTRUCTURE ',2A4,' PSEUDO-EXISTS ONLY.')
112 FORMAT (' ')
120 FORMAT (' NO MESSAGE FOR MESSAGE NO.',i5,' PARAMETERS = ',2I10, 10X,2A10)
130 FORMAT (//,' FATAL ERROR')
200 FORMAT (' *** SYSTEM ',2A4,' MESSAGE',i5)
201 FORMAT (' *** USER ',2A4,' MESSAGE',i5)
END SUBROUTINE smsg
