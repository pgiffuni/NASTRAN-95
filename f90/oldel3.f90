SUBROUTINE oldel3
     
!     ANY ELEMENT (NEW OR OLD) WHICH HAS NOT BEEN CONVERTED TO USE
!     EMGPRO SHOULD HAVE AN ENTRY POINT IN OLDEL1, OLDEL2, OR OLDEL3
!     ***************************************************************
 
 ENTRY tetras
 GO TO 10
 ENTRY tetrad
 GO TO 10
 ENTRY traprs
 GO TO 10
 ENTRY traprd
 GO TO 10
 ENTRY triars
 GO TO 10
 ENTRY triard
 GO TO 10
 ENTRY tria1s
 GO TO 10
 ENTRY tria1d
 GO TO 10
 ENTRY tria2s
 GO TO 10
 ENTRY tria2d
 GO TO 10
 ENTRY trplts
 GO TO 10
 ENTRY trpltd
 GO TO 10
 ENTRY wedges
 GO TO 10
 ENTRY wedged
 
 10   CALL emgold
 RETURN
END SUBROUTINE oldel3
