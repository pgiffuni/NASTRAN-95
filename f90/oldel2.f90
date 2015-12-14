SUBROUTINE oldel2
     
!     ANY ELEMENT (NEW OR OLD) WHICH HAS NOT BEEN CONVERTED TO USE
!     EMGPRO SHOULD HAVE AN ENTRY POINT IN OLDEL1, OLDEL2, OR OLDEL3
!     ***************************************************************
 
 ENTRY hexa1s
 GO TO 10
 ENTRY hexa1d
 GO TO 10
 ENTRY hexa2s
 GO TO 10
 ENTRY hexa2d
 GO TO 10
 ENTRY plotls
 GO TO 10
 ENTRY plotld
 GO TO 10
 ENTRY qdmems
 GO TO 10
 ENTRY qdmemd
 GO TO 10
 ENTRY qdplts
 GO TO 10
 ENTRY qdpltd
 GO TO 10
 ENTRY quad1s
 GO TO 10
 ENTRY quad1d
 GO TO 10
 ENTRY quad2s
 GO TO 10
 ENTRY quad2d
 GO TO 10
 ENTRY slot3s
 GO TO 10
 ENTRY slot3d
 GO TO 10
 ENTRY slot4s
 GO TO 10
 ENTRY slot4d
 
 10   CALL emgold
 RETURN
END SUBROUTINE oldel2
