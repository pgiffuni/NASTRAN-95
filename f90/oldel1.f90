SUBROUTINE oldel1
     
!     ANY ELEMENT (NEW OR OLD) WHICH HAS NOT BEEN CONVERTED TO USE
!     EMGPRO SHOULD HAVE AN ENTRY POINT IN OLDEL1, OLDEL2, OR OLDEL3
!     ***************************************************************
 
 ENTRY axif2s
 GO TO 10
 ENTRY axif2d
 GO TO 10
 ENTRY axif3s
 GO TO 10
 ENTRY axif3d
 GO TO 10
 ENTRY axif4s
 GO TO 10
 ENTRY axif4d
 GO TO 10
 ENTRY cones
 GO TO 10
 ENTRY coned
 GO TO 10
 ENTRY elbows
 GO TO 10
 ENTRY elbowd
 GO TO 10
 ENTRY flmass
 GO TO 10
 ENTRY flmasd
 GO TO 10
 ENTRY flud2s
 GO TO 10
 ENTRY flud2d
 GO TO 10
 ENTRY flud3s
 GO TO 10
 ENTRY flud3d
 GO TO 10
 ENTRY flud4s
 GO TO 10
 ENTRY flud4d
 
 10   CALL emgold
 RETURN
END SUBROUTINE oldel1
