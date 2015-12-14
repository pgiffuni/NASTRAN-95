SUBROUTINE dbmrlb ( INDEX )
!********************************************************************
!  DBMRLB  -   RELEASES AN IN-MEMORY BLOCK THAT IS CURRENTLY
!              ALLOCATED AS THE LAST BLOCK OF AN IN-MEMORY FILE.
!              THIS IS USED TO RELEASE THE NEXT ALLOCATED BLOCK FOR A
!              FILE OPENED FOR WRITE BUT WAS NEVER USED BECAUSE THE
!              FILE WAS CLOSED--I.E., THE LAST BLOCK ALLOCATED FOR A
!              FILE OPENED FOR WRITE IS NEVER USED BUT IT MUST HAVE
!              BEEN ALLOCATED JUST IN CASE THE FILE IS NOT TO BE CLOSED.
!********************************************************************
 INCLUDE   'DSIOF.COM'
 
 INTEGER, INTENT(IN OUT)                  :: INDEX
 COMMON / zzzzzz / mem(4)
 COMMON / system / isysbf, iwr
 
 indexl = INDEX
! CHECK IF OTHER BLOCKS ARE CHAINED TO THE END OF THIS BLOCK
 IF ( mem( INDEX+1 ) /= 0 ) GO TO 100
! SET "NEXT" OF PREVIOUS BLOCK TO ZERO, IF IT EXISTS
 5     lindex = mem( INDEX )
 IF ( lindex == 0 ) GO TO 10
 mem( lindex+1 ) = 0
 10    IF ( idbfre /= 0 ) GO TO 20
! FREE CHAIN IS EMPTY, THIS BLOCK BECOMES FREE CHAIN
 idbfre = INDEX
! SET "NEXT" AND "PREVIOUS" OF THIS CHAIN TO ZERO
 mem( INDEX    ) = 0
 mem( indexl+1 ) = 0
 GO TO 700
! SET BLOCKS TO BE FREED AT FIRST OF FREE CHAIN AND
! THEN CONNECT FREE CHAIN TO THIS BLOCK
 20    isave  = idbfre
 idbfre = INDEX
 mem( isave    ) = indexl
 mem( INDEX    ) = 0
 mem( indexl+1 ) = isave
 GO TO 700
! MORE THAN ONE BLOCK IN THIS CHAIN TO RELEASE BACK TO FREE CHAIN
 100   CONTINUE
 110   IF ( mem( indexl+1 ) == 0 ) GO TO 5
!WKBR SPR94012 10/94      INDEXL = MEM( INDEX+1 )
 indexl = mem( indexl+1 )
 GO TO 110
 700   CONTINUE
 RETURN
END SUBROUTINE dbmrlb
