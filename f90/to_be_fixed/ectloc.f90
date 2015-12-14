SUBROUTINE ectloc(*,ect,buf,ielem)
!*****
! ECTLOC IS A SPECIAL PURPOSE VERSION OF SUBROUTINE LOCATE.  ITS
! PURPOSE IS TO PASS THE ECT FILE SEQUENTIALLY POSITIONING EACH LOGICAL
! RECORD AFTER THE 3-WORD HEADER AND PROVIDING A POINTER TO THE
! APPROPRIATE ENTRY IN THE ELEM TABLE IN /GPTA1/. PLOTEL
! ELEMENTS ARE IGNORED.
!     NOTE---THE ECT FILE MUST BE OPEN ON EACH CALL.
 
!  ARGUMENTS
 
!     ECT   ---INPUT ---EINO FILE NAME OF THE ECT
!     BUF   ---IN/OUT---ADDRESS OF A 3-WORD ARRAY INTO WHICH
!                       THE FIRST 3 WORDS OF THE RECORD ARE READ
!     IELEM ---OUTPUT---POINTER TO 1ST WORD OF ENTRY IN ELEM
!                       TABLE IN /GPTA1/
 
! NON-STANDARD RETURN---GIVEN WHEN EOF HIT. ECT IS CLOSED BEFORE RETURN.
!*****
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: ect
 INTEGER, INTENT(IN)                      :: buf(3)
 INTEGER, INTENT(OUT)                     :: ielem
 INTEGER :: elem, plotel
 
 COMMON/ gpta1 / nelem, last, incr, elem(1)
 
 DATA plotel/ 4HPLOT /
 
! READ A 3-WORD RECORD HEADER. IF NOT 3 WORDS, TRY NEXT RECORD
 
 10 CONTINUE
 
 CALL READ(*90,*10,ect,buf,3,0,nread)
 
! SEARCH FOR MATCH OF FIRST WORD OF RECORD WITH ECT-ID WORD IN /GPTA1/
! IF FOUND AND NOT PLOTEL, RETURN POINTER.
 
 DO  i=1,last,incr
   IF( buf(1) == elem(i+3) ) GO TO 30
 END DO
 25 CALL fwdrec(*90,ect)
 GO TO 10
 30 IF( elem(i) == plotel ) GO TO 25
 ielem = i
 RETURN
 
! EOF ENCOUNTERED--CLOSE FILE AND RETURN.
 
 90 CALL CLOSE( ect, 1 )
 ielem = 0
 RETURN 1
END SUBROUTINE ectloc
