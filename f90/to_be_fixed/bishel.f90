SUBROUTINE bishel(*,list,nent,nterm,array)
!-----
!   BISHEL IS A MERGE/SORT/DUPLICATE ENTRY ELIMINATOR.  GIVEN A SORTED
! -ARRAY- AND A -LIST- TO MERGE, BISHEL ADDS THE -LIST- IN THE SORTED
! LOCATION.  SORT IS ONLY ON THE FIRST WORD OF LIST.
 
!   ARGUMENTS...
 
!     LIST  -- IN/OUT - LIST OF LENGTH NTERM TO MERGE INTO ARRAY.
!     NENT  -- IN/OUT - LENGTH OF LIST BEFORE/AFTER MERGE.
!     NTERM -- IN     - LENGTH OF ARRAY (AND LIST) ENTRIES.
!     ARRAY -- IN/OUT - ARRAY TO MERGE LIST INTO.
!     NONSTANDARD RETURN -- WHEN ARRAY(ITERM) IS A DUPLICATE.
!-----
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: list(1)
 INTEGER, INTENT(IN OUT)                  :: nent
 INTEGER, INTENT(IN)                      :: nterm
 INTEGER, INTENT(IN OUT)                  :: array(1)
 
 
 k = 1
 l = nent + 1
 m = nent - nterm + 1
 
!  . LOCATE DUPLICATES...
 
 IF (nent < nterm) GO TO 50
 IF (list(1) - array(m)  < 0) THEN
   GO TO    10
 ELSE IF (list(1) - array(m)  == 0) THEN
   GO TO    20
 ELSE
   GO TO    60
 END IF
 10 kid = list(1)
 CALL bisloc (*30, kid, array, nterm, nent/nterm, k)
 20 RETURN 1
 
!  . CREATE A HOLE IN THE LIST BY MOVING THE END OF THE LIST...
 
 30 CONTINUE
 j = l-k
 n = nent+nterm
 DO  i = 1,j
   m = l-i
   array(n) = array(m)
   n = n-1
 END DO
 
!  . LOAD LIST INTO HOLE...
 
 GO TO 70
 50 nent = 0
 GO TO 70
 60 k = l
 70 CONTINUE
 DO  i = 1,nterm
   array(k)=list(i)
   k = k+1
 END DO
 nent = nent + nterm
 RETURN
END SUBROUTINE bishel
