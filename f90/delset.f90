SUBROUTINE delset
!*****
!  THIS ROUTINE SETS VARIABLES FOR DUMMY ELEMENTS IN /GPTA1/
 
!  ALL MODULES USING /GPTA1/ SHOULD BE SURE TO CALL THIS ROUTINE
!  SO AS TO INSURE THAT DATA FOR ANY DUMMY ELEMENTS PRESENT GETS
!  INSERTED INTO /GPTA1/.
!*****
 INTEGER :: dumtyp(9)
 
 COMMON/system/ nskip(45), idum(9)
 
 COMMON/gpta1 / nelem, last, incr, NE(1)
 
 DATA dumtyp/53,54,55,56,57,58,59,60,61/
 
 DO  i = 1,9
   ngrids = idum(i) / 10000000
   nc = MOD( idum(i),10000000) / 10000
   np = MOD( idum(i),10000) / 10
   
!     ND IS DECODE AND USED IN ROUTINES DS1 AND DS1A
   
   nd = MOD(idum(i), 10)
   izero = (dumtyp(i) - 1)*incr
   NE(izero+6) = nc + ngrids + 2
   NE(izero+9) = np +2
   IF(np == 0) NE(izero+9) =  0
   NE(izero+10) = ngrids
   n = 5*ngrids + 3 + np + nc
   NE(izero+12) = n
   NE(izero+15) = ngrids + 2
 END DO
 RETURN
END SUBROUTINE delset
