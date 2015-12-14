SUBROUTINE tubed
     
!***
!  THE TUBE BEING SO SIMILAR TO THE ROD, WE ALTER THE EST FOR THE TUBE
!  SO THAT IT IS IDENTICAL TO THE ONE FOR THE ROD AND THEN CALL RODD
! DOUBLE PRECISION VERSION
! SINGLE AND DOUBLE PRECISION VERSIONS OF THIS ROUTINE ARE IDENTICAL
! APART FROM THE NAME AND THE CALL TO RODD (RODS)
!***
 
 
! EST( 1) - ELEMENT ID.
! EST( 2) - SCALAR INDEX NUMBER FOR GRID POINT A
! EST( 3) - SCALAR INDEX NUMBER FOR GRID POINT B
! EST( 4) - MATERIAL ID.
! EST( 5) - OUTSIDE DIAMETER
! EST( 6) - THICKNESS
! EST( 7) - NON-STRUCTURAL MASS
! EST( 8) - COOR. SYS. ID. FOR GRID POINT A
! EST( 9) - BASIC COORDINATES OF GRID POINT A
! EST(10) -                ...
! EST(11) -                ...
! EST(12) - COOR. SYS. ID. FOR GRID POINT B
! EST(13) - BASIC COORDINATES OF GRID POINT B
! EST(14) -               ...
! EST(15) -               ...
! EST(16) - ELEMENT TEMPERATURE
 
 COMMON   /emgest/  est(100)
 COMMON   /condas/  pi
 
! ----------------------------------------------------------------------
 
 temp = est(5) - est(6)
 a = temp * est(6) * pi
 fj = .25 * a * ( temp**2 + est(6)**2 )
 c = .5 * est(5)
 m = 18
 DO  i = 1,10
   m = m - 1
   est(m) = est(m-1)
 END DO
 est(5) = a
 est(6) = fj
 est(7) = c
 CALL rodd
 RETURN
END SUBROUTINE tubed
