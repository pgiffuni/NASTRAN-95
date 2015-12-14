SUBROUTINE edit (NAME,iopt,itest)
     
!     REMOVES SELECTED ITEMS OF THE SUBSTRUCTURE NAME FROM THE SOF.
!     THE VALUE OF IOPT IS THE SUM OF THE FOLLOWING INTEGERS REFLECTING
!     WHICH ITEMS ARE TO BE REMOVED.
 
!              1 = STIFFNESS MATRIX
!              2 = MASS MATRIX
!              4 = LOAD DATA
!              8 = SOLUTION DATA
!             16 = TRANSFORMATION DATA
!             32 = ALL ITEMS OF SUBSTRUCTURE
!             64 = APPENDED LOADS DATA
!            128 = DAMPING MATRICES
!            256 = MODES DATA
 
!     THE OUTPUT VARIABLE ITEST TAKES ON ONE OF THE FOLLOWING VALUES
!              1   NORMATL RETURN
!              4   IF NAME DOES NOT EXIST
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        andf
 INTEGER :: andf, nmsbr(2)
 COMMON /itemdt/ nitem,item(7,1)
 DATA    nmsbr / 4HEDIT,4H    /
 
 CALL chkopn (nmsbr(1))
 itest = 1
 IF (iopt <= 0) GO TO 20
 CALL fdsub (NAME(1),INDEX)
 IF (INDEX == -1) GO TO 30
 
!     REMOVE SELECTED ITEMS ACCORDING TO IOPT S VALUE.
 
 DO  i = 1,nitem
   mask = item(7,i)
   IF (andf(iopt,mask) /= 0) CALL DELETE (NAME,item(1,i),it)
 END DO
 20 RETURN
 
!     NAME DOES NOT EXIST.
 
 30 itest = 4
 RETURN
END SUBROUTINE edit
