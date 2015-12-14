SUBROUTINE pakcol(terms,nterms)
     
!     PACKS OUT A COLUMN OF AF OR DKGG MATRIX - DATA IS IN THE
!     FOLLOWING SAMPLE FORMATS.
 
!                  ---------------------
!                  I  NEGATIVE ROWSIL  I
!                  I-------------------I
!                  I  1 MATRIX TERM    I
!                  I-------------------I
!                  I  POSITIVE ROWSIL  I
!                  I-------------------I
!                  I                   I
!                  I  3 MATRIX TERMS   I
!                  I                   I
!                  ---------------------
 
!     MATRIX TERMS ARE IN DOUBLE PRECISION
 
 
 
 INTEGER, INTENT(IN OUT)                  :: terms(1)
 INTEGER, INTENT(IN OUT)                  :: nterms
 DOUBLE PRECISION :: val      ,tval
 
 INTEGER :: a        ,temp(7)
 
!     PACK COMMON BLOCK
 
 COMMON / zblpkx /       a(4)     ,irow
 
 EQUIVALENCE  ( val , a(1) )
 EQUIVALENCE  ( tval , a(3) )
 
!***********************************************************************
 
!     SORT THE MATRIX ENTRIES BY ABSOULUTE SIL VALUES
 
 iloc = 1
 10 isil = terms(iloc)
 jloc = iloc
 jsil = isil
 20 jloc = jloc + 3
 IF(jsil > 0) jloc = jloc + 4
 IF(jloc >= nterms) GO TO 60
 jsil = terms(jloc)
 IF(IABS(jsil) >= IABS(isil)) GO TO 20
 
 nt = 3
 IF(jsil > 0) nt = 7
 DO  i=1,nt
   temp(i) = terms(jloc+i-1)
 END DO
 
 kloc = jloc - 1
 DO  i=iloc,kloc
   j = kloc - i + iloc
   terms(j+nt) = terms(j)
 END DO
 
 DO  i=1,nt
   terms(iloc+i-1) = temp(i)
 END DO
 isil = jsil
 GO TO 20
 
 60 iloc = iloc + 3
 IF(isil > 0) iloc = iloc + 4
 IF(iloc < nterms) GO TO 10
 
!     PACK OUT TERMS - ADDING ANY IDENTICAL SIL
 
 iloc = 1
 70 irow = IABS(terms(iloc))
 nt = 2
 IF(terms(iloc) > 0) nt = 6
 
 DO  i=1,nt,2
   a(1) = terms(iloc+i)
   a(2) = terms(iloc+i+1)
   jloc = iloc
   80 j = jloc
   jloc = j + 3
   IF(terms(j) > 0) jloc = j + 7
   IF(jloc >= nterms) GO TO 90
   IF(terms(jloc) /= terms(iloc)) GO TO 90
   
!     DUPLICATE SILS - ADD THEM
   
   a(3) = terms(jloc+i)
   a(4) = terms(jloc+i+1)
   val = val + tval
   j = jloc
   GO TO 80
   
!     PACK OUT TERM
   
   90 CONTINUE
   CALL zblpki
   irow = irow + 1
 END DO
 
 iloc = jloc
 IF(iloc < nterms) GO TO 70
 
 RETURN
END SUBROUTINE pakcol
