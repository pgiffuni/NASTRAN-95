SUBROUTINE filswi (name1,name2)
     
!     FILSWI SWITCHES THE UNITS ASSIGNED TO THE SPECIFIED DATA BLOCKS.
 
 
 INTEGER, INTENT(IN OUT)                  :: name1
 INTEGER, INTENT(IN OUT)                  :: name2
 EXTERNAL        complf,andf,orf
 INTEGER :: fist,fiat,complf,sys,andf,unit1,unit2,orf,UNIT
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /xfist / nfist,lfist,fist(1) /xfiat / fiat(4)  &
     /system/ sys,nout,skip(21),icfiat
 DATA    mask1 / 32767/
!                     '7FFF'X
 
!     SEARCH FIST FOR POINTERS TO FIAT.
 
 IF (name1 == name2) RETURN
 k1 = 0
 k2 = 0
 n  = 2*lfist - 1
 DO  i = 1,n,2
   IF (fist(i) == name1) k1 = fist(i+1)
   IF (fist(i) == name2) k2 = fist(i+1)
 END DO
 IF (k1 > 0 .AND. k2 > 0) GO TO 10
 WRITE  (nout,9) sfm
 9 FORMAT (a23,' 2178, GINO REFERENCE NAMES, IMPROPER FOR ',  &
     'SUBROUTINE FILSWI.')
 CALL mesage (-61,0,0)
 
!     SWITCH UNIT REFERENCE NUMBERS IN FIAT.
 
 10 mask2 = complf(mask1)
 unit1 = andf(fiat(k1+1),mask1)
 unit2 = andf(fiat(k2+1),mask1)
 n     = icfiat*fiat(3) - 2
 DO  i = 4,n,icfiat
   UNIT  = andf(fiat(i),mask1)
   IF (UNIT == unit1) fiat(i) = orf(andf(fiat(i),mask2),unit2)
   IF (UNIT == unit2) fiat(i) = orf(andf(fiat(i),mask2),unit1)
 END DO
 RETURN
END SUBROUTINE filswi
