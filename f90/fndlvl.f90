SUBROUTINE fndlvl (NAME,newnm)
     
!     THIS SUBROUTINE LOOKS FOR A LOWER LEVEL SUBSTRUCTUE TO THE
!     SUBSTRUCTURE NAME.  IF NAME DOES HAVE A LOWER LEVEL SUBSTRUCTURE,
!     THE NAME OF ONE OF THESE LOWER LEVEL SUBSTRUCTURES WILL BE
!     RETURNED IN NEWNM.  IF NAME DOES NOT HAVE A LOWER LEVEL
!     SUBSTRUCTURE, NAME WILL BE RETURNED IN NEWNM.  IF NAME IS NOT
!     KNOWN TO THE SYSTEM, BLANKS WILL BE RETURNED IN NEWNM.
 
 
 INTEGER, INTENT(IN)                      :: NAME(2)
 INTEGER, INTENT(OUT)                     :: newnm(2)
 EXTERNAL        rshift,andf
 INTEGER :: rshift,andf,buf
 DIMENSION  nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 DATA    ll    / 2 /
 DATA    iempty/ 4H     /, nmsbr / 4HFNDL,4HVL   /
 
!     CHECK IF NAME EXISTS
 
 CALL chkopn (nmsbr(1))
 CALL fdsub (NAME(1),k)
 IF(k /= -1) GO TO 10
 newnm(1) = iempty
 newnm(2) = iempty
 RETURN
 
!     FIND THE LOWER LEVEL SUBSTRUCTURE
 
 10 CALL fmdi (k,imdi)
 ill = andf(rshift(buf(imdi+ll),20),1023)
 IF(ill == 0) GO TO 20
 
!     NAME DOES HAVE A LOWER LEVEL SUBSTRUCTURE
 
 CALL fdit (ill,jdit)
 newnm(1) = buf(jdit)
 newnm(2) = buf(jdit+1)
 RETURN
 
!     NAME DOES NOT HAVE A LOWER LEVEL SUBSTRUCTURE
 
 20 newnm(1) = NAME(1)
 newnm(2) = NAME(2)
 RETURN
END SUBROUTINE fndlvl
