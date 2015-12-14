SUBROUTINE fndnxl (NAME,newnm)
     
!     THE SUBROUTINE LOOKS FOR A HIGHER LEVEL SUBSTRUCTURE TO THE
!     SUBSTRUCTURE NAME.  IF NAME DOES HAVE A HIGHER LEVEL SUBSTRUCTURE,
!     THE NAME OF THE HIGHER LEVEL SUBSTRUCTURE WILL BE RETURNED IN
!     NEWNM.  IF NAME DOES NOT HAVE A HIGHER LEVEL SUBSTRUCTURE, NAME
!     WILL BE RETURNED IN NEWNM.  IF NAME IS NOT KNOWN TO THE SYSTEM,
!     BLANKS WILL BE RETURNED IN NEWNM.
 
 
 INTEGER, INTENT(IN)                      :: NAME(2)
 INTEGER, INTENT(OUT)                     :: newnm(2)
 EXTERNAL        andf
 LOGICAL :: ditup,mdiup
 INTEGER :: andf,buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     mdi,mdipbn,mdilbn,mdibl,hl
 DIMENSION  nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,iodum(8),  &
     mdi,mdipbn,mdilbn,mdibl,nxtdum(15),ditup,mdiup
 DATA    hl    / 2     /
 DATA    iempty/ 4H    /, nmsbr /4HFNDN,4HXL  /
 
 CALL chkopn (nmsbr(1))
 CALL fdsub  (NAME(1),k)
 IF (k /= -1) GO TO 10
 newnm(1) = iempty
 newnm(2) = iempty
 RETURN
 
!     FIND THE HIGHER LEVEL SUBSTRUCTURE TO NAME.
 
 10 CALL fmdi (k,imdi)
 i = andf(buf(imdi+hl),1023)
 IF (i == 0) GO TO 20
 
!     NAME DOES HAVE A HIGHER LEVEL SUBSTRUCTURE.
 
 CALL fdit (i,jdit)
 newnm(1) = buf(jdit  )
 newnm(2) = buf(jdit+1)
 RETURN
 
!     NAME DOES NOT HAVE A HIGHER LEVEL SUBSTRUCTURE.
 
 20 newnm(1) = NAME(1)
 newnm(2) = NAME(2)
 RETURN
END SUBROUTINE fndnxl
