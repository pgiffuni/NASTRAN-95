SUBROUTINE crsub (NAME,i)
     
    !     THE SUBROUTINE CREATES AN ENTRY FOR THE SUBSTRUCTURE NAME IN THE
    !     DIT THE OUTPUT PARAMETER I INDICATES THAT THE SUBSTRUCTURE NAME
    !     IS THE ITH SUBSTRUCTURE IN THE DIT.
 
 
    INTEGER, INTENT(IN)                      :: NAME(2)
    INTEGER, INTENT(OUT)                     :: i
    LOGICAL :: ditup
    INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl
    DIMENSION  iempty(2),nmsbr(2)
    COMMON /zzzzzz/ buf(1)
    COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
        iodum(8),mdidum(4),nxtdum(15),ditup
    DATA    iempty/ 2*4H    /
    DATA    indsbr/ 1 /, nmsbr /4HCRSU,4HB   /
 
    CALL chkopn (nmsbr(1))
    IF (ditsiz == ditnsb*2) GO TO 10
 
    !     THERE IS AN EMPTY INTERNAL DIRECTORY SPACE IN THE MDI.
 
    CALL fdsub (iempty(1),i)
    IF (i /= -1) GO TO 20
    GO TO 30
 
    !     NO INTERNAL EMPTY SPACE IN THE MDI.  DIRECTORY FOR THE NEW
    !     SUBSTRUCTURE
 
10  ditsiz = ditsiz + 2
    i = ditsiz/2
 
    !     UPDATE DIT.
 
20  ditnsb = ditnsb + 1
    CALL fdit (i,jdit)
    buf(jdit  ) = NAME(1)
    buf(jdit+1) = NAME(2)
    ditup = .true.
    RETURN
 
    !     ERROR MESSAGES.
 
30  CALL errmkn (indsbr,5)

    RETURN
END SUBROUTINE crsub
