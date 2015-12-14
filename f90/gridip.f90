SUBROUTINE gridip (grid,seqss,LEN,ipset,cset,no,z,lloc)
     
!     THIS SUBROUTINE FINDS SETS OF IP NUMBERS AND DEGREE OF FREEDOM
!     COMPONENT NUMBERS FOR GRID POINTS DEFINED IN A BASIC
!     SUBSTRUCTURE THAT IS A COMPONENT OF A PSEUDO-STRUCTURE.
 
!     ARGUMENTS
!               GRID   - GRID POINT ID NUMBER
!               SEQSS  - THE STARTING ADDRESS IN OPEN CORE OF THE
!                        PSEUDO-STRUCTURE EQSS RECORD
!               LEN    - LENGTH OF THE EQSS
!               IPSET  - THE SET OF IP NUMBERS FOR GRID
!               CSET   - COMPONENTS OF GIVEN IP NUMBER
!               NO     - THE NUMBER OF IP DEFINED BY GRID
 
 
 
 INTEGER, INTENT(IN OUT)                  :: grid
 INTEGER, INTENT(IN)                      :: seqss
 INTEGER, INTENT(IN)                      :: LEN
 INTEGER, INTENT(OUT)                     :: ipset(6)
 INTEGER, INTENT(OUT)                     :: cset(6)
 INTEGER, INTENT(OUT)                     :: no
 INTEGER, INTENT(IN)                      :: z(1)
 INTEGER, INTENT(OUT)                     :: lloc
 EXTERNAL orf,rshift
 INTEGER :: orf,rshift, posno
 COMMON  /cmbfnd/ inam(2),ierr
 
 ierr = 0
 nent = LEN/3
 
!     SEARCH FOR THE GRID ID IN THE EQSS
 
!     NOTE --- FOR RAPID LOCATION OF ALL IP FOR A GIVEN GRID,
!              THE COMPONENT WORD OF THE EQSS HAS HAD ITS FIRST
!              SIX BITS PACKED WITH A CODE-  THE FIRST THREE
!              BITS GIVE THE NUMBER OF THE IP AND THE SECOND
!              THREE THE TOTAL NO. OF IP.  E.G. 011101 MEANS
!              THE CURRENT IP IS THE THIRD OF FIVE FOR THIS
!              GRID ID.
 
 
 CALL bisloc (*30,grid,z(seqss),3,nent,loc)
 k = seqss + loc - 1
 icode = rshift(z(k+2),26)
 
!     ICODE CONTAINS SIX BIT CODE
 
 posno = icode/8
 noapp = icode - 8*posno
 
!     POSNO IS THE POSITION NUMBER OF THE GRID WE HAVE FOUND,
!     NOAPP IS THE TOTAL NUMBER OF APPEARANCES OF THAT GRID.
 
 IF (noapp == 0) posno = 1
 IF (noapp == 0) noapp = 1
 istart = k - 3*(posno-1)
 lloc = istart
 
!     PICK UP RIGHT 26 BITS BY MASK26 FOR CSET(I), INSTEAD OF R/LSHIFT
 
 mask26 = maskn(26,0)
 
 DO  i = 1,noapp
   kk = istart + 3*(i-1)
   ipset(i) = z(kk+1)
   cset(i)  = orf(z(kk+2),mask26)
 END DO
 
 no = noapp
 GO TO 40
 30 ierr = 1
 40 RETURN
END SUBROUTINE gridip
