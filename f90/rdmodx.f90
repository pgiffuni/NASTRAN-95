SUBROUTINE rdmodx (FILE,mode,word)
     
!     ENTRY POINTS - RDMODX (FILE ,MODE,WORD)
!                    RDMODY (A    ,MODE,WORD)
!                    RDMODE (*,*,*,MODE,WORD)
!                    RDWORD (      MODE,WORD)
!     RDMODX, RDMODE AND RDWORD CALLED BY PLOT, FIND, PARAM AND SETINP
!     RDMODY CALLED ONLY BY PLOT
 
!     REVISED 10/10/92 BY G.CHAN/UNISYS
!     THE ORIGINAL WAY PASSING 'FILE' AND ARRAY 'A' FROM RDMODX AND
!     RDMODY ARE NOT ANSI FORTRAN77 STANDARD. THERE IS NO GUARANTY THAT
!     RDMODE AND RDWORD WILL PICK THEM UP CORRECTLY. MODIFICATIONS HERE
!     ARE (1) SAVE 'FILE' IN /XRDMOD/, AND (2) COMPUTE A REFERENCE
!     POINTER, REFPTR, SUCH THAT ARRAY A IS ACCESSIBLE VIA ARRAY Z
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: mode(1)
 INTEGER, INTENT(OUT)                     :: word(2)
 INTEGER :: filex,check1,check2,bitson,ENTRY,complf,eor,BLANK,  &
      refptr,z,a(1), NAME(2),next(2)
 COMMON /xrdmod/ filex,refptr,check1,check2,bitson,ENTRY
 COMMON /zzzzzz/ z(1)
 DATA    BLANK , eor,NAME / 1H ,1000000, 4HRDMO,4HDX  /
 
!     -RDMODX- IS CALLED IF -MODE- IS TO BE READ FROM DATA SET -FILE-
 
 ENTRY  = 0
 filex  = FILE
 check1 = 13579
 GO TO 10
 
 
 ENTRY rdmody (a,mode,word)
!     ==========================
 
!     -RDMODY- IS CALLED IF -MODE- IS TO BE READ FROM THE -A- ARRAY
 
!     COMPUTE THE REFERENCE POINTER FROM Z(1) TO A(1), AND NEXT TIME
!     WHEN A ARRAY IS USED, USE Z ARRAY WITH THE REFERENCE POINTER
 
 ENTRY  = 1
 refptr = locfx(a(1)) - locfx(z(1))
 check2 = 24680
 10 bitson = complf(0)
 RETURN
 
 
 ENTRY rdmode (*,*,*,mode,word)
!     ==============================
 
!     -RDMODE- IS CALLED TO READ -MODE-
!     IF MODE = -4, THE NEXT TWO WORDS ARE READ INTO -WORD-
!     IF MODE IS NEGATIVE AND NOT = -4, ONLY THE NEXT ONE WORD IS READ
!     INTO -WORD-
!     RETURN 1 - NUMERIC MODE (-MODE- NEGATIVE)
!                -MODE- = -1, -WORD- IS INTEGER
!                -MODE- = -2, -WORD- IS REAL NUMBER
!                -MODE- = -3, -WORD- IS ZERO ?
!                -MODE- = -4, -WORD- IS D.P.REAL
!     RETURN 2 - ALPHABETIC MODE (-MODE- POSITIVE)
!     RETURN 3 - END OF LOGICAL CARD (RECORD TERMINATED),
!                -MODE- = 1000000
 
 IF (ENTRY /= 0) GO TO 80
 IF (check1 /= 13579) CALL mesage (-37,0,NAME)
 
 20 CALL fread (filex,mode,1,0)
 IF (mode(1) < 0) THEN
   GO TO    70
 ELSE IF (mode(1) == 0) THEN
   GO TO    30
 ELSE
   GO TO    40
 END IF
 30 CALL fread (filex,0,0,1)
 GO TO 20
 40 IF (mode(1) >= eor) GO TO 60
 50 CALL fread (filex,next,2,0)
 IF (next(1) /= bitson .AND. next(1) /= BLANK) RETURN 2
 mode(1) = mode(1) - 1
 IF (mode(1) > 0) THEN
   GO TO    50
 ELSE
   GO TO    20
 END IF
 60 CALL fread (filex,0,0,1)
 RETURN 3
 
 70 i = 1
 IF (mode(1) == -4) i = 2
 CALL fread (filex,word,i,0)
 RETURN 1
 
 80 IF (check2 /= 24680) CALL mesage (-37,0,NAME)
 mode(1) = z(ENTRY+refptr)
 ENTRY   = ENTRY + 1
 IF (mode(1) < 0) THEN
   GO TO   120
 ELSE IF (mode(1) == 0) THEN
   GO TO    80
 END IF
 90 IF (mode(1) >= eor) GO TO 110
 100 next(1) = z(ENTRY+0+refptr)
 next(2) = z(ENTRY+1+refptr)
 ENTRY   = ENTRY + 2
 IF (next(1) /= bitson .AND. next(1) /= BLANK) RETURN 2
 mode(1) = mode(1) - 1
 IF (mode(1) > 0) THEN
   GO TO   100
 ELSE
   GO TO    80
 END IF
 110 ENTRY   = ENTRY + 1
 RETURN 3
 
 120 word(1) = z(ENTRY+refptr)
 ENTRY   = ENTRY + 1
 IF (mode(1) /= -4) RETURN 1
 word(2) = z(ENTRY+refptr)
 ENTRY   = ENTRY + 1
 RETURN 1
 
 
 ENTRY rdword (mode,word)
!     ========================
 
!     -RDWORD- IS CALLED TO READ TWO BCD WORDS INTO -WORD-
!     NOTE - ALL DATA FIELD DELIMITERS ARE SKIPPED
 
 word(1) = next(1)
 word(2) = next(2)
 130 mode(1) = mode(1) - 1
 IF (mode(1) <= 0) GO TO 160
 IF (ENTRY   /= 0) GO TO 140
 IF (check1  /= 13579) CALL mesage (-37,0,NAME)
 CALL fread (filex,next,2,0)
 GO TO 150
 
 140 IF (check2 /= 24680) CALL mesage (-37,0,NAME)
 next(1) = z(ENTRY  +refptr)
 next(2) = z(ENTRY+1+refptr)
 ENTRY   = ENTRY + 2
 150 IF (next(1) == bitson .OR. next(1) == BLANK) GO TO 130
 160 RETURN
END SUBROUTINE rdmodx
