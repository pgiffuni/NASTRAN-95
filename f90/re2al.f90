SUBROUTINE re2al (re,alph)
     
 
 REAL, INTENT(IN OUT)                     :: re
 INTEGER, INTENT(OUT)                     :: alph(2)
 EXTERNAL        lshift
 
 COMMON /system/ ibuf,nout,dummy(37),nbpw
 
 CALL fp2a8 (*40,re,alph)
 IF (nbpw-60 < 0) THEN
   GO TO    30
 ELSE IF (nbpw-60 == 0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 
!     FOR 60- OR 64- BIT MACHINES, SAVE THE SECOND HALF OF REAL NUMBER
!     IN THE SECOND ALPH WORD. THAT IS -
!     THE FULL REAL NUMBER IS IN ALPH(1), ALL 8 BYTES, OR
!     FIRST 4 BYTES IN ALPH(1), AND LAST 4 BYTES IN ALPH(2)
 
 10   alph(2) = lshift(alph(1),24)
 GO TO 30
 20   alph(2) = lshift(alph(1),32)
 30   RETURN
 
 40   WRITE  (nout,50)
 50   FORMAT (99X,'(IN FP2A8, CALLED FROM RE2AL)')
 CALL mesage (-61,0,0)
 GO TO 30
END SUBROUTINE re2al
