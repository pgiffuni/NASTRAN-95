SUBROUTINE decode (code,list,n)
     
!     DECODE DECODES THE BITS IN A WORD AND RETURNS A LIST OF INTEGERS
!     CORRESPONDING TO THE BIT POSITIONS WHICH ARE ON. NUMBERING
!     CONVENTION IS RIGHT (LOW ORDER) TO LEFT (HIGH ORDER) 00 THRU 31.
 
!     ARGUMENTS
 
!     CODE - INPUT  - THE WORD TO BE DECODED
!     LIST - OUTPUT - AN ARRAY OF DIMENSION .GE. 32 WHERE THE INTEGERS
!                     CORRESPONDING TO BIT POSITIONS ARE STORED
!     N    - OUTPUT - THE NUMBER OF ENTRIES IN THE LIST  I.E. THE NO.
!                     OF 1-BITS IN THE WORD
 
 
 
 INTEGER, INTENT(IN OUT)                  :: code
 INTEGER, INTENT(OUT)                     :: list(1)
 INTEGER, INTENT(OUT)                     :: n
 EXTERNAL     andf
 INTEGER :: andf,two
 COMMON /two/ two(32)
 
 n = 0
 DO  i = 1,32
   IF (andf(two(33-i),code) == 0) CYCLE
   n = n + 1
   list(n) = i - 1
 END DO
 
 RETURN
END SUBROUTINE decode
