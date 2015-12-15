SUBROUTINE xrgdev
     
!     PURPOSE - XRGDEV PROCESSES A FIELD FROM A ****CARD, ****FILE,
!               ****SBST, OR A ****RFMT CARD FROM THE RIGID FORMAT
!               DATA BASE
 
!     AUTHOR  - RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!      /SYSTEM/
!       NOUT    UNIT NUMBER FOR OUTPUT PRINT FILE
!      /XRGDXX/
!       ICOL    COLUMN CONTAINING THE FIRST CHARACTER OF THE FIELD
!       LIMIT   2 WORD ARRAY CONTAINING THE LOWER/UPPER LIMITS FOR
!               VALUES GIVEN IN THE FIELD
!       NUMBER  INTEGER VALUE FOR A ALPHA NUMBER WITHIN THE FIELD
!       RECORD  ARRAY IN 20A4 FORMAT CONTAINING THE CARD IMAGE
 
!     OUTPUT
!      /XRGDXX/
!       IERROR  ERROR FLAG IS NON-ZERO IF AN ERROR OCCURRED
!       NUM     2 WORD ARRAY CONTAINING THE VALUE(S) WITHIN THE CURRENT
!               FIELD
 
!     LOCAL VARIABLES
!       IND     INDEX TO THE ARRAY NUM
!       ISTATE  NEXT STATE (ROW = IN THE ABOVE DATA STATEMENT) TO BE
!               USED FOR SYNTAX VALIDATION BASED ON THE TYPE OF THE NEXT
!               CHARACTER IN THE FIELD
!       ISTR    COLUMN CONTAINING THE FIRST CHARACTER WITHIN THE FIEL
!       K       DO LOOP INDEX FOR SCANING CHARACTERS WITHIN THE FIELD
!       STATE   TABLE USED TO VALIDATE THE SYNTAX OF THE FIELD.  THE
!               NUMBER IN EACH ENTRY INDICATES THE ROW TO BE USED FOR
!               VALIDATING THE SYNTAX OF THE NEXT CHARACTER.  IF THE
!               VALUE IS 0 THEN A SYNTAX ERROR OCCURRED.
 
!     FUNCTIONS
!     XRGDEV SCANS THE FIELD FOR SYNTAX ERRORS AND FOR PLACING THE NUMBE
!     INTO THE NUM ARRAY.  VALID FIELDS ARE OF THE FORM 'NNN,' OR
!     'NNN-NNN,' WITH EMBEDDED BLANKS ALLOWED AND NUMBERS MAY BE OF
!     ANY VALUE THAT IS WITHIN THE LIMITS OF THE ARRAY LIMIT.
 
!     SUBROUTINES CALLED - XRGDTP
 
!     CALLING SUBROUTINES - XRGSUB,XRGDCF
 
!     ERRORS
!       ERROR MESSAGES 8021 AND 8022 ARE GIVEN FOR SYNTAX OR VALUE RANGE
!       ERRORS.
 
 INTEGER :: record, state(5,7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /xrgdxx/ irestr, nsubst, iphase, icol  , NUMBER, itype ,  &
     istate, ierror, num(2), ind   , nument        ,  &
     record(20)    , ICHAR(80)     , limit(2)      ,  &
     icount, idmap , iscr  , NAME(2),member(2)     , ignore
 COMMON /system/ isysbf, nout  , dum(98)
!                   NUMBER  ,      -    BLANK    OTHER
 DATA    state / 1,   2,      3,      6,      0,  &
     1,   0,      0,      2,      0, 4,   0,      0,      3,      0,  &
     4,   2,      0,      5,      0, 0,   2,      0,      5,      0,  &
     0,   2,      3,      6,      0, 1,   0,      0,      7,      0 /
 
 IF (icol > 80) GO TO 110
 istate = 7
 ind    = 1
 num(1) = 0
 istr   = icol
 DO  k = istr,80
   icol = k
   CALL xrgdtp
   istate = state(itype,istate)
   IF (istate /= 0) GO TO 20
   ierror = 1
   j = 0
   WRITE  (nout,10) ufm,k,record,j,(i,i=1,8),ierror,(j,i=1,8)
   10   FORMAT (a23,' 8020, SYNTAX ERROR NEAR COLUMN ',i3,  &
       ' IN THE FOLLOWING CARD- ',/20X,20A4, /,(20X,i1,i9,7I10))
   GO TO 110
   20   SELECT CASE ( istate )
     CASE (    1)
       GO TO 30
     CASE (    2)
       GO TO 60
     CASE (    3)
       GO TO 40
     CASE (    4)
       GO TO 30
     CASE (    5)
       GO TO 50
     CASE (    6)
       GO TO 50
     CASE (    7)
       GO TO 50
   END SELECT
   30   num(ind) = num(ind)*10 + NUMBER
   CYCLE
   40   ind    = 2
   num(2) = 0
   50   CONTINUE
 END DO
 60   IF (ind == 2) GO TO 70
 num(2) = num(1)
 GO TO 90
 70   IF (num(2) > num(1)) GO TO 90
 ierror = 1
 WRITE  (nout,80) ufm,num(1),num(2),record
 80   FORMAT (a23,' 8021, NON-INCREASING RANGE ',i3,1H-,i3,  &
     ' IN THE FOLLOWING CARD -', /20X,20A4)
 90   CONTINUE
 IF (num(1) >= limit(1) .AND. num(2) <= limit(2)) GO TO 110
 WRITE  (nout,100) ufm,limit,record
 100  FORMAT (a23,' 8022, NUMBERS ARE OUT OF THE RANGE ',i3,1H-,i3,  &
     ' IN THE FOLLOWING CARD - ', /20X,20A4)
 ierror = 1
 110  CONTINUE
 RETURN
END SUBROUTINE xrgdev
