SUBROUTINE ffhelp (*,*,j)
     
 INTEGER, INTENT(IN OUT)                  :: j
 CHARACTER (LEN=1) :: qmark
 CHARACTER (LEN=4) :: STOP,     yes,      help,     xx
 COMMON /machin/ mach
 COMMON /system/ dummy(3), in
 COMMON /xreadx/ nout
 COMMON /xechox/ skip(2),  iechos
 COMMON /qmarkq/ qmark
 DATA    STOP,   yes,      help  /   'STOP',   'Y   ',   'HELP' /
 
!     THIS ROUTINE IS CALLED ONLY BY FF
 
 SELECT CASE ( j )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 50
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 120
   CASE (    5)
     GO TO 140
 END SELECT
 10   WRITE (nout,20)
 20   FORMAT (///1X,  &
 'GENERATED OUTPUT CARDS ARE SAVED ONLY IF FILE NAME IS GIVEN.',  &
     //,' YOU MAY ENTER NASTRAN EXECUTIVE CONTROL AND CASE CONTROL',  &
     ' CARDS FIRST',/,' (NO INPUT ECHO ON SCREEN)', //,  &
     ' ADDITIONAL INPUT INFORMATION WILL BE GIVEN WHEN YOU ENTER ',  &
     12H'BEGIN BULK', //,  &
     ' YOU MAY QUIT FREE-FIELD PROGRAM AT ANY TIME BY ENTERING ',  &
     6H'STOP', /,' NORMALLY, JOB TERMINATES BY ',9H'ENDDATA', //,  &
     ' YOU MAY USE ',10H'READFILE',' COMMAND TO READ ANY FILE WHICH',  &
     14H was 'STOPPED', /,  &
     ' BEFORE, AND CONTINUE FROM WHERE THE PREVIOUS JOB ENDED', //,  &
     ' FREE-FIELD INPUT IS AVAILABLE ONLY IN BULK DATA SECTION', /,  &
     ' AND IS ACTIVATED BY A COMMA OR EQUAL SIGN IN COLS. 1 THRU 10',  &
     //,' BOTH UPPER-CASE AND LOWER-CASE LETTERS ARE ACCEPTABLE',//,  &
     ' REFERENCE - G.CHAN: ',1H','cosmic/nastran free-field INPUT',  &
     2H',, /13X,'12TH NASTRAN USERS',1H',' colloquium, may 1984')
 WRITE (nout,30) qmark
 30   FORMAT (/,' MORE',a1,' (Y,N) - ')
 
 READ (in,40,END=80) xx
 40   FORMAT (a4)
 CALL upcase (xx,4)
 IF (xx /= yes) GO TO 80
 
 50   WRITE (nout,60)
 60   FORMAT (///,' THE FOLLOWING SYMBOLS ARE USED FOR FREE-FIELD INPUT'  &
     , //10X,'SYMBOL', 12X,'FUNCTION',/,9X,2('----'),5X,10('----'),  &
     /10X,', OR BLANK  FIELD SEPERATORS',  &
     /10X,'  =         DUPLICATES ONE CORRESPONDING FIELD',  &
     /10X,'  ==        DUPLICATES THE REMAINING FIELDS',  &
     /10X,'  *(N)      INCREMENT BY N', /10X,'  %(E)      ENDING VALUE BY E',  &
     /10X,'  /         THIS INPUT FIELD IS SAME AS PREVIOUS FIELD',  &
     /10X,'  J)        FIELD INDEX, J-TH FIELD (MUST FOLLOWED BY A V ALUE)',  &
     /10X,')+ OR 10)   INDEX FOR CONTINUATION FIELD',  &
     /10X,'  )         (IN COL. 1 ONLY) DUPLICATES THE CONTINUATION  &
     ID',/22X,'OF PREVIOUS CARD INTO FIELD 1 OF CURRENT CARD',  &
     /10X,'  ,         COL.1 ONLY, AUTO-CONTINUATION ID GENERATION',  &
     /10X,'  =(N)      1ST FIELD ONLY, DUPLICATES N CARDS WITH PROPE  &
     R',/22X,' INCREMENTS',  &
     /12X,'+A-I',6X,'CONTINUATION ID CAN BE DUPLICATED AUTOMATICALLY  &
 ', /22X,'ONLY IF IT IS IN PLUS-ALPHA-MINUS-INTEGER FORM',  &
     //1X,'EXAMPLES:',  /1X,'GRID, 101,,  0.  0. ,  7. 8)2  )+ABC-2',  &
     /1X,'=(11),*(1)  ,,  *(1.), /  %(23.45),==')
 IF (j == 1 .OR. iechos /= -2) GO TO 170
 WRITE (nout,30) qmark
 READ (in,40,END=80) xx
 CALL upcase (xx,4)
 IF (xx == yes) GO TO 140
 80   IF (xx == STOP) RETURN 2
 IF (mach == 4 .AND. in == 5) REWIND in
 GO TO 190
 
 100  WRITE (nout,110)
 110  FORMAT (//,24H enter 'N' for no punch,, /7X,  &
     38H'Y' for punch in free-field FORMAT, OR, /7X,  &
     43H'X' for punch in nastran fixed-field FORMAT,/)
 GO TO 190
 
 120  WRITE (nout,130)
 130  FORMAT (/,' MIFYLE - IS A RESERVED WORD.  TRY ANY OTHER NAME')
 GO TO 190
 
 140  WRITE (nout,150)
 150  FORMAT (//,' *** FREE-FIELD INPUT IS OPTIONAL.',//5X,'FOUR (4)',  &
     ' CONTROL OPTIONS ARE AVAILABLE - CAN BE ENTERED AT ANY TIME',  &
     /7X,'1.  PROMPT=ON, PROMPT=OFF, OR PROMPT=YES(DEFAULT)',  &
     /7X,'2.  SCALE/10,  OR SCALE/8',  &
     /7X,'3.  CANCEL=N,  (TO CANCEL N PREVIOUSLY GENERATED CARDS)',  &
     /7X,'4.  LIST  =N,  (TO   LIST N PREVIOUSLY GENERATED CARDS)',  &
 //7X,'ENTER ''HELP'' IF YOU NEED ADDITIONAL INSTRUCTIONS',  &
     //7X,'INTEGER INPUT SHOULD BE LIMITED TO 8 DIGITS',  &
     /7X,'UP TO 12 DIGITS ARE ALLOWED FOR FLOATING PT. NUMBER INPUT',  &
     /7X,'HOWEVER, ONLY UP TO 8 DIGIT ACCURACY IS KEPT',  &
     /7X,'             INPUT           RESULT ',  &
     /7X,'         ------------       --------',  &
     /7X,'E.G.     123.456789         123.4567',  &
     /7X,'         123.456789+6       .12345+9',  &
     /7X,'         -123.4567D+5       -.1234+8',  &
     /7X,'         123.45678E+4       1234567.',  &
     /7X,'         0.00123456-3       .12345-5',  &
     /7X,'         0.0123456789       .0123456',  &
     /7X,'         .00000123456       .12345-5')
 IF (iechos /= -2) WRITE (nout,160)
 160  FORMAT (/7X,'(3 AND 4 ARE AVAILABLE ONLY IN THE FREE-FIELD STAND',  &
     '-ALONE VERSION)')
 170  WRITE (nout,180)
 180  FORMAT (/4X,'UP TO 94 CHARATERS ALLOWABLE ON AN INPUT LINE. ',  &
     ' C/R TO CONTINUE')
 READ (in,40,END=80) xx
 CALL upcase (xx,4)
 IF (xx == help) GO TO 50
 IF (j /= 1) GO TO 80
 190  RETURN 1
END SUBROUTINE ffhelp
