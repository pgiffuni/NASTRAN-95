SUBROUTINE dumod4
     
!*****
 
!     DUMMY DECK FOR MODULE DUMMOD4 - SEE USER'S MANUAL SECTION 5.6.
!                                     FOR MODULE PROPERTIES, CHECK
!                                     SUBROUTINE XMPLDD OR USE DIAG 31.
 
!*****
 
 COMPLEX :: parm9
 
 DOUBLE PRECISION :: parm8,parm10
 
 INTEGER :: parm1,parm2,parm3,parm4,parm7
 INTEGER :: infile(8),outfil(8),scrfil(10)
 
 COMMON /BLANK /  parm1,parm2,parm3,parm4,parm5,parm6,  &
     parm7(2),parm8,parm9,parm10(2)
 COMMON /zzzzzz/  x(1)
 
 DATA    infile/  101,102,103,104,105,106,107,108/
 DATA    outfil/  201,202,203,204,205,206,207,208/
 DATA    scrfil/  301,302,303,304,305,306,307,308,309,310/
 
 RETURN
END SUBROUTINE dumod4
