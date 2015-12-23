SUBROUTINE page2 (lines)
     
!     2ND MASTER PAGING ROUTINE FOR NASTRAN
 
!     IABS(LINES) = NO. OF LINES TO BE ADDED FOR OUTPUT
!     IF CURRENT PAGE CAN NOT ACCOMODATE THE INCOMING LINES, A NEW PAGE
!     IS INITIATED WITH PROPER HEADINGS.
 
!     IF LINES IS NEGATIVE, A 6-LINE HEADER IS PRINTED.
!     IF LINES IS POSITIVE, A 3-LINE HEADER IS PRINTED AND FOLLOWED BY
!        3 BLANK LINES.
 
!     ENTRY POINT PAGE3 -
!     A 3-LINE HEADER IS PRINTED, NO BLANK LINES FOLLOWED. LINES CAN BE
!     NEGATIVE OR POSITIVE.
 
!     SIMPLIFIED BY G.CHAN/UNISYS, AND PAGE3 ADDED  12/92
 
 IMPLICIT INTEGER (a-z) 
 INTEGER, INTENT(IN OUT)                  :: lines
 INTEGER :: titlex(18),NAME(2),fchar,date,page
 CHARACTER (LEN=7) ::  machos
 CHARACTER (LEN=11) :: mchnam
 CHARACTER (LEN=30) :: ahead
 CHARACTER(LEN=3)   :: month(12)
 COMMON /chmach/ mchnam, machos
 COMMON /machin/ mach(4)
 COMMON /system/ ibuf,nout,dum6(6),sym,st,page,line,tline,maxlin,  &
                 date(3),dum15(15),ofp,dum8(8),crdate(3)
 COMMON /output/ title(32),subtit(32),label(32),head1(32), head2(32),head3(32)
 EQUIVALENCE     (titlex(1),title(1))
 DATA    month /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',  &
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/
 DATA    NAME  / 4H PAG, 4HE2  /
 
 flag  = 2
 
 10 IF (lines == 0) GO TO 100
 ll    = IABS(lines)
 IF (sym-line < ll .OR. ofp /= 0) GO TO 30
 20 line  = line + ll
 GO TO 100
 
 30 page  = page  + 1
 tline = tline + line
 line  = 0
 IF (tline > maxlin) GO TO 90
 in    = date(1)
 
!   ASSEMBLE PAGE HEADING
 
 ahead = ' '
 ncmnam = INDEX(mchnam,' ') - 1
 IF (ncmnam <= -1) ncmnam = 11
 ncmos  = INDEX(machos,' ') - 1
 IF (ncmos <= -1) ncmos = 7
 fchar = (18 - ncmnam - ncmos)/2 + 1
 WRITE (ahead(fchar:fchar+1),35) crdate(3)
 35 FORMAT (a2)
 fchar = fchar + 3
 ahead(fchar:30) = mchnam(1:ncmnam) // ' ' // machos(1:ncmos) // ' NASTRAN'
 WRITE  (nout,40) titlex, ahead, month(in),date(2),date(3),page
 40 FORMAT (1H1, 4X, 17A4, a2, ' /', a30, '/ ',a3, 1X, i2, ', ', i2, &
            ' / PAGE', i6)
 WRITE  (nout,50) subtit
 50 FORMAT ( 5X,31A4,a3)
 WRITE  (nout,60) label
 60 FORMAT (/5X,31A4,a3)
 line  = line + 4
 IF (flag  < 0) GO TO 20
 IF (lines > 0) GO TO 70
 
 WRITE (nout,60) (head1(i),i=1,32)
 WRITE (nout,50) (head2(i),i=1,32)
 WRITE (nout,50) (head3(i),i=1,32)
 line  = line + 4
 GO TO 20
 
 70 WRITE  (nout,80)
 80 FORMAT (///)
 line  = line + 4
 GO TO 20
 
!     MAX LINES EXCEEDED.  BUMP MAXLINES BY 3000 AND CALL MESAGE
 
 90 maxlin = maxlin + 3000
 CALL mesage (-19,tline,NAME)
 
 100 ofp  = 0
 RETURN
 
 
 ENTRY page3 (lines)
!     ===================
 
 flag = -3
 GO TO 10
 
END SUBROUTINE page2
