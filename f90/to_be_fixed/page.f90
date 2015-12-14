SUBROUTINE page
     
!     MASTER PAGING ROUTINE FOR NASTRAN.
 
 INTEGER :: otpe,date,crdate,sym,titlex(18),NAME(2),fchar
 CHARACTER (LEN=7) :: amachos
 CHARACTER (LEN=11) :: mchnam
 CHARACTER (LEN=30) :: head
 COMMON /chmach/ mchnam, machos
 COMMON /machin/ mach(4)
 COMMON /system/ sysbuf,otpe,mpcn(3),spcn,method,loadn,sym,st,  &
     ipage,line,itline,maxlin,date(3),dum15(15),iofp, x(8),crdate(3)
 COMMON /output/ title(32),subtit(32),label(32),head1(32), head2(32),head3(32)
 EQUIVALENCE     (titlex(1),title(1))
 DATA    month /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',  &
     'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/
 DATA    NAME  / 4HPAGE, 4H    /
 
 iout  = 1
 10 ipage = ipage  + 1
 itline= itline + line
 line  = 0
 IF (itline > maxlin) GO TO 70
 in = date(1)
 
!   ASSEMBLE PAGE HEADING
 
 ahead = ' '
 ncmnam = INDEX(mchnam,' ') - 1
 IF (ncmnam <= -1) ncmnam = 11
 ncmos  = INDEX(machos,' ') - 1
 IF (ncmos <= -1) ncmos = 7
 fchar = (18 - ncmnam - ncmos)/2 + 1
 WRITE (ahead(fchar:fchar+1),15) crdate(3)
 15 FORMAT (a2)
 fchar = fchar + 3
 ahead(fchar:30) = mchnam(1:ncmnam) // ' ' // machos(1:ncmos) // ' NASTRAN'
 
 WRITE  (otpe,20) titlex, ahead, month(in),date(2),date(3),ipage
 20 FORMAT (1H1,4X,17A4,a2,' /',a30,'/ ',a3,1X,i2,', ',i2, ' / PAGE',i6)
 WRITE  (otpe,30) subtit
 30 FORMAT (5X,31A4,a3)
 WRITE  (otpe,40) label
 40 FORMAT (1H0,4X,31A4,a3)
 line = line + 4
 IF (iout == 0) GO TO 60
 WRITE (otpe,40) (head1(i),i=1,32)
 WRITE (otpe,30) (head2(i),i=1,32)
 WRITE (otpe,30) (head3(i),i=1,32)
 line = line + 4
 60 RETURN
 
!     MAX LINES EXCEEDED.  BUMP MAXLINES BY 3000 AND CALL MESAGE
 
 70 maxlin = maxlin + 3000
 CALL mesage (-19,itline,NAME)
 GO TO 60
 
 
 ENTRY page1
!     ===========
 
 iout = 0
 GO TO 10
END SUBROUTINE page
