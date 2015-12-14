SUBROUTINE xygraf (graph)
     
 
 REAL, INTENT(OUT)                        :: graph(3,8)
 LOGICAL :: exceed
 INTEGER :: z,titlec,xtitle,titlel,titler,sysbuf,m(5), nu(5,10)
 
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /system/ sysbuf,l,idum(6),nlpp,idum2(2),lines,itlns
 COMMON /xypppp/ iframe,titlec(32),titlel(14),titler(14),  &
     xtitle(32),id(300),maxplt,xmin,xinc,exceed,i123, maxrow
 COMMON /zzzzzz/ z(1)
 
 DATA    nu(1, 1),nu(1, 2),nu(1, 3) /  4H****   ,4H **   ,4H**** /,  &
     nu(2, 1),nu(2, 2),nu(2, 3) /  4H*  *   ,4H  *   ,4H   * /,  &
     nu(3, 1),nu(3, 2),nu(3, 3) /  4H*  *   ,4H  *   ,4H  *  /,  &
     nu(4, 1),nu(4, 2),nu(4, 3) /  4H*  *   ,4H  *   ,4H *   /,  &
     nu(5, 1),nu(5, 2),nu(5, 3) /  4H****   ,4H****  ,4H**** /
 
 DATA    nu(1, 4),nu(1, 5),nu(1, 6) /  4H****   ,4H* *   ,4H**** /,  &
     nu(2, 4),nu(2, 5),nu(2, 6) /  4H   *   ,4H* *   ,4H*    /,  &
     nu(3, 4),nu(3, 5),nu(3, 6) /  4H ***   ,4H****  ,4H**** /,  &
     nu(4, 4),nu(4, 5),nu(4, 6) /  4H   *   ,4H  *   ,4H   * /,  &
     nu(5, 4),nu(5, 5),nu(5, 6) /  4H****   ,4H  *   ,4H**** /
  nu(1, 7),nu(1, 8),nu(1, 9) /  4H****   ,4H****  ,4H**** /,  &
     nu(2, 7),nu(2, 8),nu(2, 9) /  4H*      ,4H   *  ,4H*  * /,  &
     nu(3, 7),nu(3, 8),nu(3, 9) /  4H****   ,4H  *   ,4H**** /,  &
     nu(4, 7),nu(4, 8),nu(4, 9) /  4H*  *   ,4H *    ,4H*  * /,  &
     nu(5, 7),nu(5, 8),nu(5, 9) /  4H****   ,4H*     ,4H**** /
 
 DATA    nu(1,10)                   /  4H****   /,  &
     nu(2,10)                   /  4H*  *   /,  &
     nu(3,10)                   /  4H****   /,  &
     nu(4,10)                   /  4H   *   /,  &
     nu(5,10)                   /  4H****   /
 
 CALL page1
 
!     GRAPH HEADING DATA
 
 IF (iframe < 0 .OR. iframe > 99999) iframe = 0
 n  = 100000
 DO  i = 1,5
   n  = n/10
   m(i) = iframe/n
   iframe = iframe - m(i)*n
   m(i) = m(i) + 1
 END DO
 n1 = m(1)
 n2 = m(2)
 n3 = m(3)
 n4 = m(4)
 n5 = m(5)
 lines = lines + 21
 itlns = itlns + 21
 WRITE  (l,20) (nu(i,n1),nu(i,n2),nu(i,n3),nu(i,n4),nu(i,n5),i=1,5)
 20 FORMAT (1H0,60X,25HF     r     a     m     e, //,  &
     5(59X,a4,2X,a4,2X,a4,2X,a4,2X,a4,/))
 WRITE  (l,30) titlec,(xtitle(i),i=1,28)
 30 FORMAT (1H0,4X,31A4,a3, /1H0,4X,15HX-axis title = ,28A4,/1H0)
 
 IF (i123 == 1) GO TO 70
 
!     DUAL FRAME TITLE FRAME
 
 WRITE  (l,60)
 WRITE  (l,40) titlel, titler
 40 FORMAT (13X,1HI,57X,3HI i,57X,1HI, /13X,2HI ,14A4,4HI i ,14A4,1HI,  &
     /13X,1HI,57X,3HI i,57X,1HI)
 WRITE  (l,50) (graph(i,6),graph(i,7),graph(i,8),i=2,3)
 50 FORMAT (12X,2(2H i,1P,e14.6,1P,e21.6,1P,e21.6,2H i))
 60 FORMAT (13X,1H+,57(1H-),3H+ +,57(1H-),1H+)
 WRITE  (l,60)
 GO TO 110
 
!     WHOLE FRAME TITLE FRAME
 
 70 WRITE  (l,80)
 WRITE  (l,90) titlel
 80 FORMAT (13X,1H+,117(1H-),1H+)
 90 FORMAT (13X,1HI,117X,1HI/13X,2HI ,14A4,60X,1HI, /13X,1HI,117X,1HI)
 WRITE  (l,100) graph(1,6),graph(1,7),graph(1,8)
 100 FORMAT (13X,1HI,1P,e14.6,37X,1P,e14.6,37X,1P,e14.6,2H i)
 WRITE  (l,80)
 
!     DUMP GRAPH
 
 110 f     = xmin - xinc
 DO  i = 1,maxplt
   temp  = f + FLOAT(i)*xinc
   i1    = (i-1)*30 + 1
   i2    = i1 + 29
   lines = lines + 1
   itlns = itlns + 1
   IF (lines-nlpp > 0) THEN
     GO TO   140
   END IF
   120 CONTINUE
   WRITE  (l,130) temp,(z(j),j=i1,i2)
   130 FORMAT (1X,1P,e11.4,1X,29A4,a3)
   CYCLE
   140 lines = 1
   WRITE  (l,150) temp,(z(j),j=i1,i2)
   150 FORMAT (1H1,1P,e11.4,1X,29A4,a3)
 END DO
 
 IF (i123 == 1) GO TO 170
 WRITE (l,60)
 GO TO 180
 170 WRITE (l,80)
 
 180 IF (exceed) WRITE (l,190) uim
 exceed = .false.
 190 FORMAT (a29,'. THERE WERE MORE POINTS BELOW THIS POINT WHICH WE',  &
     'ARE NOT PLOTTED HERE',/5X,'DUE TO CORE RESTRICTION')
 RETURN
END SUBROUTINE xygraf
