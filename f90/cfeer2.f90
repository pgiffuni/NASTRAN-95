SUBROUTINE cfeer2 (iret)
     
!     CFEER2 INITIALIZES AND CALLS CDCOMP FOR CFCNTL
 
 
 INTEGER, INTENT(OUT)                     :: iret
 LOGICAL :: qpr
 INTEGER :: filea    ,filel    ,fileu    ,scr1     ,  &
     scr2     ,scr3     ,scr4     ,scr5     ,  &
     scr6     ,scr7     ,scr8     ,scr9     ,  &
     sr1fil   ,sr2fil   ,sr3fil   ,dumm     , typout   ,bbbbar
 DOUBLE PRECISION :: det      ,mindia   ,dz(1)
 COMMON  /cdcmpx/   filea(7) ,filel(7) ,fileu(7) ,sr1fil   ,  &
     sr2fil   ,sr3fil   ,det(2)   ,power    , nz       ,mindia   ,bbbbar(5)
 COMMON  /feeraa/   dumm(36) ,scr1     ,scr2     ,scr3     ,  &
     scr4     ,scr5     ,scr6     ,scr7     ,  &
     scr8     ,scr9     ,dumq(72) ,mcblt(7) , mcbut(7)
 COMMON  /feerxc/   dumxc(21),qpr
 COMMON  /zzzzzz/   z(1)
 COMMON  /unpakx/   typout   ,irow     ,nlast    ,incr
 COMMON  /system/   ksystm(65)
 EQUIVALENCE        (z(1), dz(1)     ) ,(nout,ksystm(2))   ,  &
     (iprec,ksystm(55))
 
 itype    = iprec + 2
 iret     = 0
 filea(1) = scr1
 filel(1) = scr3
 fileu(1) = scr4
 sr1fil   = scr5
 sr2fil   = scr6
 sr3fil   = scr7
 filea(2) = dumm(3)
 filea(3) = dumm(3)
 filea(4) = dumm(4)
 filea(5) = itype
 filea(6) = 0
 filea(7) = 0
 filel(5) = itype
 nz       = korsz(z)
 bbbbar(1)= 0
 CALL cdcomp (*110,z,z,z)
 
!     ---------- SPECIAL PRINT -------------------------------
 
 IF (.NOT.qpr) GO TO 80
 WRITE  (nout,10)
 10 FORMAT (//,7H cfeer2,//)
 WRITE  (nout,20)
 20 FORMAT (1H ,13(10H----------))
 typout = itype
 irow   = 1
 nlast  = dumm(2)
 limit  = 2*nlast
 incr   = 1
 ibuf   = nz - ksystm(1) - 2
 ifilxx = scr3
 30 CALL gopen (ifilxx,z(ibuf),0)
 DO  i = 1,nlast
   WRITE  (nout,40) i
   40 FORMAT (1H ,6HCOLUMN,i4)
   CALL unpack (*50,ifilxx,z)
   IF (iprec == 2) WRITE (nout,60) (dz(j),j=1,limit)
   IF (iprec /= 2) WRITE (nout,70) ( z(j),j=1,limit)
   50 CONTINUE
 END DO
 CALL CLOSE (ifilxx,1)
 WRITE (nout,20)
 IF (ifilxx == scr4) GO TO 80
 ifilxx = scr4
 GO TO 30
 60 FORMAT (1H ,13(10H----------)/(1H ,4D25.16))
 70 FORMAT (1H ,13(10H----------)/(1H ,4E25.16))
 80 CONTINUE
 
!     --------------------------------------------------------
 
 90 DO  i = 1,7
   mcbut(i) = fileu(i)
   mcblt(i) = filel(i)
 END DO
 RETURN
 
 110 iret = 1
 GO TO 90
END SUBROUTINE cfeer2
