SUBROUTINE cyct2a (kaa,kxx,g1,g2,scr1,scr2,scr3)
     
 INTEGER, INTENT(IN)                      :: kaa
 INTEGER, INTENT(IN)                      :: kxx
 INTEGER, INTENT(IN)                      :: g1
 INTEGER, INTENT(IN)                      :: g2
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER :: mcb(7)
 COMMON /system/idum(54),iprec
 
 mcb(1)=kaa
 CALL rdtrl(mcb)
 IF(mcb(1) <= 0)GO TO 30
 mcb(1)=kxx
 CALL rdtrl(mcb)
 IF(mcb(1) <= 0)GO TO 30
 isc2=scr2
 mcb(1)=g1
 CALL rdtrl (mcb)
 IF (mcb(2) <= 0)isc2=0
 mcb(1)=g2
 CALL rdtrl(mcb)
 IF(mcb(2) <= 0)GO TO 10
 isc=scr2
 iout=1
 isc1=kxx
 GO TO 20
 10 isc=kxx
 isc1=scr2
 iout=0
 
!     NO FIRST TERM IF ISC2=0, NO SECOND TERM IF IOUT=0
 
 20 IF (isc2 == 0) GO TO 25
 
!     COMPUTE FIRST TERM
 
 CALL ssg2b(kaa,g1,0,scr1,0,iprec,1,scr2)
 CALL ssg2b(g1,scr1,0,isc,1,iprec,1,isc1)
 
!     COMPUTE SECOND TERM
 
!     COMPUTE SECOND TERM
 
 25 IF(iout == 0) GO TO 29
 CALL ssg2b(kaa,g2,0,scr1,0,iprec,1,kxx)
 CALL ssg2b(g2,scr1,isc2,kxx,1,iprec,1,scr3)
 29 mcb(1) = kxx
 30 RETURN
END SUBROUTINE cyct2a
