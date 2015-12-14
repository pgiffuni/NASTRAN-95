SUBROUTINE cfbsor (ll,ul,bx,xx,iopt)
     
 
 INTEGER, INTENT(IN)                      :: ll
 INTEGER, INTENT(IN)                      :: ul
 INTEGER, INTENT(IN OUT)                  :: bx
 INTEGER, INTENT(IN)                      :: xx
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER :: mcb(7)
 COMMON /zzzzzz/ iz(1)
 COMMON /gfbsx / jfl(7), jfu(7),jfb(7),jfx(7),jx,jprec,jsign
 COMMON /fbsx  / mfl(7),mflt(7),mfb(7),mfx(7),mx,mprec,msign
 COMMON /system/ idum(54),iprec
 
 nz = korsz(iz)
 mcb(1) = IABS(bx)
 CALL rdtrl (mcb)
 IF (iopt == 1) GO TO 100
 
!     SYMETRIC FBS
 
 mfl(1) = ll
 CALL rdtrl (mfl)
 DO  i = 1,7
   mfb(i) = mcb(i)
   mfx(i) = mcb(i)
 END DO
 mfx(1) = xx
 mx     = nz
 mfx(5) = MAX0(mfl(5),mfb(5))
 mprec  = iprec
 msign  = +1
 IF (bx < 0) msign = -1
 CALL fbs (iz,iz)
 CALL wrttrl (mfx)
 20 RETURN
 
!     UNSYMETRIC FBS
 
 100 jfl(1) = ll
 CALL rdtrl (jfl)
 jfu(1) = ul
 CALL rdtrl (jfu)
 DO  i = 1,7
   jfb(i) = mcb(i)
   jfx(i) = mcb(i)
 END DO
 jfx(1) = xx
 jx     = nz
 jprec  = iprec
 jfx(5) = MAX0(jfl(5),jfb(5))
 jsign  = +1
 IF (bx < 0) jsign = -1
 CALL gfbs (iz,iz)
 CALL wrttrl (jfx)
 GO TO 20
END SUBROUTINE cfbsor
