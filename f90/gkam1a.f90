SUBROUTINE gkam1a (mi,phidh,sdt,scr1,scr2,iopt,iout,nopp,w,nw,  &
        nosdt,lhset,i2dd,iws,scr3)
     
 
 INTEGER, INTENT(IN)                      :: mi
 INTEGER, INTENT(IN OUT)                  :: phidh
 INTEGER, INTENT(IN OUT)                  :: sdt
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN)                      :: iout
 INTEGER, INTENT(IN OUT)                  :: nopp
 REAL, INTENT(IN)                         :: w(1)
 INTEGER, INTENT(IN)                      :: nw
 INTEGER, INTENT(IN)                      :: nosdt
 INTEGER, INTENT(IN)                      :: lhset
 INTEGER, INTENT(IN)                      :: i2dd
 INTEGER, INTENT(IN)                      :: iws
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER :: sysbuf,mcb(7), FILE, NAME(2),ihh(3)
 DOUBLE PRECISION :: md,mc,ma(2),zero(2)
 DIMENSION  itab(2),itabt(13)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /system/  ksystm(65)
 COMMON /BLANK /  xx(9), kdamp
 COMMON /condas/  pi,twophi,radeg,degra,s4pisq
 COMMON /packx /  it1,it2,ii,jj,incr
 COMMON /unpakx/  it11,iii,jjj,incr1
 EQUIVALENCE      (ksystm(1),sysbuf), (ksystm(2),nout),  &
     (ksystm(55),iprec), (md,ma(1)), (mc,ma(2))
 DATA    zero  /  0.0D0,0.0D0  /
 DATA    NAME  /  4HGKAM,4H1A  /
 DATA    ihh   /  4HMHH ,4HBHH ,4HKHH  /
 DATA    g     /  0.0          /
 DATA    itabt ,  itab(1) / 4,15,21,1,25,22,2,35,23,3,45,24,4,0 /
 
 
 mc = 0.0
 IF (nopp < 0) GO TO 10
 
!     COMPUTE PHIDH(T)*I2DD*PHIDH ONTO SCR2
 
 CALL ssg2b (i2dd,phidh,0,scr1,0,2,1,iout)
 CALL ssg2b (phidh,scr1,0,scr2,1,2,1,iout)
 mcb(1) = i2dd
 CALL rdtrl (mcb)
 IF (mcb(4) /= 6) GO TO 11
 mcb(1) = scr2
 CALL rdtrl (mcb)
 mcb(4) = 6
 CALL wrttrl (mcb)
 11 CONTINUE
 mii = scr1
 10 IF (nopp < 0) mii = iout
 
!     BUILD  MII  DATA  BLOCK  = MIXF(W)
 
 lc = korsz(w(nw+1))
 nz = lc - sysbuf
 
!     RESTORE MODES
!     FILE = SCR3
 
 CALL OPEN (*130,scr3,w(nz+1),0)
 CALL fread (scr3,w,nw,1)
 CALL CLOSE (scr3,1)
 FILE = mi
 CALL OPEN (*170,mi,w(nz+1),0)
 imi  = 0
 CALL skprec (mi,iws)
 21 CONTINUE
 nz   = nz - sysbuf
 ibuf = nz - sysbuf
 icrq = -ibuf
 IF (icrq > 0) GO TO 150
 CALL gopen (mii,w(nz+1),1)
 CALL makmcb (mcb,mii,lhset,6,iprec)
 IF (kdamp == 1) mcb(5) = mcb(5) + 2
 
!     SET UP FOR  PACK  AND  UNPACK
 
 it1  = 2
 IF (kdamp == 1) it1 = 4
 it2  = mcb(5)
 incr = 1
 it11 = 2
 incr1= 1
 DO  i = 1,nw
   mc   = 0.0
   k    = iws + i - 1
   ii   = i
   jj   = i
   iii  = k
   jjj  = k
   IF (imi /= 0) GO TO 85
   CALL unpack (*160,mi,md)
   22 CONTINUE
   SELECT CASE ( iopt )
     CASE (    1)
       GO TO 30
     CASE (    2)
       GO TO 50
     CASE (    3)
       GO TO 40
   END SELECT
   
!     BUILDING  MHH
   
   30 CALL pack (md,mii,mcb)
   CYCLE
   
!     BUILDING  KHH
   
   40 md = md*w(i)*w(i)
   IF (kdamp /= 1) GO TO 30
   ASSIGN 45 TO iret
   IF (nosdt > 0) GO TO 70
   45 mc = g*md
   GO TO 30
   
!     BUILDING  BHH
   
   50 CONTINUE
   IF (kdamp == 1) GO TO 61
   ASSIGN 60 TO iret
   IF (nosdt > 0) GO TO 70
   60 md = md*w(i)*g
   GO TO 30
   61 md = 0.0
   GO TO 30
   
!     LOOK UP G(W)  IN  SDT
   
   70 IF (itab(1) > 0) GO TO 80
   itab(1) = 1
   itab(2) = nosdt
   CALL pretab (sdt,w(nw+1),w(nw+1),w(ibuf),ibuf-1,iz,itab(1),itabt)
   80 CALL tab (itab(2),w(i)/twophi,g)
   GO TO iret, (60,45)
   
!     PICK UP MODAL MASS FROM LAMA
   
   85 CALL fread (mi+1,0,-5,0)
   CALL fread (mi+1,xmass,1,0)
   CALL fread (mi+1,0,-1,0)
   md = xmass
   GO TO 22
   
!     ADD  INTERPOLATION HERE
   
 END DO
 CALL CLOSE (mi  ,1)
 CALL CLOSE (mi+1,1)
 NE = lhset - nw
 IF (NE <= 0) GO TO 110
 DO  i = 1,NE
   CALL pack (zero,mii,mcb)
 END DO
 110 CALL  wrttrl (mcb)
 CALL  CLOSE (mii,1)
 RETURN
 
!     ERROR MESAGES
 
 130 ip1  = -1
 140 CALL mesage (ip1,FILE,NAME)
 RETURN
 150 ip1  = -8
 FILE = icrq
 GO TO 140
 160 WRITE  (nout,9001) sfm,ihh(iopt)
 9001 FORMAT (a25,' 2203, NULL COLUMN FOUND IN MI FILE DURING ASSEMBLY',  &
     ' OF ',a4,' MATRIX BY GKAM MODULE.')
 ip1  = -37
 GO TO 140
 
!     USE LAMA RATHER THAN MI
 
 170 CONTINUE
 CALL gopen (mi+1,w(nz+1),0)
 CALL skprec (mi+1,1)
 CALL fread (mi+1,mcb,-7*(iws-1),0)
 imi = 1
 GO TO 21
END SUBROUTINE gkam1a
