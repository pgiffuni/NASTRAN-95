SUBROUTINE opt2b (ipr,pr,pl,rr)
     
 
 INTEGER, INTENT(IN)                      :: ipr(1)
 REAL, INTENT(IN OUT)                     :: pr(1)
 REAL, INTENT(IN)                         :: pl(1)
 REAL, INTENT(IN)                         :: rr(1)
 INTEGER :: count, outtap,sysbuf,iy(1)
 REAL :: y(1),z(8)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / skp1(2),count,skp2(6),nwdsp, skp3(6),nprw,nklw,ntotl,conv
 COMMON /zzzzzz/ core(1)
 COMMON /system/ sysbuf,outtap
 EQUIVALENCE     (core(1),z(1),MAX), (eps,z(2)), (gama,z(3)),  &
     (iprnt,z(7)), (iy(1),y(1),z(8))
!     EQUIVALENT ARE  (IPR,PR)
 
 nmes = 0
 ch   = 1.0
 
 DO  np = 1,nprw,nwdsp
   alph= pr(np+4)
   i   = 1
   icp = ntotl - 4
   3 icp = icp+4
   IF (iy(icp) <=  0) GO  TO 5
   IF (iy(icp) /= np) GO  TO 3
   
!     SPECIAL HANDLING OF TRIM6
   
   4 alph = y(icp+i)
   
   5 IF (alph < 0.0) THEN
     GO TO    70
   ELSE IF (alph == 0.0) THEN
     GO TO    40
   END IF
   
!     POSITIVE ALPHA, CALCULATE PNEW
   
   10 irr = (np+nwdsp)/nwdsp
   IF (ABS(gama-1.0) < 1.0E-4) ch = 0.25*rr(irr) + 0.75
   pnew = pr(np+3)*((alph/(alph+(1.0-alph)*gama))**ch)
   IF (ipr(np+5) == 0) GO TO 30
   
!     COMPARE TO LIMIT DATA
   
   kpl  = ipr(np+5)
   delp = pnew/pr(np+2)
   IF (delp < pl(kpl)) GO TO 20
   kpl = kpl + 1
   IF (delp <= pl(kpl) .OR. pl(kpl) == 0) GO TO 30
   
!     RECALCULATE ALPHA, PNEW  BASED ON THE LIMIT
   
   20 pnew = pr(np+2)*pl(kpl)
   alph =-pnew*gama/(pnew*(1.0-gama)-pr(np+3))
   
   30 pr(np+4) = alph
   IF (np == iy(icp)) y(icp+i) = alph
   GO TO 80
   
!     ZERO STRESS INPUT, CHANGE ALPH TO 0.0001
   
   40 IF (iprnt == 0 .OR. nmes >= 100) GO TO 60
   nmes = nmes + 1
   CALL page2 (-2)
   WRITE  (outtap,50) uwm,ipr(np)
   50 FORMAT (a25,' 2303, FULLY-STRESSED DESIGN DETECTED ZERO STRESS ',  &
       'FOR PROPERTY',i9, /5X,'CHECK PROPERTY CARD OR UNLOADED ', 'ELEMENT(S)')
   60 alph = 1.0E-4
   GO TO 10
   
!     NO CHANGE IN ALPH (-1.0 DETECTED)
   
   70 alph = -1.0E0
   IF (np == iy(icp)) GO TO 30
   
   80 IF (np /= iy(icp)) CYCLE
   i = i + 1
   IF (i <= 3) GO TO 4
   icp = icp + 4
   
 END DO
 
 RETURN
END SUBROUTINE opt2b
