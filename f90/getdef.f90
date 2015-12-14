SUBROUTINE getdef (dfrm,ph,mag,conv,plttyp,buf,gpt,d)
     
 
 INTEGER, INTENT(IN)                      :: dfrm
 REAL, INTENT(IN OUT)                     :: ph
 INTEGER, INTENT(IN OUT)                  :: mag
 REAL, INTENT(IN)                         :: conv
 INTEGER, INTENT(IN)                      :: plttyp
 INTEGER, INTENT(IN OUT)                  :: buf(1)
 INTEGER, INTENT(IN OUT)                  :: gpt(1)
 REAL, INTENT(OUT)                        :: d(3,1)
 INTEGER :: siln,rew,sp,gp,gpx,sil1,sil2, trl(7),TYPE
 REAL :: maxdef
 COMMON /BLANK / ngp,lsil,skp11(3),ngpset,skp12(4),skp2(6),msil
 COMMON /xxparm/ pbufsz,ploter(5),penpap(30),scale(4),maxdef
 COMMON /zntpkx/ defc(4),siln,last
 EQUIVALENCE     (defval,defc(1))
 DATA    inprew, rew / 0,1 /
 
 last = 0
 k    = 3*ngpset
 DO  i = 1,k
   d(i,1) = 0.0
 END DO
 trl(1) = dfrm
 CALL rdtrl (trl(1))
 IF (trl(5) <= 0) RETURN
 sp = trl(5)
 ASSIGN 140 TO TYPE
 
!     NOTE TRANSIENT RESPONSE HAS SP = 1
 
 IF (sp < 3) GO TO 30
 ASSIGN 130 TO TYPE
 IF (mag /= 0) GO TO 30
 ASSIGN 120 TO TYPE
 sn = SIN(ph)*conv
 cn = COS(ph)*conv
 IF (plttyp == 2) GO TO 20
 
!     DISPLACEMENT OR ACCELERATION
 
 i1 = 1
 i2 = sp - 1
 IF (plttyp == 3 .OR. plttyp == 4) cn = -cn
 GO TO 30
 
!     VELOCITY
 
 20 i1 = sp - 1
 i2 = 1
 30 CONTINUE
 maxdef = 0.
 CALL intpk (*170,dfrm,0,sp,0)
 gp   = 0
 siln = 0
 CALL gopen (msil,buf(1),inprew)
 CALL fread (msil,sil2,1,0)
 
!     -GP- = PREVIOUS EXISTENT GRID POINT IN THIS SET. FIND NEXT ONE.
 
 40 k = gp + 1
 DO  gpx = k,ngp
   IF (gpt(gpx) /= 0) GO TO 60
 END DO
 sil1 = lsil + 1
 GO TO 100
 60 IF (gpx /= gp+1) GO TO 70
 sil1 = sil2
 GO TO 80
 70 gp = gp + 1
 CALL fread (msil,sil2,1,0)
 GO TO 60
 
!     -SIL1- = SIL NUMBER OF NEXT EXISTENT GRID POINT. READ SIL NUMBER
!              OF NEXT GRID POINT.
 
 80 gp = gpx
 gpx = IABS(gpt(gp))
 IF (gp == ngp) sil2 = lsil + 1
 IF (gp /= ngp) CALL fread (msil,sil2,1,0)
 
!     READ NEXT DEFORMATION VALUE AT THIS EXISTING GRID POINT.
 
 90 IF (siln <= lsil .AND. siln >= sil1) GO TO 150
 100 IF (last /= 0) GO TO 160
 110 CALL zntpki
 GO TO TYPE, (120,130,140)
 120 defval = defc(i1)*cn - defc(i2)*sn
 GO TO 140
 130 defval = conv*SQRT(defc(1)**2 + defc(sp-1)**2)
 140 IF (ABS(defval) > maxdef) maxdef = ABS(defval)
 GO TO 90
 150 IF (siln > sil1+2 .OR. siln >= sil2) GO TO 40
 k = siln - sil1 + 1
 d(k,gpx) = defval
 IF (last == 0) THEN
   GO TO   110
 END IF
 
 160 CALL CLOSE (msil,rew)
 170 RETURN
END SUBROUTINE getdef
