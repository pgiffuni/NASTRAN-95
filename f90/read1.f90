SUBROUTINE read1 (dm,mr,scr1,scr2,scr3,phia,uset,nr1,lama,scr4)
     
 
 INTEGER, INTENT(IN OUT)                  :: dm
 INTEGER, INTENT(IN)                      :: mr
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN OUT)                  :: scr2
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER, INTENT(IN OUT)                  :: phia
 REAL, INTENT(IN OUT)                     :: uset
 INTEGER, INTENT(OUT)                     :: nr1
 INTEGER, INTENT(IN OUT)                  :: lama
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER :: imr(7),sysbuf, iscr1(7),  nam(2)
 DOUBLE PRECISION :: dcore(1),si,term
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /system/  sysbuf,nout,ksystm(63)
 COMMON /zzzzzz/  core(1)
 COMMON /unpakx/  itb,ii,jj,incur
 COMMON /packx /  ita1,itb1,ii1,jj1,incur1
 COMMON /bitpos/  um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 EQUIVALENCE      (dcore(1),core(1))
 DATA    nam   /  4HREAD,4H1    /
 
!     BRING MR INTO CORE
 
 lc = korsz(core) - sysbuf
 CALL gopen (mr,core(lc+1),0)
 imr(1) = mr
 CALL rdtrl (imr)
 nr   = imr(2)
 nr1  = nr
 ii   = 1
 jj   = nr
 incur= 1
 itb  = imr(5)
 nr2  = itb*nr
 ivi  = nr*nr
 iphi = ivi
 ivi2 = itb*ivi
 ialph= 2*ivi
 iloop= 0
 k    = 0
 DO  i = 1,nr
   CALL unpack (*12,mr,core(k+1))
   GO TO 16
   
!     NULL COLUMN
   
   12 DO  j = 1,nr2
     core(j+k) = 0.0
   END DO
   16 kkk = k + ivi2
   DO  j = 1,nr2
     core(j+kkk) = 0.0
   END DO
   IF (itb == 1) GO TO 18
   kkk = kkk/2
   dcore(kkk+i) = 1.0D0
   GO TO 19
   18 core(kkk+i) = 1.0
   19 k = k + nr2
 END DO
 CALL CLOSE (mr,1)
 
!     COMPUTE SI
 
 IF (itb /= 2) GO TO 35
 30 si = 0.0D0
 DO  i = 1,nr
   term = 0.0D0
   DO  j = 1,nr
     k  = (j-1)*nr + i
     kk = ivi + j
     term = term + dcore(k)*dcore(kk)
   END DO
   k  = ivi + i
   si = si + term*dcore(k)
 END DO
 IF (si > 0.0D0) GO TO 51
 53 WRITE  (nout,52) ufm
 52 FORMAT (a23,' 2200, INCONSISTENT RIGID BODY SYSTEM.')
 CALL mesage (-61,0,nam)
 51 CONTINUE
 si = 1.0D0/DSQRT(si)
 
!     CONVERT VI INTO PHI
 
 DO  i = 1,nr
   k = ivi + i
   dcore(k) = dcore(k) *si
 END DO
 iloop = iloop + 1
 IF (iloop == nr) GO TO 120
 
!     CALCULATE ALPHAJ
 
 DO  j = 1,iloop
   k = ialph + j
   dcore(k) = 0.0D0
   DO  i = 1,nr
     term = 0.0D0
     DO  l = 1,nr
       kk  = (l-1)*nr + i
       kkk = ivi + nr + l
       term = term + dcore(kk)*dcore(kkk)
     END DO
     kk = iphi + (j-1)*nr + i
     dcore(k) = dcore(k)+term*dcore(kk)
   END DO
 END DO
 
!     COMPUTE NEXT V VECTOR
 
 DO  i = 1,nr
   term = 0.0D0
   DO  j = 1,iloop
     kk = ialph + j
     k  = iphi + (j-1)*nr + i
     term = term + dcore(kk)*dcore(k)
   END DO
   k = ivi + nr + i
   dcore(k) = dcore(k) - term
 END DO
 ivi = ivi + nr
 GO TO 30
 35 ssi = 0.0
 DO  i = 1,nr
   sterm = 0.0
   DO  j = 1,nr
     k  = (j-1)*nr + i
     kk = ivi + j
     sterm = sterm + core(k)*core(kk)
   END DO
   k   = ivi + i
   ssi = ssi + sterm*core(k)
 END DO
 IF (ssi <= 0.0) GO TO 53
 ssi = 1.0/SQRT(ssi)
 
!     CONVERT VI INTO PHI
 
 DO  i = 1,nr
   k = ivi + i
   core(k) = core(k)*ssi
 END DO
 iloop = iloop + 1
 IF (iloop == nr) GO TO 120
 
!     CALCULATE ALPHAJ
 
 DO  j = 1,iloop
   k = ialph + j
   core(k) = 0.0
   DO  i = 1,nr
     sterm = 0.0
     DO  l = 1,nr
       kk  = (l-1)*nr + i
       kkk = ivi + nr + l
       sterm = sterm + core(kk)*core(kkk)
     END DO
     kk = iphi + (j-1)*nr + i
     core(k) = core(k) + sterm*core(kk)
   END DO
 END DO
 
!     COMPUTE NEXT V VECTOR
 
 DO  i = 1,nr
   sterm = 0.0
   DO  j = 1,iloop
     kk = ialph + j
     k  = iphi + (j-1)*nr + i
     sterm = sterm + core(kk)*core(k)
   END DO
   k = ivi + nr + i
   core(k) = core(k) - sterm
 END DO
 ivi = ivi + nr
 GO TO 35
 
!     PACK PHIRO
 
 120 ita1 = itb
 itb1 = itb
 ii1  = 1
 jj1  = nr
 incur1 = 1
 CALL gopen (scr1,core(lc+1),1)
 CALL makmcb (iscr1,scr1,nr,1,itb)
 DO  i = 1,nr
   k = ivi2 + (i-1)*nr2
   CALL pack (core(k+1),scr1,iscr1)
 END DO
 CALL CLOSE (scr1,1)
 CALL wrttrl (iscr1(1))
 
!     COMPUTE PHILO = DM*PHIRO
 
 CALL ssg2b (dm,scr1,0,scr2,0,itb,1,scr4)
 
!     MERGE PHIRP AND PHILO TO FORM PHIA
 
 CALL sdr1b (scr3,scr2,scr1,scr4,ua,ul,ur,uset,0,0)
 CALL gopen (scr4,core(lc+1),0)
 lc = lc - sysbuf
 CALL gopen (phia,core(lc+1),1)
 imr(1) = scr4
 CALL rdtrl (imr(1))
 nprob = imr(3)
 dcore(1) = 0.d0
 jj = nprob
 incur = 1
 i3 = 3
 DO  j = 1,nr
   ii = 0
   CALL unpack (*150,scr4,core(i3))
   ii1 = ii
   jj1 = jj
   CALL pack (core(i3),phia,iscr1)
   CYCLE
   
!     NULL COLUMN
   
   150 ii1 = 1
   jj1 = 1
   CALL pack (core,phia,iscr1)
 END DO
 CALL CLOSE (scr4,1)
 CALL CLOSE (phia,1)
 lc = lc + sysbuf
 
!     PUT NR ZEROS ON LAMA
 
 CALL gopen (lama,core(lc+1),1)
 dcore(1) = 0.d0
 DO  i = 1,nr
   CALL WRITE (lama,core,itb,1)
 END DO
 CALL CLOSE (lama,2)
 RETURN
END SUBROUTINE read1
