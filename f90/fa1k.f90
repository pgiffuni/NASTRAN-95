SUBROUTINE fa1k (imeth,k,rho,outfil,ico)
     
!     FA1K BUILDS AN INTERPOLATED MATRIX ON OUTFIL FROM QHHL OR FSAVE
 
 
 INTEGER, INTENT(IN OUT)                  :: imeth
 REAL, INTENT(OUT)                        :: k
 REAL, INTENT(OUT)                        :: rho
 INTEGER, INTENT(IN)                      :: outfil
 INTEGER, INTENT(IN)                      :: ico
 LOGICAL :: NEW
 INTEGER :: sysbuf,out,buff,buff1,floop,ns(2),TYPE,trl(7),  &
      fsave,qhhl,scr2,scr3,scr4,mcb(7)
 
 DIMENSION       z(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /system/ sysbuf,out,dum(52),iprec
 COMMON /BLANK / floop
 EQUIVALENCE     (iz(1),z(1))
 DATA    fsave / 201/, qhhl /104/, scr2,scr3,scr4 /302,303,304/
 DATA    ns    / 4HFA1K,4H    /
 
 ncore = korsz(iz) - ico
 buff  = ncore - sysbuf
 buff1 = buff  - sysbuf
 trl(1)= fsave
 CALL rdtrl (trl)
 
!     READ IN DEPENDENT POINTS AND SET K AND RHO
 
 jj    = trl(3)*3
 ifil  = fsave
 CALL gopen (fsave,iz(buff+1),0)
 CALL READ (*430,*5,fsave,z,jj,1,nwr)
 5 CONTINUE
 i     = (floop-1)*3 + 1
 cmach = z(i  )
 k     = z(i+1)
 rho   = z(i+2)
 incr1 = 1
 incr  = 1
 ii    = 1
 inn   = 1
 SELECT CASE ( imeth )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 100
 END SELECT
 
!     SURFACE SPLINE INTERPOLATION
 
 10 IF (floop /= 1) GO TO 70
 
!     SET UP CALL TO SPLINE INTERPOLATOR
 
 nrho = trl(7)
 ji   = nrho*3
 j    = 1
 DO  i = 1,jj,ji
   z(j  ) = z(i  )
   z(j+1) = z(i+1)
   j    = j + 2
 END DO
 nd   = jj/ji
 ni   = trl(4)
 TYPE = 1
 idp  = 1
 iip  = nd*2 + idp
 ig   = iip + 2*ni
 ni2  = ni*2
 CALL READ (*430,*60,fsave,z(iip),ni2,1,nwr)
 60 CALL fwdrec (*430,fsave)
 CALL CLOSE (fsave,2)
 
!     REWRITE QHHL SO EACH LIST MATRIX IS A COLUMN
 
 GO TO 300
 15 j  = 1
 20 ji = ig
 GO TO 350
 30 IF (j == ncol) GO TO 40
 j  = j + 1
 GO TO 20
 40 CALL CLOSE (qhhl,1)
 CALL CLOSE (outfil,1)
 CALL wrttrl (trl)
 GO TO 200
 
!     GET A COLUMN FROM FSAVE AND BUILD QHH ON OUTFIL
 
 70 nf = 2 + (floop-1)/trl(7)
 DO  i = 1,nf
   CALL fwdrec (*430,fsave)
 END DO
 85 iout = trl(6)
 iti  = iout
 ito  = iout
 nwc  = 1
 IF (ito == 2 .OR. ito == 3) nwc = 2
 IF (ito == 4) nwc = 4
 mcb(1) = qhhl
 CALL rdtrl (mcb)
 nc   = mcb(3)
 nn   = nc
 nnn  = nc*nc
 CALL unpack (*410,fsave,z)
 ij   = 1
 CALL CLOSE (fsave,1)
 CALL gopen (outfil,iz(buff+1),1)
 mcb(1) = outfil
 mcb(2) = 0
 mcb(3) = nc
 mcb(4) = 1
 mcb(5) = iout
 mcb(6) = 0
 mcb(7) = 0
 DO  i = 1,nc
   CALL pack (z(ij),outfil,mcb)
   ij = ij + nc*nwc
 END DO
 CALL CLOSE (outfil,1)
 CALL wrttrl (mcb)
 GO TO 450
 
!     LINEAR SPLINE INTERPOLATION
 
 
!     IS A GOOD MATRIZ ON FSAVE
 
 100 eps = .001
 NEW = .true.
 ni  = trl(4)
 IF (floop == 1) GO TO 110
 ok    = z(i-2)
 omach = z(i-3)
 IF (ABS(cmach-omach) < eps) NEW = .false.
 
!     REWRITE QHHL IF NEW IS TRUE
 
 IF (.NOT.NEW) GO TO 180
 IF (floop /= 1) GO TO 120
 
!     TEST TO SEE IF QHHL HAS ENOUGH MACH NUMBERS
 
 110 nip  = ni*2
 nogo = 0
 iip  = jj + 1
 CALL READ (*430,*111,fsave,z(iip),nip,1,nwr)
 111 CALL bckrec (fsave)
 temp = 0.0
 DO  i = 1,jj,3
   IF (temp == z(i)) CYCLE
   temp = z(i)
   nf   = 0
   DO  j = 1,nip,2
     IF (temp-z(iip+j-1) < eps) nf = nf + 1
   END DO
   IF (nf > 1) CYCLE
   WRITE (out,400) ufm,temp
   nogo = 1
 END DO
 IF (nogo == 1 ) GO TO 410
 120 j   = 1
 nrd = 0
 DO  i = 1,jj,3
   IF (ABS(cmach-z(i)) < eps) GO TO 126
   CYCLE
   126 IF (z(i+2) /= rho) CYCLE
   z(j  ) = z(i  )
   z(j+1) = z(i+1)
   j   = j + 2
   nrd = nrd + 1
 END DO
 idp = 1
 iip = nrd*2 + idp
 ni2 = ni*2
 CALL READ (*430,*130,fsave,z(iip),ni2,1,nwr)
 130 CALL fwdrec (*430,fsave)
 CALL CLOSE (fsave,2)
 GO TO 300
 135 ig  = iip + ni*2
 nf  = 0
 ik  = 1
 ifil= qhhl
 jj  = 2*ni + 1
 i   = 1
 138 IF (ABS(cmach-z(iip+i-1)) < eps) GO TO 140
 
!     SKIP MATRIX
 
 DO  j = 1,ncm
   CALL fwdrec (*430,qhhl)
 END DO
 GO TO 150
 140 z(iip+ik) = z(iip+i)
 ik = ik + 2
 nf = nf + 1
 ji = ig
 GO TO 350
 150 CONTINUE
 i  = i + 2
 IF (i == jj) GO TO 160
 GO TO 138
 160 CALL CLOSE (qhhl,1)
 CALL CLOSE (outfil,1)
 CALL wrttrl (trl)
 
!     SET UP CALL TO SPLINE INTERPOLATION
 
 TYPE = -1
 nd   = nrd
 ni   = nf
 GO TO 200
 
!     GET COLUMN FROM FSAVE AND BUILD QHH
 
 170 CALL gopen (fsave,iz(buff+1),0)
 ij = 3
 171 DO  i = 1,ij
   CALL fwdrec (*430,fsave)
 END DO
 GO TO 85
 180 IF (ok-k == 0.0) GO TO 190
 trl(7) = trl(7) + 1
 190 ij = trl(7) + 1
 CALL wrttrl (trl)
 GO TO 171
 
!     CALL MINTRP
 
 200 ig   = iip + 2*ni
 nc   = ncore - ig
 nogo = 0
 CALL mintrp (ni,z(iip),nd,z(idp),TYPE,0,0,0.0,outfil,scr2,scr3,  &
     scr4,z(ig),nc,nogo,iprec)
 IF (nogo == 1) GO TO 410
 
!     INTERPOLATED MATRIX IS ON SCR2 MOVE TO FSAVE
 
 CALL OPEN (*430,fsave,iz(buff+1),3)
 CALL gopen (scr2,iz(buff1+1),0)
 trl(1) = scr2
 CALL rdtrl (trl)
 ncol = trl(2)
 nn   = trl(3)
 nnn  = nn
 iti  = trl(5)
 ito  = iti
 iout = iti
 trl(1) = fsave
 trl(2) = 0
 trl(6) = 0
 trl(7) = 0
 i = 1
 210 CALL unpack (*410,scr2,z)
 CALL pack (z,fsave,trl)
 IF (i == ncol) GO TO 230
 i = i + 1
 GO TO 210
 230 CALL CLOSE (scr2,1)
 CALL CLOSE (fsave,1)
 CALL rdtrl (trl)
 trl(6) = ito
 IF (imeth == 2) trl(7) = 1
 CALL wrttrl (trl)
 GO TO 170
 
!     SET UP COLUMN - MATRIX COPY
 
 300 CALL gopen (qhhl,iz(buff+1),0)
 trl(1) = qhhl
 CALL rdtrl (trl)
 ncol = trl(2)/trl(3)
 ncm  = trl(3)
 CALL gopen (outfil,iz(buff1+1),1)
 nnn  = ncm
 nn   = ncm*ncm
 iti  = trl(5)
 ito  = iti
 iout = iti
 nwc  = 1
 IF (ito == 2 .OR. ito == 3) nwc = 2
 IF (ito == 4) nwc = 4
 trl(1) = outfil
 trl(2) = 0
 trl(3) = nn
 trl(6) = 0
 trl(7) = 0
 SELECT CASE ( imeth )
   CASE (    1)
     GO TO 15
   CASE (    2)
     GO TO 135
 END SELECT
 
!     MAKE A COLUMN INTO MATRIX
 
 350 DO  ilop = 1,ncm
   CALL unpack (*360,qhhl,z(ji))
   GO TO 380
   360 n = ncm*nwc
   DO  ij = 1,n
     z(ji+ij-1) = 0.0
   END DO
   380 ji = ji + ncm*nwc
 END DO
 CALL pack (z(ig),outfil,trl)
 SELECT CASE ( imeth )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 150
 END SELECT
 
!     ERROR MESSAGES
 
 400 FORMAT (a23,' 2270, LINEAR INTERPOLATION WITHOUT ENOUGH IND. ',  &
     'MACH NUMBERS EQUAL TO DEP. MACH ',f10.4)
 410 WRITE  (out,420) ufm
 420 FORMAT (a23,' 2271, INTERPOLATION MATRIX IS SINGULAR')
 GO TO  440
 430 CALL mesage (-3,ifil,ns)
 440 CALL mesage (-61,0,ns)
 450 RETURN
END SUBROUTINE fa1k
