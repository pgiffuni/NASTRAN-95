SUBROUTINE gkam
     
!     ROUTINE TO ASSEMBLE MODAL MATRICES
 
!     INPUTS = 9
 
!     USETD,PHIA,MI,LAMA,SDT,M2DD,B2DD,K2DD,CASECC
 
!     OUTPUTS = 4
 
!     MHH,BHH,KHH,PHIDH
 
!     SCRATCHES = 4
 
 INTEGER :: usetd,b2dd,sdt,phia,phidh,bhh,scr1,scr2,scr3,  &
     phidh1,sysbuf,casecc,NAME(2)
 REAL :: lfreq
 DIMENSION       mcb(7),icore(2),BLOCK(11),iblock(11)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / noue,nlmode,lfreq,hfreq,nom2dd,nob2dd,nok2dd,  &
     noncup,nmode,kdamp
 COMMON /packx / it1,it2,ii,jj,incr
 COMMON /unpakx/ it3,ii1,jj1,incr1
 COMMON /condas/ pi,twophi,radeg,degra,s4pisq
 COMMON /zzzzzz/ core(1)
 COMMON /system/ sysbuf,nout
 
 EQUIVALENCE     (core(1),icore(1)),(iblock(1),BLOCK(1))
 
 DATA    NAME  / 4HGKAM,4H    /
 DATA    iblock(1),iblock(7),BLOCK(2),BLOCK(8)   / 1,1,1.0,1.0 /
 DATA    usetd , phia,mi, lama,sdt,m2dd,b2dd,k2dd/  &
     101   , 102, 103,104, 105,106, 107, 108 /
 DATA    mhh   , bhh,khh,phidh/ 201   , 202,203,204  /
 DATA    scr1  , scr2,scr3,phidh1,casecc / 301   , 302 ,303 ,304   ,109    /
 
 
!     PICK UP AND STORE SELECTED MODES, SAVING EIGENVECTORS
 
 lc1  = korsz(core)
 nz   = lc1 - sysbuf
 icrq = 2*sysbuf - nz
 IF (icrq > 0) GO TO 220
 
!     FIND SELECTED SDT INTO CASECC
 
 CALL gopen (casecc,core(nz+1),0)
 CALL fread (casecc,icore,166,1)
 CALL CLOSE (casecc,1)
 i149  = 149
 nosdt = icore(i149)
 
!     OPEN  LAMA, PHIA, AND PHI0H
 
 CALL gopen (lama,core(nz+1),0)
 CALL skprec (lama,1)
 nz = nz - sysbuf
 CALL gopen (phia,core(nz+1),0)
 icore(1) = phia
 CALL rdtrl (icore)
 nvect = icore(2)
 nz = nz - sysbuf
 IF (noue < 0) phidh1 = phidh
 CALL gopen (phidh1,core(nz+1),1)
 mcb(1) = phia
 CALL rdtrl (mcb)
 mcb(1)= phidh1
 it1   = mcb(5)
 it2   = it1
 it3   = it1
 incr  = 1
 incr1 = 1
 ii    = 1
 ii1   = 1
 jj    = mcb(3)
 jj1   = jj
 mcb(2)= 0
 mcb(6)= 0
 mcb(7)= 0
 isw   = 1
 modes = 1
 DO  i = 1,nvect
   CALL READ (*190,*40,lama,core(nz-6),7,0,iflag)
   
!     PICK UP FREQUENCY
   
   f = core(nz-2)
   IF (nlmode == 0) GO TO 50
   
!     ACCEPT LAMA
   
   20 core(modes) = f*twophi
   modes = modes + 1
   CALL unpack (*210,phia,core(modes))
   GO TO 30
   
!     FREQUENCY RANGE SPECIFICATION
   
   50 IF (f > hfreq) GO TO 40
   IF (f >= lfreq) GO TO 20
   CALL skprec (phia,1)
   isw = isw + 1
   CYCLE
   30 CALL pack (core(modes),phidh1,mcb)
   IF (nlmode == 0) CYCLE
   IF (modes > nlmode) GO TO 40
 END DO
 
!     DONE
 
 40 CALL CLOSE (lama,1)
 CALL CLOSE (phia,1)
 CALL CLOSE (phidh1,1)
 CALL wrttrl (mcb)
 GO TO 60
 
!     BUILD PHIDH
 
 60 lhset = modes - 1
 nmode = isw
 IF (lhset <= 0) GO TO 230
 IF (noue  < 0) GO TO 70
 CALL gkam1b (usetd,scr1,scr2,phidh,phidh1,modes,core,lhset,noue, scr3)
 
!     FORM H MATRICES
 
 70 modes = modes - 1
 
!     SAVE MODES ON SCRATCH3 IN CASE DMI WIPES THEM OUT
 
 nz = lc1 - sysbuf
 CALL OPEN  (*250,scr3,core(nz+1),1)
 CALL WRITE (scr3,core(1),modes,1)
 CALL CLOSE (scr3,1)
 noncup = 1
 
!     FORM  MHH
 
 CALL  gkam1a (mi,phidh,sdt,scr1,scr2,1,mhh,nom2dd,core(1),modes,  &
     nosdt,lhset,m2dd,isw,scr3)
 IF (nom2dd < 0) GO TO 80
 CALL ssg2c (scr1,scr2,mhh,1,iblock(1))
 80 CONTINUE
 
!     FORM  BHH
 
 IF (nosdt == 0 .AND. nob2dd < 0) GO TO 90
 CALL gkam1a (mi,phidh,sdt,scr1,scr2,2,bhh,nob2dd,core(1),modes,  &
     nosdt,lhset,b2dd,isw,scr3)
 IF (nob2dd < 0) GO TO 90
 CALL ssg2c (scr1,scr2,bhh,1,iblock(1))
 90 CONTINUE
 
!     FORM  KHH
 
 CALL gkam1a (mi,phidh,sdt,scr1,scr2,3,khh,nok2dd,core(1),modes,  &
     nosdt,lhset,k2dd,isw,scr3)
 IF (nok2dd < 0) GO TO 100
 CALL ssg2c (scr1,scr2,khh,1,iblock(1))
 100 CONTINUE
 IF (nob2dd < 0 .AND. nom2dd < 0 .AND. nok2dd < 0) noncup = -1
 RETURN
 
!     ERROR MESAGES
 
 120 ip1 = -1
 130 CALL mesage (ip1,ip2,NAME)
 190 ip2 = lama
 ip1 = -3
 GO TO 130
 210 WRITE  (nout,215) sfm
 215 FORMAT (a25,' 2204, UNPACK FOUND NULL COLUMN IN PHIA FILE IN ',  &
     'GKAM MODULE.')
 ip1 = -37
 GO TO 130
 220 ip1 = -8
 FILE = icrq
 GO TO 130
 
!     NO MODES SELECTED
 
 230 ip1 = -47
 GO TO 130
 250 ip2 = scr3
 GO TO 120
END SUBROUTINE gkam
