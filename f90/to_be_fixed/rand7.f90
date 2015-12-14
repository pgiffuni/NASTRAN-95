SUBROUTINE rand7 (ifile,nfile,psdl,dit,icoup,nfreq,npsdl,ntau,  &
        ltab,casecc,xycdb)
     
!     STORES STUFF IN CORE FOR LATER RANDOM ANALSIS
 
 
 INTEGER, INTENT(IN)                      :: ifile(1)
 INTEGER, INTENT(IN)                      :: nfile
 INTEGER, INTENT(IN)                      :: psdl
 INTEGER, INTENT(IN OUT)                  :: dit
 INTEGER, INTENT(OUT)                     :: icoup
 INTEGER, INTENT(OUT)                     :: nfreq
 INTEGER, INTENT(OUT)                     :: npsdl
 INTEGER, INTENT(OUT)                     :: ntau
 INTEGER, INTENT(OUT)                     :: ltab
 INTEGER, INTENT(IN OUT)                  :: casecc
 INTEGER, INTENT(IN)                      :: xycdb
 INTEGER :: itlist(7), FILE,sysbuf, NAME(2),  ipsdl(6)
 REAL :: z(1)
 
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ iz(1)
 
 EQUIVALENCE (z(1),iz(1))
 
 DATA NAME  /4HRAND,1H7/
 DATA itlist/2,  55,25,1,  56,26,5 /
! *****
!     IDENTIFICATION OF VARIABLES
! *****
!     IFILE    ARRAY OF INPUT FILES
!     NFILE    LENGTH OF IFILE ARRAY
!     PSDL     POWER SPECTRAL DENSITY LISTS FROM DPD
!     DIT      DIRECT INPUT TABLES
!     ICOUP    COUPLED,UNCOUPLED, OR NOGO FLAG
!     NFREQ    NUMBER OF FRENQUICIES
!     NPSDL    NUMBER OF PSDL  SETS
!     NTAU     NUMBER OF TAUS
!     LTAB     LENGTH OF DATA FOR TAB ROUTINE
!     CASECC   CASECONTROL FILE
!     SYSBUF   LENGTH OF ONE GINO BUFFER
!     NTABL    NUMBER OF UNIQUE TABLE ID-S
!     ITABL    POINTER TO LIST OF TABLE ID-S
 
 
!     BUILD FREQUENCY LIST
 
 icoup = 0
 lcore = korsz(iz)
 ibuf1 = lcore-sysbuf
 
!     XYCDB MUST BE PRESENT
 
 FILE = xycdb
 CALL OPEN (*700,xycdb,iz(ibuf1),0)
 CALL CLOSE (xycdb,1)
 lcore = ibuf1-1
 
!     EXTRACT  SET NO FROM CASECC
 
 CALL gopen (casecc,iz(ibuf1),0)
 CALL fread (casecc,iz,166,1)
 i163  = 163
 irand = iz(i163)
 CALL CLOSE (casecc,1)
 IF (irand == 0) GO TO 700
 
!     FIND DATA FILE
 
 DO  i = 1,nfile
   FILE = ifile(i)
   CALL OPEN (*10,FILE,iz(ibuf1),0)
   CALL skprec (FILE,1)
   CALL fread (FILE,iz,10,1)
   i10 = 10
   LEN = iz(i10)-1
   nfreq = 0
   
!     EXTRACT FREQUENCIES
   
   5 CALL READ (*910,*30,FILE,f,1,0,j)
   CALL fread (FILE,iz,-LEN,0)
   nfreq = nfreq +1
   z(nfreq) = f
   GO TO  5
   30 CALL CLOSE (FILE,1)
   GO TO 40
 END DO
 
!     NO DATA FILES--EXIT
 
 700 icoup = -1
 GO TO 90
 
!     BRING IN PSDL CARDS
 
 40 lcore = lcore -nfreq
 FILE  = psdl
 CALL OPEN (*700,psdl,z(ibuf1),0)
 l = nfreq+1
 npsdl = 0
 itau  =-1
 CALL READ (*910,*41,psdl,iz(nfreq+1),lcore,0,j)
 GO TO 980
 41 k = nfreq +3
 IF (j == 2) GO TO 45
 j = k+j-1
 
!     DETERMINE RECORD THAT RANDOM TAU-S ARE IN
 
 DO  i = k,j
   IF (iz(i) == irand) GO TO 43
 END DO
 itau = -1
 GO TO 45
 
!     FOUND RANDT CARDS
 
 43 itau = i-k
 
!     FIND SELECTED PSDL CARDS
 
 45 CALL READ (*910,*47,psdl,ipsdl(1),6,0,j)
 IF (ipsdl(1) /= irand) GO TO 45
 npsdl   = npsdl+1
 iz(l  ) = ipsdl(2)
 iz(l+1) = ipsdl(3)
 iz(l+2) = ipsdl(4)
 iz(l+3) = ipsdl(5)
 iz(l+4) = ipsdl(6)
 l = l+5
 GO TO 45
 47 IF (npsdl /= 0) GO TO 48
 
!     UNABLE TO FIND SELECTED PSDL CARDS
 
 CALL CLOSE (psdl,1)
 GO TO 700
 
!     POSITION TAPE FOR TAUS
 
 48 IF (itau <= 0) GO TO 49
 CALL skprec (psdl,itau)
 49 lcore = lcore-npsdl*5
 
!     EXTRACT LIST OF TABLES  AND CHECK FOR COUPLED SYSTEM
 
 jj = nfreq +1
 k  = nfreq +5*npsdl
 ntabl = 0
 itabl = ibuf1-1
 loop60:  DO   i = jj,k,5
   IF (iz(i) == iz(i+1)) GO TO 61
   
!     COUPLED
   
   icoup =1
   61 IF (ntabl == 0) GO TO 62
   DO  j=1,ntabl
     l = itabl +j
     IF (iz(l) == iz(i+4)) CYCLE loop60
   END DO
   
!     STORE TABLE ID
   
   62 ntabl = ntabl +1
   iz(itabl) = iz(i+4)
   itabl = itabl -1
 END DO loop60
 iz(itabl) = ntabl
 
!     BRING IN  TAU-S
 
 ntau  = 0
 lcore = lcore- ntabl-1
 IF(itau == -1) GO TO 70
 CALL READ (*70,*70,psdl,z(k+1),lcore,0,ntau)
 GO TO 980
 70 CALL CLOSE (psdl,1)
 
!     SETUP FOR TABLES
 
 lcore = lcore -ntau
 ltab  = 0
 IF(ntabl == 0) GO TO 90
 l =  k +ntau+1
 CALL pretab (dit,iz(l),z(l),iz(ibuf1),lcore,ltab,iz(itabl),itlist (1))
 90 RETURN
 
!     FILE ERRORS
 
 901 CALL mesage (ip1,FILE,NAME)
 910 ip1 =-2
 GO TO 901
 980 ip1= -8
 GO TO 901
END SUBROUTINE rand7
