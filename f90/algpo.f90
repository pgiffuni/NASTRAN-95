SUBROUTINE algpo (scr1)
     
 
 INTEGER, INTENT(IN)                      :: scr1
 EXTERNAL        orf
 INTEGER :: apress,atemp,strml,NAME(2),corwds,itrl(7),two(32),  &
     orf,casecc,caseca,geom3a, lend(3),sysbuf,  &
     pgeom,rd,rdrew,wrt,wrtrew,rew,norew,pload2(3), temp(3),tempd(3),lrec(5)
 DIMENSION       rrec(5)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / apress,atemp,strml,pgeom,iprtk,ifail
 COMMON /system/ sysbuf,nout
 COMMON /names / rd,rdrew,wrt,wrtrew,rew,norew
 COMMON /two   / two
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE     (lrec(1),rrec(1))
 DATA    lend  / 3*2147483647 /
 DATA    NAME  / 4HALG ,4H    /
 DATA    labp  / 4HPLOA/      , labt /4HTEMP/
 DATA    pload2/ 6809,68,199 /, temp/5701,57,27/, tempd/5641,65,98/
 DATA    casecc, caseca,geom3a /101,201,202/
 
!     ALG WILL USE OPEN CORE AT IZ
!     ALLOCATE OPEN CORE
 
 nz    = korsz(iz)
 ibuf1 = nz - sysbuf
 ibuf2 = ibuf1 - sysbuf - 1
 last  = ibuf2 - 1
 
!     CHECK FOR SUFFICIENT CORE
 
 IF (last <= 0) CALL mesage (-8,0,NAME)
 left   = corwds(iz(1),iz(last))
 kaperr = 0
 katerr = 0
 ifail  = 1
 
!     OPEN GEOM3A FOR OUTPUT OF PLOAD2 AND TEMP DATA
 
 CALL gopen (geom3a,iz(ibuf1),wrtrew)
 
!     AERODYNAMIC PRESSURE SECTION
 
 IF (apress < 0) GO TO 20
 ifile = scr1
 CALL OPEN (*901,scr1,iz(ibuf2),rdrew)
 8 CALL READ (*11,*10,scr1,lrec,5,1,nwar)
 10 IF (lrec(1) == labp) GO TO 12
 GO TO 8
 
!     NO PLOAD2 CARDS ON SCR1 FILE
 
 11 kaperr = 1
 idp    = 0
 CALL REWIND (scr1)
 WRITE (nout,2001) uwm
 GO TO 20
 12 CONTINUE
 
!     CREATE PLOAD2 RECORD
 
 CALL WRITE (geom3a,pload2,3,0)
 idp = lrec(3)
 GO TO 16
 14 CALL READ (*18,*15,scr1,lrec,5,1,nwar)
 15 IF (lrec(1) /= labp) GO TO 18
 16 CALL WRITE (geom3a,lrec(3),3,0)
 GO TO 14
 18 CALL WRITE (geom3a,iz,0,1)
 CALL REWIND (scr1)
 
!     AERODYNAMIC TEMPERATURE SECTION
 
 20 IF (atemp < 0) GO TO 35
 IF (apress < 0) CALL OPEN (*901,scr1,iz(ibuf2),rdrew)
 21 CALL READ (*23,*22,scr1,lrec,5,1,nwar)
 22 IF (lrec(1) == labt) GO TO 24
 GO TO 21
 
!     NO TEMP CARDS ON SCR1 FILE
 
 23 katerr = 1
 idt = 0
 WRITE (nout,2002) uwm
 GO TO 35
 24 CONTINUE
 
!     CREATE TEMP RECORD
 
 CALL WRITE (geom3a,temp,3,0)
 idt   = lrec(3)
 dtemp = rrec(5)
 itpd  = 1
 GO TO 28
 26 CALL READ (*30,*27,scr1,lrec,5,1,nwar)
 27 IF (lrec(1) /= labt) GO TO 30
 28 CALL WRITE (geom3a,lrec(3),3,0)
 itpd = itpd + 1
 IF (itpd <= 3) dtemp = dtemp + rrec(5)
 GO TO 26
 30 CALL WRITE (geom3a,iz,0,1)
 
!     CREATE TEMPD RECORD. AVERAGE FIRST THREE TEMPS. ON BLADE ROOT.
 
 CALL WRITE (geom3a,tempd,3,0)
 CALL WRITE (geom3a,idt,1,0)
 dtemp = dtemp/3.0
 CALL WRITE (geom3a,dtemp,1,1)
 
!     CLOSE GEOM3A
 
 35 CALL WRITE (geom3a,lend,3,1)
 CALL CLOSE (geom3a,1)
 itrl(1) = geom3a
 itrl(2) = 0
 itrl(3) = 0
 itrl(4) = 0
 itrl(5) = 0
 itrl(6) = 0
 itrl(7) = 0
 IF (apress < 0 .OR. kaperr == 1) GO TO 40
 ibit = 68
 i1 = (ibit-1)/16 + 2
 i2 = ibit - (i1-2)*16 + 16
 itrl(i1) = orf(itrl(i1),two(i2))
 40 IF (atemp < 0 .OR. katerr == 1) GO TO 50
 ibit = 57
 i1 = (ibit-1)/16 + 2
 i2 = ibit - (i1-2)*16 + 16
 itrl(i1) = orf(itrl(i1),two(i2))
 ibit = 65
 i1 = (ibit-1)/16 + 2
 i2 = ibit - (i1-2)*16 + 16
 itrl(i1) = orf(itrl(i1),two(i2))
 50 CALL wrttrl (itrl)
 
!     CLOSE SCR1
 
 IF (apress >= 0 .OR. atemp >= 0) CALL CLOSE (scr1,1)
 IF (kaperr == 1) apress = -1
 IF (katerr == 1) atemp  = -1
 
!     SET IFAIL TO INDICATE ALG MODULE FAILED. CONDITIONAL JUMP BASED
!     ON VALUE OF IFAIL IS PERFORMED AFTER EXITING FROM ALG MODULE.
 
 IF (apress == -1 .AND. atemp == -1) ifail = -1
 
!     NEW CASE CONTROL DATA BLOCK
!     OPEN CASECC AND COPY ALL SUBCASES WITH CHANGES MADE TO
!     STATIC AND THERMAL LOAD ID-S
 
 ifile = casecc
 CALL OPEN (*901,casecc,iz(ibuf1),rdrew)
 CALL fwdrec (*902,casecc)
 CALL gopen (caseca,iz(ibuf2),wrtrew)
 60 CALL READ (*70,*65,casecc,iz,left,1,nwds)
 65 izx = 4
 iz(izx) = idp
 izx = 7
 iz(izx) = idt
 CALL WRITE (caseca,iz,nwds,1)
 GO TO 60
 70 CALL CLOSE (casecc,1)
 CALL CLOSE (caseca,1)
 itrl(1) = casecc
 CALL rdtrl (itrl)
 itrl(1) = caseca
 CALL wrttrl (itrl)
 GO TO 999
 901 CALL mesage (-1,ifile,NAME)
 GO TO 999
 902 CALL mesage (-2,ifile,NAME)
 GO TO 999
 999 RETURN
 
 2001 FORMAT (a25,' - ALG MODULE - AERODYNAMIC PRESSURES REQUESTED VIA',  &
     ' PARAM APRESS, BUT NOUT3=0 IN AERODYNAMIC INPUT', /41X,  &
     'OR AERODYNAMIC CALCULATION FAILED. REQUEST IGNORED.')
 2002 FORMAT (a25,' - ALG MODULE - AERODYNAMIC TEMPERATURES REQUESTED ',  &
     'VIA PARAM ATEMP, BUT NOUT3=0 IN AERODYNAMIC INPUT' ,/41X,  &
     'OR AERODYNAMIC CALCULATION FAILED. REQUEST IGNORED.')
END SUBROUTINE algpo
