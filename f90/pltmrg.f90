SUBROUTINE pltmrg
     
!     MODULE PLTMRG WRITES GINO DATA BLOCKS WHICH ARE USED AS INPUT TO
!     THE PLOT MODULE FOR PLOTTING A SUBSTRUCTURE.
 
!     APRIL 1974
 
 LOGICAL :: ident
 INTEGER :: buf      ,sysbuf   ,z(3)     ,casess   ,pcdb    ,  &
     pltp     ,gps      ,els      ,bgp      ,casep   ,  &
     eqex     ,scr1     ,srd      ,plts     ,FILE    ,  &
     eqss     ,subr(2)  ,casecc(2),buf1     ,elid    ,  &
     buf2     ,buf3     ,buf4     ,buf5     ,rc      ,  &
     bar      ,quad4    ,tria3    ,offset
 REAL :: rz
 COMMON /BLANK /   NAME(2)  ,ngptot   ,lsil     ,npset   , nm(2)    ,buf(7)
 COMMON /system/   sysbuf
 COMMON /names /   rd       ,rdrew    ,wrt      ,wrtrew  , rew      ,norew
 COMMON /zzzzzz/   rz(1)
 EQUIVALENCE       (z(1),rz(1))
 DATA    plts  ,   eqss     ,subr               ,casecc          /  &
     4HPLTS,   4HEQSS   ,4HPLTM   ,4HRG     ,4HCASE   ,4HCC  /
 DATA    casess,   pcdb     ,pltp     ,gps      ,els      /  &
     101   ,   102      ,201      ,202      ,203      /,  &
     bgp   ,   casep    ,eqex     ,scr1     ,srd      /  &
     204   ,   205      ,206      ,301      ,1        /,  &
     bar   ,   quad4    ,tria3                        /  &
     2HBR  ,   2HQ4     ,2HT3                         /
 
!     INITIALIZE
 
 ncore = korsz(z)
 buf1  = ncore- sysbuf + 1
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf4  = buf3 - sysbuf
 buf5  = buf4 - sysbuf
 ncore = buf5 - 1
 ngptot= 0
 lsil  = 0
 npset =-1
 IF (ncore <= 0) GO TO 9008
 CALL sofopn (z(buf3),z(buf4),z(buf5))
 
!     STRIP SUBSTRUCTURE RECORDS FROM CASESS AND WRITE CASEP (CASECC)
 
 FILE = casess
 CALL OPEN (*9001,casess,z(buf1),rdrew)
 FILE = casep
 CALL OPEN  (*9001,casep,z(buf2),wrtrew)
 CALL fname (casep,buf)
 CALL WRITE (casep,buf,2,1)
 FILE = casess
 10 CALL READ (*9002,*9003,casess,z,2,1,nwds)
 IF (z(1) /= casecc(1) .OR. z(2) /= casecc(2)) GO TO 10
 20 CALL READ (*40,*30,casess,z,ncore,1,nwds)
 GO TO 9008
 30 CALL WRITE (casep,z,nwds,1)
 GO TO 20
 40 CALL clstab (casep,rew)
 CALL CLOSE  (casess,rew)
 
!     BASIC GRID POINT DATA
 
 nm(1) = NAME(1)
 nm(2) = NAME(2)
 item  = plts
 CALL sfetch (NAME,plts,srd,rc)
 IF (rc /= 1) GO TO 6100
 
!     READ SUBSTRUCTURE NAMES AND TRANSFORMATION DATA INTO OPEN CORE.
 
 CALL suread (z,3,nwds,rc)
 IF (rc /= 1) GO TO 6106
 nss = z(3)
 IF (14*nss > ncore) GO TO 9008
 CALL suread (z,14*nss,nwds,rc)
 IF (rc /= 1) GO TO 6106
 icore = 14*nss + 1
 
!     READ THE BASIC GRID POINT DATA FROM THE PLTS ITEM OF EACH BASIC
!     SUBSTRUCTURE COMPRISING THE PSEUDOSTRUCTURE TO BE PLOTTED.
!     TRANSFORM THE COORDINATES TO THE BASIC COORDINATE SYSTEM OF THE
!     PSEUDOSTRUCTURE AND WRITE THEM ON BGP (BGPDT).
 
 FILE = bgp
 CALL OPEN (*9001,bgp,z(buf1),wrtrew)
 CALL fname (bgp,buf)
 CALL WRITE (bgp,buf,2,1)
 j = 1
 120 nm(1) = z(j  )
 nm(2) = z(j+1)
 ngp = 0
 CALL sfetch (nm,plts,srd,rc)
 IF (rc == 1) GO TO 130
 CALL smsg (rc-2,plts,nm)
 GO TO 170
 130 i = 1
 CALL sjump (i)
 ident = .false.
 DO  i = 1,3
   IF (z(j+i+1) /= 0) GO TO 150
   IF (z(j+i+5) /= 0) GO TO 150
   IF (z(j+i+9) /= 0) GO TO 150
   IF (ABS(rz(j+4*i+1)-1.0) > 1.0E-4) GO TO 150
 END DO
 ident = .true.
 150 CALL suread (buf,4,nwds,rc)
 IF (rc == 2) GO TO 170
 ngp = ngp + 1
 IF (ident .OR. buf(1) < 0) GO TO 160
 buf(5) = z(j+2)
 buf(6) = z(j+3)
 buf(7) = z(j+4)
 CALL gmmats (z(j+5),3,3,-2,buf(2),3,1,0,buf(5))
 CALL WRITE (bgp,buf,1,0)
 CALL WRITE (bgp,buf(5),3,0)
 GO TO 150
 160 CALL WRITE (bgp,buf,4,0)
 GO TO 150
 170 ngptot = ngptot+ngp
 z(j+2) = ngp
 j = j + 14
 IF (j < icore) GO TO 120
 CALL WRITE (bgp,0,0,1)
 CALL CLOSE (bgp,rew)
 buf(1) = bgp
 buf(2) = ngptot
 DO  i = 3,7
   buf(i) = 0
 END DO
 CALL wrttrl (buf)
 
!     ALLOCATE 5 WORDS PER COMPONENT BASIC SUBSTRUCTURE AT THE TOP OF
!     OPEN CORE.  THIS ARRAY IS HEREINAFTER REFERRED TO AS *SDATA*
 
!     SAVE THE BASIC SUBSTRUCTURE NAMES AND THE NUMBER OF STRUCTURAL
!     GRID POINTS IN EACH IN SDATA.  DO NOT SAVE SUBSTRUCTURES FOR
!     WHICH NO PLTS ITEM WAS FOUND.
 
 j = 1
 DO  i = 1,nss
   IF (z(14*i-11) == 0) CYCLE
   z(j  ) = z(14*i-13)
   z(j+1) = z(14*i-12)
   z(j+2) = z(14*i-11)
   j = j + 5
 END DO
 IF (j <= 1) GO TO 9200
 nss = j/5
 isx = nss*5
 icore = j
 lcore = ncore - j + 1
 
!     COMPUTE EQEX (EQEXIN)
 
 
!     READ THE EQEXIN DATA FROM THE PLTS ITEM OF EACH BASIC SUBSTRUCTURE
!     USE THREE WORDS IN OPEN CORE FOR EACH GRID POINT   (1) EXTERNAL
!     ID, (2) INTERNAL ID, (3) SUBSTRUCTURE SEQUENCE NUMBER IN SDATA.
!     INCREMENT THE INTERNAL IDS BY THE NUMBER OF GRID POINTS ON THE
!     PRECEDING SUBSTRUCTURES.
 
 k   = icore
 ngp = 0
 DO  i = 1,nss
   nm(1) = z(5*i-4)
   nm(2) = z(5*i-3)
   CALL sfetch (nm,plts,srd,rc)
   n = 2
   CALL sjump (n)
   rc = 3
   IF (n < 0) GO TO 6106
   n = z(5*i-2)
   DO  j = 1,n
     CALL suread (z(k),2,nwds,rc)
     IF (rc /= 1) GO TO 6106
     z(k+1) = z(k+1) + ngp
     z(k+2) = i
     k = k + 3
     IF (k+2 > ncore) GO TO 9008
   END DO
   ngp = ngp + n
 END DO
 
!     SORT ON EXTERNAL IDS AND WRITE RECORD 1 OF EQEX.
 
 CALL sort (0,0,3,1,z(icore),3*ngp)
 FILE = eqex
 CALL OPEN (*9001,eqex,z(buf1),wrtrew)
 CALL fname (eqex,buf)
 CALL WRITE (eqex,buf,2,1)
 DO  i = 1,ngp
   CALL WRITE (eqex,z(icore+3*i-3),2,0)
 END DO
 CALL WRITE (eqex,0,0,1)
 
!     SAVE THE TABLE IN OPEN CORE ON SCR1 TO USE IN COMPUTING RECORD 2
!     OF EQEX
 
 FILE = scr1
 CALL OPEN  (*9001,scr1,z(buf2),wrtrew)
 CALL WRITE (scr1,z(icore),3*ngp,1)
 CALL CLOSE (scr1,rew)
 CALL OPEN  (*9001,scr1,z(buf2),rdrew)
 
!     READ GROUP 0 OF THE EQSS ITEM OF THE SUBSTRUCTURE TO BE PLOTTED
!     INTO OPEN CORE AT ICORE.  READ THE EXTERNAL AND INTERNAL IDS FOR
!     EACH CONTRIBUTING BASIC SUBSTRUCTURE INTO OPEN CORE FOLLOWING
!     GROUP 0.  SAVE THE CORE POINTERS FOR EACH GROUP IN SDATA.
 
 nm(1) = NAME(1)
 nm(2) = NAME(2)
 item  = eqss
 CALL sfetch (NAME,eqss,srd,rc)
 IF (rc /= 1) GO TO 6100
 CALL suread (z(icore),lcore,nwds,rc)
 IF (rc /= 2) GO TO 9008
 k   = icore + nwds
 n   = z(icore+2)
 iss = 1
 DO  i = 1,n
   IF (iss > isx) GO TO 240
   IF (z(icore+2*i+2) /= z(iss) .OR. z(icore+2*i+3) /= z(iss+1)) GO TO 240
   z(iss+3) = k
   230 IF (k+2 > ncore) GO TO 9008
   CALL suread (z(k),3,nwds,rc)
   k = k + 2
   IF (rc == 1) GO TO 230
   z(iss+4) = (k-z(iss+3))/2
   iss = iss + 5
   CYCLE
   240 j = 1
   CALL sjump (j)
 END DO
 
!     READ SIL NUMBERS INTO OPEN CORE.
 
 ksil = k - 1
 n = z(icore+3)
 IF (ksil+n+1 > ncore) GO TO 9008
 DO  i = 1,n
   CALL suread (z(ksil+i),2,nwds,rc)
   IF (rc /= 1) GO TO 6106
 END DO
 lsil = z(ksil+n)
 
!     READ THE TABLE OF EXTERNAL ID (GP), INTERNAL ID (IP), AND SUB-
!     STRUCTURE NUMBER (SSN) FROM SCR1 ONE ENTRY AT A TIME.  LOCATE
!     THE GP IN THE EQSS DATA INDICATED BY SSN AND LOOK UP THE SIL
!     NUMBER.  WRITE GP AND SIL ON EQEX.  IF GP NOT FOUND, THEN SIL=-1.
 
 270 CALL READ (*9002,*290,scr1,buf,3,0,n)
 i  = buf(3)
 j  = z(5*i-1)
 i5 = 5*i
 CALL bisloc (*280,buf(1),z(j),2,z(i5),k)
 i = z(j+k) + ksil
 buf(2) = 10*z(i) + 1
 CALL WRITE (eqex,buf,2,0)
 GO TO 270
 280 buf(2) = -1
 CALL WRITE (eqex,buf,2,0)
 GO TO 270
 290 CALL WRITE (eqex,0,0,1)
 CALL CLOSE (eqex,rew)
 CALL CLOSE (scr1,rew)
 buf(1) = eqex
 buf(2) = ngptot
 DO  i = 3,7
   buf(i) = 0
 END DO
 CALL wrttrl (buf)
 
!     INTERPRET PLOT SETS AND GENERATE PLTP (PLTPAR)
 
 
!     AT PRESENT, ONLY ONE PLOT SET (DEFINED IN PHASE 1) IS ALLOWED.
 
!     PHASE 2 PLOT SET DEFINITIONS ARE IGNORED.
 
!     COPY PCDB TO PLTP
 
 FILE = pcdb
 CALL OPEN (*9001,pcdb,z(buf1),rdrew)
 CALL fwdrec (*9002,pcdb)
 FILE = pltp
 CALL OPEN  (*9001,pltp,z(buf2),wrtrew)
 CALL fname (pltp,buf)
 CALL WRITE (pltp,buf,2,1)
 310 CALL READ  (*330,*320,pcdb,z(icore),lcore,1,nwds)
 GO TO 9008
 320 CALL WRITE (pltp,z(icore),nwds,1)
 GO TO 310
 330 CALL CLOSE (pcdb,rew)
 CALL CLOSE (pltp,rew)
 buf(1) = pcdb
 CALL rdtrl (buf)
 buf(1) = pltp
 CALL wrttrl (buf)
 DO  i = 1,nss
   z(5*i-1) = 0
   z(5*i  ) = 1
 END DO
 npset = 1
 
!     GPSETS
 
 
!     LOCATE THE GPSETS DATA OF THE PLTS ITEM OF EACH BASIC SUBSTRUCTURE
!     AND READ THE NUMBER OF GRID POINTS IN THE ELEMENT SET.  STORE THIS
!     AS THE FOURTH ENTRY IN SDATA
 
 n = 3
 ngpset = 0
 item   = plts
 DO  i = 1,nss
   nm(1) = z(5*i-4)
   nm(2) = z(5*i-3)
   CALL sfetch (nm,plts,srd,rc)
   CALL sjump (n)
   rc = 3
   IF (n < 0) GO TO 6106
   CALL suread (z(5*i-1),1,nwds,rc)
   IF (rc /= 1) GO TO 6106
   ngpset = ngpset + z(5*i-1)
 END DO
 
!     WRITE RECORDS 0 AND 1 OF GPS AND FIRST WORD OF RECORD 2.
 
 FILE = gps
 CALL OPEN  (*9001,gps,z(buf1),wrtrew)
 CALL fname (gps,buf)
 CALL WRITE (gps,buf,2,1)
 CALL WRITE (gps,1,1,1)
 CALL WRITE (gps,ngpset,1,0)
 
!     READ GPSETS DATA FROM THE PLTS ITEM OF EACH BASIC SUBSTRUCTURE.
!     INCREMENT THE ABSOLUTE VALUE OF THE POINTERS IN IT BY THE NUMBER
!     OF GRID POINTS IN THE ELEMENT SETS OF THE PRECEDING BASIC
!     SUBSTRUCTURES.  WRITE THE RESULT ON GPS (GPSETS).
 
 n = 3
 ngpset = 0
 DO  i = 1,nss
   CALL sfetch (z(5*i-4),plts,srd,rc)
   CALL sjump (n)
   CALL suread (z(icore),lcore,nwds,rc)
   IF (rc /= 2) GO TO 9008
   nwds = nwds - 1
   DO  j = 1,nwds
     IF (z(icore+j) < 0.0) THEN
       GO TO  1020
     ELSE IF (z(icore+j) == 0.0) THEN
       GO TO  1040
     ELSE
       GO TO  1030
     END IF
     1020 z(icore+j) = z(icore+j) - ngpset
     CYCLE
     1030 z(icore+j) = z(icore+j) + ngpset
     1040 CONTINUE
   END DO
   CALL WRITE (gps,z(icore+1),nwds,0)
   ngpset = ngpset + z(5*i-1)
 END DO
 CALL clstab (gps,rew)
 
!     ELSETS
 
 
!     READ THE ELSETS DATA FROM THE PLTS ITEM OF EACH BASIC SUBSTRUCTURE
!     INCREMENT ALL NON-ZERO GRID POINT CONNECTION INDICES BY THE NUMBER
!     OF STRUCTURAL GRID POINTS OF THE PRECEDING SUBSTRUCTURES.  WRITE
!     THE RESULT ON ELS (ELSETS).
 
!     NOTE   THE ELEMENT TYPES WILL BE SCRAMBLED.  LIKE ELEMENT TYPES
!            FROM THE CONTRIBUTING BASIC SUBSTRUCTURES WILL NOT BE
!            GROUPED TOGETHER.
 
!     NOTE   THE BAR HAS ADDITIONALLY 6 OFFSET DATA VALUES. QUAD4 AND
!            TRIA3 HAS 1 OFFSET DATA EACH
 
 FILE = els
 CALL OPEN  (*9001,els,z(buf1),wrtrew)
 CALL fname (els,buf)
 CALL WRITE (els,buf,2,1)
 ngp = 0
 
!     LOOP OVER BASIC SUBSTRUCTURES
 
 DO  i = 1,nss
   nm(1) = z(5*i-4)
   nm(2) = z(5*i-3)
   CALL sfetch (nm,plts,srd,rc)
   n  = 4
   CALL sjump (n)
   rc = 3
   IF (n < 0) GO TO 6106
   
!     LOOP OVER ELEMENT TYPES
   
   2010 CALL suread (buf,2,n,rc)
   IF (rc == 2) GO TO 2040
   IF (rc /= 1) GO TO 6106
   CALL WRITE (els,buf,2,0)
   ngpel  = buf(2)
   offset = 0
   IF (buf(1) == bar) offset = 6
   IF (buf(1) == quad4 .OR. buf(1) == tria3) offset = 1
   
!     LOOP OVER ELEMENTS
   
   2020 CALL suread (elid,1,n,rc)
   IF (rc /= 1) GO TO 6106
   CALL WRITE (els,elid,1,0)
   IF (elid <= 0) GO TO 2010
   CALL suread (indx,1,n,rc)
   CALL WRITE (els,indx,1,0)
   CALL suread (z(icore),ngpel+offset,n,rc)
   IF (rc /= 1) GO TO 6106
   
!     LOOP OVER CONNECTIONS
   
   k = icore
   DO  j = 1,ngpel
     IF (z(k) /= 0) z(k) = z(k) + ngp
     k = k + 1
   END DO
   CALL WRITE (els,z(icore),ngpel+offset,0)
   GO TO 2020
   2040 ngp = ngp + z(5*i-2)
 END DO
 
 CALL WRITE  (els,0,0,1)
 CALL clstab (els,rew)
 
!     NORMAL MODULE COMPLETION
 
 CALL sofcls
 RETURN
 
!     ABNORMAL MODULE COMPLETION
 
 6100 IF (rc == 2) rc = 3
 CALL smsg (rc-2,item,nm)
 GO TO 9200
 6106 CALL smsg (rc+4,item,nm)
 GO TO 9200
 9001 n = 1
 GO TO 9100
 9002 n = 2
 GO TO 9100
 9003 n = 3
 GO TO 9100
 9008 n = 8
 9100 CALL mesage (n,FILE,subr)
 CALL CLOSE  (FILE,rew)
 9200 CALL sofcls
 npset = -1
 RETURN
END SUBROUTINE pltmrg
