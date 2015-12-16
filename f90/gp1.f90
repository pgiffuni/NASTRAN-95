SUBROUTINE gp1
     
!     GP1  BUILDS THE FOLLOWING DATA BLOCKS--
!       1. GRID POINT LIST (GPL)
!       2. EXTERNAL INTERNAL GRID POINT EQUIVALENCE TABLE (EQEXIN)
!       3. GRID POINT DEFINITION TABLE (GPDT)
!       4. COORDINATE SYSTEM TRANSFORMATION MATRICES (CSTM)
!       5. BASIC GRID POINT DEFINITION TABLE (BGPDT)
!       6. SCALAR INDEX LIST (SIL)
 
!     THE FOLLOWING CARDS ARE READ BY GP1--
!       1. GRID
!       2. CELASI, CDAMPI, CMASSI  (I=1,2,3,4)
!       3. SPOINT
!       4. SEQGP   (SEQEP IS PROCESSED IN DPD1)
!       5. CORDIJ  (I=1,2,  J=R,S,C)
 
!     IMPORTANT
!     =========
!     REVISED  7/89 BY G.CHAN/UNISYS, TO ALLOW GRID, SCALAR AND EXTRA
!     POINT EXTERNAL ID UP TO 8 DIGITS FOR ALL 32-BIT MACHINES
!     PREVIOUSLY, ID OF 2000000 IS THE UPPER LIMIT FOR IBM AND VAX
 
!     REVISED  8/89 BY G.CHAN/UNISYS, AS PART OF THE EFFORT TO ALLOW A
!     NASTRAN JOB TO EXCEED 65535 LIMIT.
!     NORMALLY, IF GRID POINTS OR SCALAR POINTS DO NOT HAVE VERY LARGE
!     EXTERNAL ID NUMBERS, THEIR ID NOS. ARE MULTIPLIED BY 1000, SO THAT
!     999 ADDITIONAL POINTS CAN SQUEEZE IN VIA SEQGP CARDS. (NOTE - A
!     7- OR 8-DIGIT ID NO., TIMES 1000, EXCEEDS A 32-BIT WORD COMPUTER
!     HARDWARE LIMIT). THIS MULTIPLY FACTOR IS NOW ADJUSTABLE, 1000,100,
!     OR 10, SO THAT ADDITIONAL DIGITS CAN BE USED FOR THE EXTERNAL GRID
!     OR SCALAR POINTS IN CASE THERE ARE LIMITTED SEQGP CARDS PRESENT.
!     THIS VARIABLE MULTIPLIER (10,100, OR 1000) IS ALSO RECORDED IN THE
!     3RD WORD OF THE HEADER RECORD OF THE GPL DATA BLOCK FOR LATER USE.
!     THE ACTUAL FACTOR OF THE MULTIPLIER IS ALSO MACHINE DEPENDENT.
!     UNIVAC, A 36-BIT MACHINE, CAN HAVE A MULTIPLIER OF 100 OR 1000.
!     OTHER 60- OR 64- BIT MACHINES, THE MULTIPLIER REMAINS AT 1000
!     IF THE MULTIPLIER IS 1000, THE SEQGP AND SEQEP CARDS, AS BEFORE,
!     CAN HAVE 4 SEQID LEVELS, SUCH AS XXX.X.X.X
!     IF THE MULTIPLIER IS 100, SEQGP AND SEQEP CARDS ARE LIMITED TO
!     3 SEQID LEVELS, XXX.X.X
!     FINALLY, IF MULTIPLIER IS 10, SEQGP AND SEQEP ARE LIMITED TO XXX.X
 
!     SPECIAL CONSIDERATION FOR THE AXISYM. AND HYDROELAS. PROBLEMS - 10
!     IS USED FOR THE MULTIPLIER, AND THEREFOR A ONE SEQID LEVEL IS
!     AVAILABLE. PREVIOUSLY, SEQGP CARDS WERE NOT USED IN AXISYM. AND
!     HYDROELAS. PROBLEMS, AND NO USER WARNING MESSAGE PRINTED
 
!     NO ADJUSTABLE MULTIPLY FACTOR FOR SUBRSTRUCTURING (MULT=1000,
!     SEE ALSO SGEN)
 
!     THE 65535 LIMITATION INVOLVES ONLY A SAMLL CHANGE IN STA 973
 
 EXTERNAL        rshift
 INTEGER :: rd,wrt,cls,FILE,elem,axic,z,sysbuf,buf1,buf2,  &
     buf3,geomp,gpl,eqexin,gpdt,cstm,bgpdt,sil,scr1,  &
     scr2,wrtrew,rdrew,a,spoint,flag,grid,clsrew,  &
     seqgp,gpfl,cord,cordij,gp1ah,geom1,geom2,ptr,  &
     solv,solvp,scalpt,TYPE,offset,rshift
 REAL :: length
 DIMENSION       a(34),aa(34),ab(3),ac(3),ai(3),aj(3),ak(3),ax(3),  &
     ar(3),spoint(2),grid(2),seqgp(2),cordij(12),  &
     cord(6),gp1ah(2),scalpt(2),zz(1),mcb(7)
 CHARACTER (LEN=29) :: lvl1,lvl2
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /BLANK / luset,nogpdt,nocstm
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /zzzzzz/ z(1)
 COMMON /system/ ksystm(100)
 COMMON /setup / nfile(6),ptr
 COMMON /gpta1 / nelem,lastx,incrx,elem(1)
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 EQUIVALENCE     (ksystm( 1),sysbuf), (ksystm( 2),iout ),  &
     (ksystm(24),icfiat), (ksystm(27),axic ),  &
     (ksystm(38),iaxif ), (ksystm(40),nbpw ),  &
     (ksystm(56),itherm), (ksystm(69),isubs)
 EQUIVALENCE     (z( 1),zz(1)), (a( 1),aa(1)), (a( 4),ab(1)),  &
     (a( 7),ac(1)), (a(10),ai(1)), (a(13),aj(1)),  &
     (a(16),ak(1)), (a(19),ax(1)), (a(22),ar(1)),  &
     (nocstm,ifl ), (geomp,geom1), (mcb(2),kn  )
 EQUIVALENCE     (igpdt,icsdt)
 DATA    geom1 / 101/, geom2 / 102/,  &
     gpl   / 201/, eqexin/ 202/, gpdt  / 203/,  &
     cstm  / 204/, bgpdt / 205/, sil   / 206/, scr1  / 301/, scr2  / 302/
 DATA    gp1ah / 4HGP1 , 4H    /, cord  / 6,6,6,13,13,13/,  &
     grid  / 4501,45  /, seqgp / 5301,53  /,  &
     cordij/ 1701,17,1801,18,1901,19,2001,20,2101,21,2201,22/, scalpt/ 5551,49  /
 DATA    mcb   / 7*0      /, large / 100000000/,  &
     lvl1  / '3  I.E.  XXX.X.X.X TO XXX.X.X' /,  &
     lvl2  / '2  I.E.  XXX.X.X.X TO XXX.X  ' /
 
!     PERFORM GENERAL INITIALIZATION
 
 CALL delset
 nz     = korsz(z)
 buf1   = nz - sysbuf - 2
 buf2   = buf1 - sysbuf
 buf3   = buf2 - sysbuf
 nogo   = 0
 nocstm = 0
 nogpdt =-1
 nogmp1 = 1
 maxa1  = 0
 mult   = 1000
 axi    = 0
 IF (axic /= 0 .OR. iaxif /= 0) axi = 1
 IF (axi   /= 0) mult = 10
 IF (isubs /= 0) mult = 1000
 imax = large
 IF (nbpw == 32) imax =  2147483
 IF (nbpw == 36) imax = 34359738
!         2147483=2**31/1000   34359738=2**35/1000
 
!     READ SCALAR ELEMENT CONNECTION CARDS (IF PRESENT).
!     EXTRACT SCALAR POINTS AND WRITE THEM ON SCR2.
 
 FILE = scr2
 CALL OPEN (*1170,scr2,z(buf2),wrtrew)
 nosclr= 0
 m8    =-8
 a(11) =-1
 DO  k = 12,16
   a(k) = 0
 END DO
 CALL preloc (*80,z(buf1),geom2)
 i = 1
 DO  i = 1,lastx,incrx
   kk = elem(i+10)
   IF (kk == 0) CYCLE
   CALL locate (*60,z(buf1),elem(i+3),flag)
   nn = elem(i+5)
   40 CALL READ (*1180,*60,geom2,a,nn,0,flag)
   DO  k = 3,4
     IF (a(k) == 0 .OR. (kk == 1 .AND. a(k+2) /= 0)) CYCLE
     a(10)  = a(k)
     nosclr = 1
     CALL WRITE (scr2,a(10),1,0)
   END DO
   GO TO 40
   60 CONTINUE
 END DO
 
!     COPY SCALAR POINTS DEFINED ON SPOINT CARDS (IF PRESENT) ONTO SCR2.
 
 CALL locate (*80,z(buf1),scalpt,flag)
 nosclr = 1
 CALL READ   (*1180,*70,geom2,z,buf2-1,1,n)
 CALL mesage (m8,0,gp1ah)
 70 CALL WRITE  (scr2,z,n,0)
 
!     CLOSE FILES. IF SCALAR POINTS PRESENT, SORT LIST.
!     THEN DISCARD DUPLICATES AND WRITE UNIQUE LIST ON SCR2.
 
 80 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr2,clsrew)
 CALL CLOSE (geom2,clsrew)
 IF (nosclr == 0) GO TO 110
 nfile(1) = gpdt
 nfile(2) = bgpdt
 nfile(3) = sil
 CALL OPEN  (*1170,scr2,z(buf1),rdrew)
 CALL sorti (scr2,0,1,1,z,buf1-1)
 CALL CLOSE (scr2,clsrew)
 FILE = nfile(6)
 CALL OPEN (*1170,FILE,z(buf1),rdrew)
 CALL OPEN (*1170,scr2,z(buf2),wrtrew)
 last = -1
 90 CALL READ (*1180,*100,FILE,a(10),1,0,flag)
 IF (a(10) == last) GO TO 90
 CALL WRITE (scr2,a(10),1,0)
 last = a(10)
 GO TO 90
 100 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr2,clsrew)
 CALL CLOSE (FILE,clsrew)
 CALL OPEN  (*1170,scr2,z(buf3),rdrew)
 
!     READ GRID ENTRIES (IF PRESENT).
!     MERGE GRID AND SCALAR NOS.
!     CREATING LIST IN CORE OF EXTERNAL NO., MULT * EXTERNAL NO.
!     WRITE 7-WORD GRID AND SCALAR ENTRIES ON SCR1.
 
 110 a(1)  = large
 a(10) = large
 FILE  = scr1
 IF (maxa1 == 0) CALL OPEN (*1170,scr1,z(buf2),wrtrew)
 i = -1
 nogrid = 0
 IF (maxa1 == 0) CALL preloc (*190,z(buf1),geom1)
 CALL locate (*200,z(buf1),grid,flag)
 nogrid = 1
 CALL READ  (*1180,*1200,geom1,a,8,0,flag)
 CALL WRITE (scr1,a,7,0)
 120 IF (nosclr == 0) GO TO 140
 CALL READ  (*1180,*1200,scr2,a(10),1,0,flag)
 CALL WRITE (scr1,a(10),7,0)
 130 IF (nogrid == 0) GO TO 160
 IF (nosclr == 0) GO TO 140
 IF (a(1) -  a(10) < 0.0) THEN
   GO TO   140
 ELSE IF (a(1) -  a(10) == 0.0) THEN
   GO TO  1250
 ELSE
   GO TO   160
 END IF
 
!     GRID NO. .LT. SCALAR NO.
 
 140 i = i + 2
 z(i) = a(1)
 
!     GRID POINT EXTERNAL ID * MULT IS LIMITED TO COMPUTER MAXIMUM
!     INTEGER SIZE
 
 IF (a(1) <= imax .OR. axi /= 0) GO TO 142
 IF (a(1) > maxa1) maxa1 = a(1)
 GO TO 146
 142 z(i+1) = mult*a(1)
 146 CALL READ  (*1180,*150,geom1,a,8,0,flag)
 CALL WRITE (scr1,a,7,0)
 GO TO 130
 150 nogrid = 0
 a(1) = large
 IF (nosclr == 0) GO TO 180
 
!     SCALAR NO. .LT. GRID NO.
 
 160 i = i + 2
 z(i) = a(10)
 
!     SCALAR POINT EXTERNAL ID * MULT IS LIMITED TO COMPUTER MAXIMUM
!     INTEGER SIZE
 
 IF (a(10) <= imax .OR. axi /= 0) GO TO 162
 IF (a(10) > maxa1) maxa1 = a(10)
 GO TO 166
 162 z(i+1) = mult*a(10)
 166 CALL READ  (*1180,*170,scr2,a(10),1,0,flag)
 CALL WRITE (scr1,a(10),7,0)
 GO TO 130
 170 nosclr = 0
 a(10)  = large
 IF (nogrid == 0) GO TO 180
 GO TO 140
 
!     LIST COMPLETE ONLY IF MAXA1 .LE. ZERO
 
!     IF MAXA1 IS .GT. ZERO, SOME LARGE GRID OR SCALAR POINTS HAD BEEN
!     LEFT OUT IN LIST. MAXA1 IS THE LARGEST GRID OR SCALAR POINT
!     EXTERNAL ID.  RESET MULT AND REPEAT COMPILING LIST
 
 180 IF (maxa1 <= 0) GO TO 185
 IF (isubs /= 0) GO TO 183
 CALL REWIND (scr1)
 CALL REWIND (geom1)
 IF (nosclr /= 0) CALL REWIND (scr2)
 mult = 100
 IF (maxa1 > imax*10) mult = 10
 imax  = (imax/mult)*1000
 maxa1 = -1
!WKBR CALL PAGE (-3)
 CALL page2(-3)
 IF (mult == 100) WRITE (iout,182) uwm,lvl1
 IF (mult ==  10) WRITE (iout,182) uwm,lvl2
 182 FORMAT (a25,' 2140A, DUE TO THE PRESENCE OF ONE OR MORE GRID OR ',  &
     'SCALAR POINTS WITH VERY LARGE EXTERNAL ID''S, THE SEQGP' ,  &
 /5X,'AND SEQEP CARDS, IF USED, ARE FORCED TO REDUCE FROM ',  &
     'ALLOWABLE 4 SEQID LEVELS TO ',a29,/)
 GO TO 110
 
 183 WRITE  (iout,184) ufm
 184 FORMAT (a23,' 2140B, EXTERNAL GRID OR SCALAR POINT ID TOO BIG')
 CALL mesage (-61,0,0)
 
 185 n    = i
 neqex= n
 n1   = n + 1
 n2   = n + 2
 igpdt= n2
 ilist= n2
 kn   = n1/2
 CALL CLOSE (scr1,clsrew)
 CALL CLOSE (scr2,clsrew)
 GO TO 210
 
!     NO GRID CARDS PRESENT-- TEST FOR ANY SCALAR PTS.
 
 190 nogmp1 = 0
 200 IF (nosclr == 0) GO TO 980
 GO TO 120
 
!     READ THE SEQGP TABLE (IF PRESENT)
!     FOR EACH ENTRY, FIND MATCH IN THE SORTED EXTERNAL GRID POINTS
!     AND REPLACE SEQUENCE NO. WITH SEQGP NO.
 
 210 noseq  = 0
 nogpdt = 1
 IF (nogmp1 == 0) GO TO 260
 ASSIGN 230 TO ndx
 spoint(2) = 0
 ierr = 1
 ASSIGN 220 TO nerr
 CALL locate (*250,z(buf1),seqgp,flag)
 noseq = 1
 ifail = 0
 2010 CALL READ (*1180,*2020,geomp,z(n2),buf1-1,1,flag)
 ifail = ifail + 1
 GO TO 2010
 2020 IF (ifail == 0) GO TO 2060
 nwds = (ifail-1)*(buf1-1) + flag
 WRITE  (iout,2040) ufm,nwds
 2040 FORMAT (a23,' 3135, UNABLE TO PROCESS SEQGP DATA IN SUBROUTINE ',  &
     'GP1 DUE TO INSUFFICIENT CORE.', //5X,  &
     'ADDITIONAL CORE REQUIRED =',i10,7H  words)
 CALL mesage (-61,0,0)
 
!     CHECK FOR MULTIPLE REFERENCES TO GRID (OR SCALAR) POINT ID NOS.
!     AND SEQUENCE ID NOS. ON SEQGP CARDS
 
 2060 k  = n2
 kk = n2 + flag - 1
 jj = kk - 2
 2080 DO  i = k,jj,2
   IF (z(i) < 0 .OR. i >= kk) GO TO 2275
   ii = i + 2
   ifail = 0
   DO  j = ii,kk,2
     IF (z(i) /= z(j)) CYCLE
     IF (ifail /=   0) GO TO 2260
     ifail = 1
     nogo  = 1
     IF (k /= n2) GO TO 2110
     WRITE  (iout,2100) ufm,z(i)
     2100 FORMAT (a23,' 3136, MULTIPLE REFERENCES TO GRID (OR SCALAR) POINT'  &
         ,      ' ID NO.',i9,'  ON SEQGP CARDS.')
     GO TO 2260
     2110 idseq1 = z(i)/1000
     irmndr = z(i) - 1000*idseq1
     IF (irmndr /= 0 .AND. mult >= 10) GO TO 2140
     IF (axi /= 0) GO TO 2130
     WRITE  (iout,2120) ufm,idseq1
     2120 FORMAT (a23,' 3137, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,6X,  &
         ' ON SEQGP CARDS.')
     GO TO 2260
     2130 IF (axi == 1) WRITE (iout,2135) ufm
     2135 FORMAT (a23,' 3137A, SEQGP CARDS WITH MORE THAN ONE SEQID LEVEL ',  &
         'ARE ILLEGAL FOR AXISYSM. OR HYDROELAS. PROBLEM')
     axi  = 2
     nogo = 1
     GO TO 2260
     2140 idseq2 = irmndr/100
     irmndr = irmndr - 100*idseq2
     IF (irmndr /= 0 .AND. mult >= 100) GO TO 2180
     WRITE  (iout,2160) ufm,idseq1,idseq2
     2160 FORMAT (a23,' 3137, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,5X,'ON SEQGP CARDS.')
     GO TO 2260
     2180 idseq3 = irmndr/10
     irmndr = irmndr - 10*idseq3
     IF (irmndr /= 0) GO TO 2220
     WRITE  (iout,2200) ufm,idseq1,idseq2,idseq3
     2200 FORMAT (a23,' 3137, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,1H.,i1,4X,'ON SEQGP CARDS.')
     GO TO 2260
     2220 WRITE  (iout,2240) ufm,idseq1,idseq2,idseq3,irmndr
     2240 FORMAT (a23,' 3137, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,1H.,i1,1H.,i1,'  ON SEQGP CARDS.')
     2260 z(j) = -z(j)
   END DO
   
   2275 IF (jj < kk .OR. mult == 1000) CYCLE
   l = z(i)
   IF (mult   <=   10) GO TO 2280
   IF (MOD(l,10) /= 0) GO TO 2276
   z(i) = l/10
   CYCLE
   2276 IF (maxa1 == 0) CYCLE
   maxa1 = 0
   nogo  = 1
   WRITE  (iout,2277) ufm
   2277 FORMAT (a23,' 2140B, ILLEGAL DATA IN SEQGP CARD, POSSIBLY CAUSED',  &
       ' BY LARGE GRID OR SCALAR POINTS')
   CYCLE
   2280 IF (mult == 1) GO TO 2282
   IF (MOD(l,100) /= 0) GO TO 2276
   z(i) = l/100
   CYCLE
   2282 IF (axi == 0) CALL mesage (-37,0,nam)
   IF (MOD(l,1000) == 0) CYCLE
   IF (axi == 1) WRITE (iout,2135) ufm
   axi  = 2
   nogo = 1
 END DO
 
 IF (k /= n2) GO TO 2290
 jj = kk
 k  = k + 1
 GO TO 2080
 
 2290 DO  i = n2,kk,2
   IF (z(i) < 0) z(i) = -z(i)
 END DO
 IF (nogo == 1) GO TO 2400
 
!     CHECK TO SEE IF ANY SEQUENCE ID NO. ON SEQGP CARDS IS THE SAME
!     AS A GRID (OR SCALAR) POINT ID NO. THAT HAS NOT BEEN RESEQUENCED
 
 loop2390:  DO  i = k,kk,2
   IF (z(i) < 0) CYCLE loop2390
   idseq1 = z(i)/mult
   irmndr = z(i) - mult*idseq1
   IF (irmndr /= 0) CYCLE loop2390
   DO  j = n2,kk,2
     IF (idseq1 == z(j)) CYCLE loop2390
   END DO
   DO  j = 1,n1,2
     IF (idseq1 == z(j)) GO TO 2360
   END DO
   CYCLE loop2390
   2360 nogo = 1
   WRITE  (iout,2380) ufm,idseq1
   2380 FORMAT (a23,' 3138, SEQUENCE ID NO.',i6,' ON SEQGP CARDS IS THE ',  &
       'SAME AS A', /5X,'GRID (OR SCALAR) POINT ID NO. THAT HAS ',  &
       'NOT BEEN RESEQUENCED.')
 END DO loop2390
 2400 CONTINUE
 i = -1
 220 i = i + 2
 IF (i > flag) GO TO 240
 a(1) = z(n2+i-1)
 a(2) = z(n2+i  )
 GO TO 1060
 230 z(2*k) = a(2)
 GO TO 220
 
!     SORT THE CORE TABLE BY INTERNAL GRID PT NO
!     THUS FORMING THE GPL (EXTERNAL GRID PT NOS IN SORT BY INTERNAL NO)
 
 240 IF (nogo /= 0) GO TO 1165
 CALL sorti (0,0,2,2,z,n1)
 
!     CLOSE GEOM1. WRITE THE GPL. FIRST RECORD IS A SINGE ENTRIED LIST
!     OF EXTERNAL GRID NOS. IN INTERNAL SORT. SECOND RECORD IS A DOUBLE
!     ENTRIED LIST OF EXTERAL GRID NO., SEQUENCE NO. (SORT IS INTERNAL).
!     ADD THE MULTIPLIER, MULT, TO THE 3RD WORD OF GPL HEADER RECORD
 
 250 IF (nogmp1 /= 0) CALL CLOSE (geom1,clsrew)
 260 CALL fname (gpl,a)
 FILE = gpl
 CALL OPEN (*1170,gpl,z(buf1),wrtrew)
 a(3) = mult
 CALL WRITE (gpl,a,3,1)
 DO  i = 1,n,2
   CALL WRITE (gpl,z(i),1,0)
 END DO
 CALL WRITE (gpl,0,0,1)
 CALL WRITE (gpl,z,n1,1)
 CALL CLOSE (gpl,clsrew)
 mcb(1) = gpl
 CALL wrttrl (mcb)
 
!     FORM INTERNAL INDEX FOR EACH EXTERNAL GRID PT. NO.
 
 i = 2
 z(i) = 1
 IF (n == 1) GO TO 310
 DO  i = 3,n,2
   z(i+1) = z(i-1) + 1
 END DO
 
!     TEST TO SEE IF EXTERNAL GRID PT NOS ARE STILL IN EXTERNAL SORT
!     I.E., IF NO SEQGP TABLE, THEN SORT IS MAINTAINED
!     OTHERWISE, SORT ON EXTERNAL GRID NO.
 
 IF (noseq /= 0) CALL sorti (0,0,2,1,z,n1)
 
!     DETERMINE IF THE GPDT CAN BE HELD IN CORE
!     NWDS= TOTAL NO OF WORDS IN THE GPDT
!     M= MAX NO OF ENTRIES CORE CAN HOLD WITH ONE BUFFER OPEN
!     IF NWDS/7.LE.M,CORE WILL HOLD THE GPDT
!     OTHERWISE THE FILE SORT ROUTINE WILL BE USED
 
 310 nwds = 7*kn
 m    = (buf1-n1)/7
 gpfl = 0
 IF (kn > m) gpfl = 7
 FILE = scr1
 
!     READ THE GRID AND SPOINT TABLES FROM SCR1
!     REPLACE THE EXTERNAL GRID PT NO WITH THE INTERNAL INDEX
!     IF CORE WILL HOLD THE GPDT, USE THE INTERNAL INDEX AS A POINTER
!     OTHERWISE, WRITE THE UNSORTED GPDT ON SCR2
 
 CALL OPEN (*1170,scr1,z(buf1),rdrew)
 FILE = scr2
 IF (gpfl /= 0) CALL OPEN (*1170,scr2,z(buf2),wrtrew)
 FILE = scr1
 ASSIGN 340 TO ndx
 ierr = 2
 ASSIGN 330 TO nerr
 330 CALL READ (*1180,*370,scr1,a,7,0,flag)
 GO TO 1060
 340 IF (gpfl /= 0) GO TO 360
 j = n1 + 7*(a(1)-1)
 DO  k = 1,7
   i = j+k
   z(i) = a(k)
 END DO
 GO TO 330
 360 CALL WRITE (scr2,a,7,0)
 GO TO 330
 370 IF (nogo /= 0) GO TO 1165
 CALL CLOSE (scr1,clsrew)
 
!     OPEN OUTPUT FILE FOR GPDT AND WRITE HEADER DATA
!     IF GPDT IS IN CORE, WRITE IT OUT
 
 FILE = gpdt
 CALL fname (gpdt,a)
 CALL OPEN  (*1170,gpdt,z(buf1),wrtrew)
 CALL WRITE (gpdt,a,2,1)
 IF (gpfl /= 0) GO TO 390
 CALL WRITE (gpdt,z(igpdt),nwds,1)
 GO TO 400
 
!     IF GPDT NOT IN CORE, CALL SORT
 
 390 nfile(1) = scr1
 nfile(2) = cstm
 nfile(3) = bgpdt
 CALL CLOSE (scr2,clsrew)
 FILE = scr2
 CALL OPEN  (*1170,scr2,z(buf2),rdrew)
 CALL sorti (scr2,gpdt,7,1,z(igpdt),buf2-igpdt)
 CALL CLOSE (scr2,clsrew)
 400 CALL CLOSE (gpdt,clsrew)
 mcb(1) = gpdt
 CALL wrttrl (mcb)
 
!     READ THE CORDIJ TABLES INTO CORE (IF PRESENT)
 
 ifl = -1
 m   = icsdt
 nolist = 0
 IF (nogmp1 == 0) GO TO 810
 ndx   = buf1 - 15
 ncore = buf1 - 15
 DO  i = icsdt,buf1
   z(i) = 0
 END DO
 FILE = geomp
 CALL preloc (*1170,z(buf1),geomp)
 DO  i = 1,6
   ij = i + i - 1
   CALL locate (*440,z(buf1),cordij(ij),flag)
   ifl = 1
   430 CALL READ (*1180,*440,geomp,z(m),cord(i),0,flag)
   m = m + 16
   IF (m > ncore) CALL mesage (-8,0,gp1ah)
   GO TO 430
   440 CONTINUE
 END DO
 CALL CLOSE (geomp,clsrew)
 m = m - 16
 ncsdt = m
 
!     TEST FOR PRESENCE OF ANY CORDIJ TABLES
 
 IF (ifl == -1) GO TO 810
 
!     REPLACE EXTERNAL GRID PT NO IN CORD1J ENTRIES (IF ANY)
!     WITH CORRESPONDING INTERNAL INDEX
!     SAVE A TABLE OF GRID PTS REFERENCED ON CORD1J ENTRIES
 
 jj   = icsdt
 ilist= ncsdt + 16
 ii   = ilist - 1
 ncore= buf1  - 3
 ierr = 3
 470 IF (z(jj+2) /= 1) GO TO 510
 nolist = 1
 ASSIGN 480 TO ndx
 ASSIGN 485 TO nerr
 a(1) = z(jj+3)
 spoint(2) =  z(jj+1)
 GO TO 1060
 480 z(jj+3) = a(1)
 z(ii+1) = a(1)
 485 ASSIGN 490 TO ndx
 ASSIGN 495 TO nerr
 a(1) = z(jj+4)
 GO TO 1060
 490 z(jj+4) = a(1)
 z(ii+2) = a(1)
 495 ASSIGN 500 TO ndx
 ASSIGN 505 TO nerr
 a(1) = z(jj+5)
 GO TO 1060
 500 z(jj+5) = a(1)
 z(ii+3) = a(1)
 505 ii = ii+3
 IF (ii > ncore) CALL mesage (-8,0,gp1ah)
 510 jj = jj + 16
 IF (jj <= ncsdt) GO TO 470
 IF (nogo /=   0) GO TO 1165
 
!     IF ANY CORD1J ENTRIES, PASS THE GPDT AND CREATE A TABLE OF THE
!     REFERENCED GRID PTS. THIS TABLE IS CALLED CSGP
 
 IF (nolist == 0) GO TO 550
 nlist = ii
 icsgp = nlist + 1
 CALL sorti (0,0,1,1,z(ilist),icsgp-ilist)
 z(icsgp) = 0
 jj = ilist
 DO  kk = ilist,nlist
   IF (z(kk+1) == z(kk)) CYCLE
   z(jj) = z(kk)
   jj = jj + 1
 END DO
 nlist = jj - 1
 icsgp = jj
 FILE  = gpdt
 CALL OPEN   (*1170,gpdt,z(buf1),rdrew)
 CALL fwdrec (*1180,gpdt)
 ncore = buf1 - 5
 i = ilist
 540 CALL READ (*1180,*1200,gpdt,z(jj),7,0,flag)
 IF (z(jj) /= z(i)) GO TO 540
 jj = jj + 5
 IF (jj > ncore) CALL mesage (-8,0,gp1ah)
 i  = i + 1
 IF (i <= nlist) GO TO 540
 ncsgp = jj - 5
 CALL CLOSE (gpdt,clsrew)
 
!     LOOP THRU THE CSDT SOLVING AS MANY COORDINATE SYSTEMS AS POSSIBLE
!     ON EACH PASS.
 
 550 nn = (ncsdt-icsdt)/16 + 1
 solv  = 0
 solvp = 0
 560 ii = icsdt
 570 IF (z(ii+2)-2 < 0.0) THEN
   GO TO   580
 ELSE IF (z(ii+2)-2 == 0.0) THEN
   GO TO   620
 ELSE
   GO TO   690
 END IF
 
!     *****  TYPE = 1 *****
!     CHECK TO SEE IF EACH OF THE 3 REFERENCE GRID PTS IS IN BASIC SYS
!     IF SO,CALCULATE THE TRANSFORMATION TO BASIC AND SET COORD SYSTEM
!     AS SOLVED, IF NOT CONTINUE TO NEXT COORDINATE SYSTEM
 
 580 i = 0
 590 k = ii + i
 j = icsgp - 1
 600 IF (z(j+1) == z(k+3)) GO TO 610
 j = j + 5
 IF (j < ncsgp) GO TO 600
 GO TO 1220
 610 IF (z(j+2) /= 0) GO TO 700
 k = i*3
 aa(k+1) = zz(j+3)
 aa(k+2) = zz(j+4)
 aa(k+3) = zz(j+5)
 i = i+1
 IF (i <= 2) GO TO 590
 GO TO 1020
 
!     ***** TYPE = 2 *****
!     CHECK THE DEFINING LOCAL COORDINATE SYSTEM
!     IF BASIC, SOLVE AS IN TYPE=1
!     IF NOT BASIC, FIND THE REFERENCED COORD SYSTEM AND TEST IF THAT
!     SYSTEM IS SOLVED. IF YES, CALCULATE THE TRANSFORMATION TO BASIC
!     IF NO, CONTINUE THRU THE CSDT
 
 620 IF (z(ii+3) /= 0) GO TO 640
 DO  i = 1,9
   k = ii + i
   aa(i) = zz(k+3)
 END DO
 GO TO 1020
 640 i = icsdt
 650 IF (z(i) == z(ii+3)) GO TO 660
 i = i + 16
 IF (i <= ncsdt) GO TO 650
 GO TO 1230
 660 IF (z(i+2) /= 3 .OR. z(i+3) /= 0) GO TO 700
 k = 0
 ASSIGN 680 TO ndx
 670 l = k + ii
 ax(1) = zz(l+4)
 ax(2) = zz(l+5)
 ax(3) = zz(l+6)
 IF (z(i+1)-2 < 0.0) THEN
   GO TO   990
 ELSE IF (z(i+1)-2 == 0.0) THEN
   GO TO  1000
 ELSE
   GO TO  1010
 END IF
 680 aa(k+1) = ar(1)
 aa(k+2) = ar(2)
 aa(k+3) = ar(3)
 k = k + 3
 IF (k <= 6) GO TO 670
 GO TO 1020
 
!     ***** TYPE = 3 *****
!     CHECK THE DEFINING LOCAL COORDINATE SYSTEM
!     IF BASIC, CONTINUE THRU CSDT
!     IF NOT BASIC, ERROR CONDITION
 
 690 IF (z(ii+3) /= 0) GO TO 1190
 
!     TEST FOR COMPLETION OF PASS THRU CSDT
 
 700 ii = ii + 16
 IF (ii <= ncsdt) GO TO 570
 
!     LOOP THRU THE CSGP (IFPRESENT) AND TRANSFORM ALL
!     POSSIBLE GRID PTS TO BASIC
 
 IF (nolist == 0) GO TO 770
 jj = icsgp
 720 IF (z(jj+1) == 0) GO TO 760
 i = icsdt
 730 IF (z(jj+1) == z(i)) GO TO 740
 i = i + 16
 IF (i <= ncsdt) GO TO 730
 ierr = 6
 spoint(1) = z(jj  )
 spoint(2) = z(jj+1)
 GO TO 1190
 740 IF (z(i+2) /= 3 .OR. z(i+3) /= 0) GO TO 760
 ax(1) = zz(jj+2)
 ax(2) = zz(jj+3)
 ax(3) = zz(jj+4)
 ASSIGN 750 TO ndx
 IF (z(i+1)-2 < 0.0) THEN
   GO TO   990
 ELSE IF (z(i+1)-2 == 0.0) THEN
   GO TO  1000
 ELSE
   GO TO  1010
 END IF
 750 zz(jj+2) = ar(1)
 zz(jj+3) = ar(2)
 zz(jj+4) = ar(3)
 zz(jj+1) = 0
 760 jj = jj + 5
 IF (jj <= ncsgp) GO TO 720
 
!     TEST TO SEE IF ALL COORDINATE SYSTEMS SOLVED
!     IF NOT, TEST TO SEE IF ANY NEW SOLUTIONS ON LAST PASS
!     IF NONE, INCONSISTANT DEFINITION OF COORDINATE SYSTEMS
!     OTHERWISE LOOP BACK THRU THE CSDT
 
 770 IF (solv ==    nn) GO TO 780
 IF (solv == solvp) GO TO 1240
 solvp = solv
 GO TO 560
 
!     WRITE THE CSTM
 
 780 CALL fname (cstm,a)
 FILE = cstm
 CALL OPEN  (*1170,cstm,z(buf1),wrtrew)
 CALL WRITE (cstm,a,2,1)
 DO  ii = icsdt,ncsdt,16
   CALL WRITE (cstm,z(ii),2,0)
   CALL WRITE (cstm,z(ii+4),12,0)
 END DO
 CALL CLOSE (cstm,clsrew)
 nocstm = nn
 mcb(3) = nn
 mcb(1) = cstm
 CALL wrttrl (mcb)
 
!     OPEN EQEXIN AND WRITE HEADER RECORD.
!     THEN WRITE FIRST RECORD (PAIRS OF EXTERNAL GRID NO., INTERNAL NO.
!     IN EXTERNAL SORT).
 
 810 FILE = eqexin
 CALL OPEN  (*1170,eqexin,z(buf1),wrtrew)
 CALL fname (eqexin,a)
 CALL WRITE (eqexin,a,2,1)
 CALL WRITE (eqexin,z,n1,1)
 CALL CLOSE (eqexin,cls)
 
!     A LIST OF DEGREES OF FREEDOM FOR EACH GRID OR SCALAR POINT IS
!     FORMED BEGINNING AT Z(ILIST) BY READING GEOM2 AND USING THE
!     CONNECTION INFORMATION IN CONJUNCTION WITH THE ELEM TABLE IN
!     /GPTA1/.
 
 FILE   = geom2
 ilist0 = ilist - 1
 nlist  = ilist + (neqex+1)/2
 IF (nlist >= buf3) CALL mesage (-8,0,gp1ah)
 DO  i = ilist,nlist
   z(i) = 0
 END DO
 jerr = 0
 CALL OPEN   (*8130,geom2,z(buf1),rdrew)
 8103 CALL fwdrec (*1180,geom2)
 8104 CALL ectloc (*8130,geom2,a,i)
 
!     ELEMENT TYPE LOCATED--PREPARE TO PROCESS EACH ELEMENT
 
 IF (elem(i+9) == 0) GO TO 8103
 j1 = elem(i+12)
 nread = j1 + elem(i+9) - 1
 nskip =-(elem(i+5 ) - nread)
 maxdof=  elem(i+24)
 itype =  elem(i+2 )
 
!     READ CONNECTION DATA FOR ELEMENT AND LOCATE EXT. GRID NBR IN
!     EQEXIN UPDATE DOF LIST FOR EACH GRID NBR
 
 8110 CALL READ (*1180,*8104,geom2,a,nread,0,m)
 DO  i = j1,nread
   IF (a(i) == 0) CYCLE
   CALL bisloc (*8122,a(i),z,2,kn,k)
   j = ilist0 + z(k+1)
   IF (itype >= 76 .AND. itype <= 79) GO TO 8115
   
!     STRUCTURE ELEMENT AND OTHERS
   
   IF (z(j) < 0) GO TO 8124
   z(j) = MAX0(z(j),maxdof)
   CYCLE
   
!     FLUID ELEMENT (CFHEX1,CFHEX2,CFWEDGE,CFTETRA)
   
   8115 IF (z(j) > 0) GO TO 8124
   
   z(j) = -1
   CYCLE
   8122 WRITE  (iout,8123) ufm,a(1),a(i)
   8123 FORMAT (a23,' 2007, ELEMENT',i8,' REFERENCES UNDEFINED GRID ',  &
       'POINT',i8)
   jerr = jerr + 1
   CYCLE
   8124 WRITE  (iout,8125) ufm,a(i)
   8125 FORMAT (a23,' 8011, GRID POINT',i8,' HAS BOTH STRUCTURE AND ',  &
       'FLUID ELEMENTS CONNECTED')
   jerr = jerr + 1
 END DO
 CALL READ (*1180,*8104,geom2,a,nskip,0,m)
 GO TO 8110
 
!     END-OF-FILE ON GEOM2---IF FATAL ERRORS, TERMINATE
 
 8130 CONTINUE
 IF (jerr /= 0) CALL mesage (-61,a,z)
 
!     OPEN BGPDT AND SIL. WRITE HEADER RECORDS. OPEN GPDT. SKIP HEADER.
 
 offset = rshift(kn,5)
 CALL fname (bgpdt,a)
 CALL fname (sil,a(3))
 FILE = bgpdt
 CALL OPEN (*1170,bgpdt,z(buf1),wrtrew)
 FILE = sil
 CALL OPEN (*1170,sil,z(buf2),wrtrew)
 FILE = gpdt
 CALL OPEN   (*1170,gpdt,z(buf3),rdrew)
 CALL fwdrec (*1180,gpdt)
 CALL WRITE  (bgpdt,a,2,1)
 CALL WRITE  (sil,a(3),2,1)
 luset = 1
 
!     READ AN ENTRY FROM THE GPDT.
!     TEST FOR DEFINING COORDINATE SYSTEM.
 
 820 CALL READ (*1180,*970,gpdt,a,7,0,flag)
 IF (a(2) < 0.0) THEN
   GO TO   910
 ELSE IF (a(2) == 0.0) THEN
   GO TO   880
 END IF
 
!     COORDINATE SYSTEM NOT BASIC--
!     USE CSDT IN CORE TO TRANSFORM TO BASIC.
 
 830 IF (nocstm == -1) GO TO 850
 i = icsdt
 840 IF (z(i) == a(2)) GO TO 860
 i = i + 16
 IF (i <= ncsdt) GO TO 840
 850 ierr = 6
 spoint(1) = a(1)
 spoint(2) = a(2)
 GO TO 1190
 860 ax(1) = aa(3)
 ax(2) = aa(4)
 ax(3) = aa(5)
 ASSIGN 870 TO ndx
 IF (z(i+1)-2 < 0.0) THEN
   GO TO   990
 ELSE IF (z(i+1)-2 == 0.0) THEN
   GO TO  1000
 ELSE
   GO TO  1010
 END IF
 870 aa(3) = ar(1)
 aa(4) = ar(2)
 aa(5) = ar(3)
 
!     GRID POINT NOW BASIC--
!     STORE DISPLACEMENT SYSTEM COORD. SYSTEM ID AND SET TYPE.
!     MAKE SURE DISPLACEMENT COORD. SYSTEM IS DEFINED.
 
 880 a(2) = a(6)
 TYPE = 1
 khr  = ilist0 + a(1)
 incr = z(khr)
 
!     IF INCR NEGATIVE - SPECIAL HYDROELASTIC GRID POINT WITH SINGLE
!     DEGREE OF FREEDOM
 
 IF (incr < 0) GO TO 905
 IF (incr == 0) incr = 6
 
!     ///////////////////////////////
 
!     TEMP PATCH
 
 incr = MAX0(incr,6)
 
!     ///////////////////////////////
 
 IF (a(2) == 0 .AND. itherm == 0) GO TO 920
 
!     IF A(2) WHICH EQUALS A(6) IS EQUAL TO -1 THEN A FLUID GRID POINT
!     AS CREATED BY IFP4 IS AT HAND AND HAS ONLY 1 DEGREE OF FREEDOM
!     ..... IF -HEAT- PROBLEM THEN ALL GRIDS HAVE 1 DEGREE OF FREEDOM.
 
 IF (a(2) == (-1) .OR. itherm > 0) GO TO 905
 IF (nocstm == -1) GO TO 900
 DO  ijk = icsdt,ncsdt,16
   IF (a(2) == z(ijk)) GO TO 920
 END DO
 900 nogo = 1
 CALL mesage (30,104,a(2))
 GO TO 920
 
!     SCALAR POINT-- SET TYPE.
 
 905 a(2) = 0
 a(6) = 0
 910 TYPE = 2
 incr = 1
 
!     WRITE ENTRY ON BGPDT AND SIL.
 
 920 CALL WRITE (bgpdt,a(2),4,0)
 CALL WRITE (sil, luset,1,0)
 
!     REPLACE INTERNAL NO. IN EQEXIN WITH CODED SIL NO.
!     THEN INCREMENT SIL NO.
 
 ncode = 10*luset + TYPE
 IF (noseq /= 0) GO TO 925
 k = 2*a(1)
 IF (z(k) - a(1) == 0.0) THEN
   GO TO   960
 ELSE
   GO TO   950
 END IF
 925 ncode = -ncode
 lmt1  = MAX0(2*(a(1)-offset),2)
 DO  k = lmt1,n1,2
   IF (z(k) == a(1)) GO TO 960
 END DO
 DO  k = 2,lmt1,2
   IF (z(k) == a(1)) GO TO 960
 END DO
 950 CALL mesage (-30,2,a)
 960 z(k)  = ncode
 luset = luset + incr
 GO TO 820
 
!     CLOSE BGPDT AND SIL. WRITE TRAILERS.
 
 970 CALL CLOSE (bgpdt,clsrew)
 CALL CLOSE (sil,clsrew)
 CALL CLOSE (gpdt,clsrew)
 luset = luset - 1
!    2147483647   = 2**31-1
 IF (luset <= 2147483647) GO TO 974
 WRITE (iout,972) ufm,luset
 972 FORMAT (a23,' 3175, TOTAL NUMBER OF DEGREES OF FREEDOM IN THE ',  &
     'PROBLEM (',i11,' ) EXCEEDS 2,147,483,647 (I.E., ' '2**31 - 1)')
 973 FORMAT (a29,' 3175, PROBLEM SIZE,',i8,' DOF''S, EXCEEDS THE OLD ',  &
     'LIMIT OF 65535.', /5X,'GOOD NEWS, JOB WILL CONTINUE')
 CALL mesage (-61,0,0)
 974 mcb(1) = bgpdt
 mcb(3) = 0
 CALL wrttrl (mcb)
 mcb(1) = sil
 mcb(3) = luset
 CALL wrttrl (mcb)
 
!     IF GRID NOS. ARE RESEQUENCED, SWITCH SIGN ON CODED SIL NO.
!     WRITE SECOND RECORD OF EQEXIN. CLOSE FILE AND WRITE TRAILER.
 
 IF (noseq == 0) GO TO 978
 DO  k = 2,n1,2
   z(k) = -z(k)
 END DO
 978 FILE = eqexin
 CALL OPEN  (*1170,eqexin,z(buf1),wrt)
 CALL WRITE (eqexin,z,n1,1)
 CALL CLOSE (eqexin,clsrew)
 mcb(1) = eqexin
 mcb(3) = 0
 CALL wrttrl (mcb)
 CALL sswtch (36,k)
 IF (k    == 1) CALL diag36 (z,buf1,gpl,sil,eqexin)
 IF (nogo /= 0) CALL mesage (-61,0,0)
 RETURN
 
!     ABNORMAL EXIT FROM GP1
 
 980 CALL CLOSE (scr1,clsrew)
 CALL CLOSE (geom1,clsrew)
 nocstm = -1
 RETURN
 
!     ===============================================================
 
!     INTERNAL SUBROUTINE TO TRANSFORM A RECTANGULAR GRID PT TO BASIC
!     I POINTS TO THE CSDT ENTRY WHERE THE TRANSFORMATION IS DEFINED
!     THE GRID PT TO BE TRANSFORMED IS STORED AT AX(1,2,3)
!     THE TRANSFORMED GRID PT WILL BE STORED AT AR(1,2,3)
 
 990 ar(1) = zz(i+ 7)*ax(1) + zz(i+ 8)*ax(2) + zz(i+ 9)*ax(3) + zz(i+4)
 ar(2) = zz(i+10)*ax(1) + zz(i+11)*ax(2) + zz(i+12)*ax(3) + zz(i+5)
 ar(3) = zz(i+13)*ax(1) + zz(i+14)*ax(2) + zz(i+15)*ax(3) + zz(i+6)
 GO TO ndx, (680,750,870)
 
!     INTERNAL SUBROUTINE TO TRANSFORM A CYLINDRICAL GRID PT TO BASIC
!     R,THETA,Z IS STORED AX(1,2,3)
 
 1000 r     = ax(1)
 ax(2) = degra*ax(2)
 ax(1) = r*COS(ax(2))
 ax(2) = r*SIN(ax(2))
 GO TO 990
 
 
!     INTERNAL SUBROUTINE TO TRANSFORM A SPHERICAL GRID PT TO BASIC
!     RHO,THETA,PHI IS STORED AT AX(1,2,3)
 
 1010 ax(2) = degra*ax(2)
 ax(3) = degra*ax(3)
 rsth  = ax(1)*SIN(ax(2))
 rcth  = ax(1)*COS(ax(2))
 ax(1) = rsth *COS(ax(3))
 ax(2) = rsth *SIN(ax(3))
 ax(3) = rcth
 GO TO 990
 
 
!     INTERNAL SUBROUTINE TO CALCULATE THE 3X3 TRANSFORMATION MATRIX
!     AND 3X1 TRANSLATION VECTOR GIVEN THREE POINTS IN THE BASIC SYSTEM
!     THE RESULTS ARE STORED BACK IN THE CSDT
 
!     STORE R0 = A IN THE CSDT
 
 1020 zz(ii+4) = aa(1)
 zz(ii+5) = aa(2)
 zz(ii+6) = aa(3)
 
!     FORM B - A
 
 DO  i = 1,3
   ak(i) = ab(i) - aa(i)
 END DO
 
!     FORM K = (B - A)/LENGTH(B - A)
!     FORM C - A
 
 length = SQRT(ak(1)**2 + ak(2)**2 + ak(3)**2)
 DO  i = 1,3
   ak(i) = ak(i)/length
   ac(i) = ac(i) - aa(i)
 END DO
 
!     FORM K X (C - A)
 
 aj(1) = ak(2)*ac(3) - ak(3)*ac(2)
 aj(2) = ak(3)*ac(1) - ak(1)*ac(3)
 aj(3) = ak(1)*ac(2) - ak(2)*ac(1)
 
!     FORM J = (K X (C-A))/LENGTH(K X (C-A))
 
 length =  SQRT(aj(1)**2 + aj(2)**2 + aj(3)**2)
 DO  i = 1,3
   aj(i) = aj(i)/length
 END DO
 
!     FORM I = J X K
 
 ai(1) = aj(2)*ak(3) - aj(3)*ak(2)
 ai(2) = aj(3)*ak(1) - aj(1)*ak(3)
 ai(3) = aj(1)*ak(2) - aj(2)*ak(1)
 
!     STORE 3X3 ROTATION MATRIX = ((IX,JX,KX),(IY,JY,KY),(IZ,JZ,KZ))
!     IN THE CSDT
 
 zz(ii+ 7) = ai(1)
 zz(ii+ 8) = aj(1)
 zz(ii+ 9) = ak(1)
 zz(ii+10) = ai(2)
 zz(ii+11) = aj(2)
 zz(ii+12) = ak(2)
 zz(ii+13) = ai(3)
 zz(ii+14) = aj(3)
 zz(ii+15) = ak(3)
 
!     SET WD 3 OF CSDT = 3 AND WD 4 = 0 TO INDICATE  SOLVED SYSTEM
!     INCREMENT SOLVED SYSTEM COUNT
 
 z(ii+2) = 3
 z(ii+3) = 0
 solv = solv + 1
 GO TO 700
 
 
!     INTERNAL SUBROUTINE TO PERFORM BINARY SEARCH ON FIRST ENTRY
!     OF A DOUBLE ENTRIED TABLE STORED AT Z(1) THRU Z(N+1)
 
 1060 klo = 1
 khi = kn
 1070 k = (klo+khi+1)/2
 1080 IF (a(1)-z(2*k-1) < 0.0) THEN
   GO TO  1090
 ELSE IF (a(1)-z(2*k-1) == 0.0) THEN
   GO TO  1150
 ELSE
   GO TO  1100
 END IF
 1090 khi = k
 GO TO 1110
 1100 klo = k
 1110 IF (khi-klo-1 < 0) THEN
   GO TO  1160
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO  1120
 ELSE
   GO TO  1070
 END IF
 1120 IF (k == klo) GO TO 1130
 k = klo
 GO TO 1140
 1130 k = khi
 1140 klo = khi
 GO TO 1080
 1150 a(1) = z(2*k)
 GO TO ndx,  (230,340,480,490,500)
 1160 CALL mesage (30,ierr,a(1))
 nogo = 1
 GO TO nerr, (220,330,485,495,505)
 1165 CALL mesage (-61,0,0)
 
 
!     FATAL ERROR MESAGES
 
 1170 ndx = -1
 GO TO 1210
 1180 ndx = -2
 GO TO 1210
 1190 CALL mesage (-30,ierr,spoint)
 1200 ndx = -3
 GO TO 1210
 1210 CALL mesage (ndx,FILE,gp1ah)
 1220 spoint(1) = z(k+3)
 spoint(2) = z(ii )
 ierr = 3
 GO TO 1190
 1230 spoint(1) = z(ii  )
 spoint(2) = z(ii+3)
 ierr = 4
 GO TO 1190
 1240 spoint(1) = 0
 spoint(2) = 0
 ierr = 5
 GO TO 1190
 1250 spoint(1) = a(1)
 spoint(2) = 0
 ierr = 12
 GO TO 1190
END SUBROUTINE gp1
