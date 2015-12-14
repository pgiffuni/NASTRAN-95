SUBROUTINE dpd1
     
!     DPD1 GENERATES THE GRID POINT LIST-DYNAMICS (GPLD),
!     USET-DYNAMICS (USETD), AND THE SCALAR INDEX LIST-DYNAMICS(SILD).
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        andf  ,orf
 LOGICAL :: nodyn ,first
 DIMENSION       buf(24)  ,epoint(2) ,seqep(2)     ,mcb(7)       ,  &
     nam(2)   ,loads(32) ,dload(2)     ,freq1(2)     ,  &
     freq(2)  ,nolin(21) ,tic(2)       ,  &
     tstep(2) ,tf(2)     ,psd(2)       ,msg(3)       ,  &
     eigr(2)  ,eigb(2)   ,eigc(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / luset ,lusetd,notfl ,nodlt ,nopsdl,nofrl ,nonlft,  &
     notrl ,noeed ,nosdt ,nbrep
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /dpdcom/ dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd ,  &
     dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed   ,  &
     scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn,  &
     loads ,dload ,freq1 ,freq  ,nolin ,nogo  ,  &
     msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  ,  &
     eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineq
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua    ,  &
     uf    ,us    ,un    ,ug    ,ue    ,up    ,une   , ufe   ,ud
 COMMON /zzzzzz/ z(1)
 COMMON /two   / two(32)
 COMMON /system/ sysbuf,iout  ,ks(37),nbpw
 EQUIVALENCE     (msg(2),ngrid)
 
 
!     SET NODYN FLAG TO TRUE IF NO DYNAMIC
 
 nodyn  = .false.
 buf(1) =  dpool
 CALL rdtrl (buf)
 IF (buf(1) /= dpool) nodyn = .true.
 
!     COMPUTE MAXIMUM EPOINT SIZE ALLOWED BY A COMPUTER WORD
 
 first = .true.
 IF (nodyn) GO TO 1020
 imax = 100000000
 IF (nbpw == 32) imax =  2147493
 IF (nbpw == 36) imax = 34359738
!         2147493=2**31/1000   34359738=2**35/1000
 maxz = imax
 mult = 1000
 
!     READ SECOND RECORD OF THE GPL INTO CORE. CREATE TABLE OF TRIPLES -
!     EXTERNAL GRID NO., SEQ. NO., AND INTERNAL GRID NO.
 
!     SEQ.NO.= EXTERNAL GIRD NO. * MULT, OR RESEQUENCED GRID PT. NO.
!     (A MULTIFICATION FACTOR WAS SAVED IN GPL HEADER RECORD BY GP1,
!      MULT = 10,100,OR 1000,  AND BY SGEN, MULT = 1000)
 
 1020 FILE = gpl
 IF (luset == 0) GO TO 1023
 CALL OPEN (*1023,gpl,z(buf1),rdrew)
 CALL READ (*2002,*2001,gpl,z(1),3,1,falg)
 CALL fwdrec (*2002,gpl)
 i = 3
 mult = z(i)
 imax = (imax/mult)*1000
 maxz = imax
 igpl = 1
 j = 1
 i = igpl
 1021 CALL READ (*2002,*1022,gpl,z(i),2,0,flag)
 z(i+2) = j
 i = i + 3
 j = j + 1
 GO TO 1021
 1022 ngpl = i - 3
 CALL CLOSE (gpl,clsrew)
 GO TO 1030
 
!     INITIALIZE FOR CASE WHERE NO GRID OR SCALAR PTS EXIST.
 
 1023 i = 1
 igpl  = 1
 luset = 0
 
!     READ EXTRA POINTS (IF ANY). ADD TO TABLE IN CORE.
!     SET INTERNAL GRID NO. OF EXTRA PTS = 0.
 
 1030 IF (nodyn) GO TO 1047
 FILE = dpool
 CALL preloc (*2001,z(buf1),dpool)
 iep  = i
 noep = 0
 CALL locate (*1045,z(buf1),epoint,flag)
 noep = 1
 1031 CALL READ (*2002,*1032,dpool,z(i),1,0,flag)
 IF (z(i) > maxz) maxz = z(i)
 z(i+1) = mult*z(i)
 z(i+2) = 0
 i = i + 3
 GO TO 1031
 1032 nep  = i - 3
 ngpl = nep
 
!     ONE OR MORE EPOINT WITH VERY LARGE EXTERNAL ID
!     FATAL IF MULTIPLIER IS 10
!     IF MULT IS 1000 OR 100, TRY TO SHRINK THE GRID POINT SEQ. NO. BY
!     10 OR 100 IF POSSIBLE, AND RESET MULT.
!     IF IT IS NOT POSSIBLE, WE HAVE A FATAL CONDITION 2140C
 
 IF (maxz == imax) GO TO 1040
 j = 0
 IF (mult ==   10) GO TO 1037
 mult = 100
 IF (maxz > 10*imax) mult = 10
 imax = (imax/mult)*1000
 j = 1000/mult
 DO  i = igpl,nep,3
   IF (MOD(z(i+1),j) /= 0) GO TO 1037
   z(i+1) = z(i+1)/j
 END DO
 GO TO 1040
 1037 WRITE  (iout,1038) ufm
 1038 FORMAT (a23,' 2140C, ONE OR MORE EPOINTS WITH  EXTERNAL ID TOO ',  &
     'LARGE.')
 IF (j /= 0) WRITE (iout,1039)
 1039 FORMAT (/5X,'SUGGESTION - RE-RUN NASTRAN JOB WITH ALL THE EPOINT',  &
     ' EXTERNAL ID''S SMALLER THAN THE LARGEST GRID POINT ID',  &
 /5X,'OR, REDUCE THE SEQID LEVEL IF SEQGP CARDS WERE USED',  &
     '.  I.E. FROM XXX.X.X TO XXX.X OR XXX')
 CALL CLOSE (dpool,clsrew)
 CALL mesage (-37,0,nam)
 
!     IF EXTRA POINTS PRESENT, READ SEQEP DATA (IF ANY).
!     REPLACE OLD SEQ NO WITH NEW SEQ NO.
 
 1040 CALL locate (*1045,z(buf1),seqep,flag)
 n1 = i
 n2 = n1 + 1
 ifail = 0
 2010 CALL READ (*2002,*2020,dpool,z(n2),buf1-1,1,flag)
 ifail = ifail + 1
 GO TO 2010
 2020 IF (ifail == 0) GO TO 2060
 nwds = (ifail-1)*(buf1-1) + flag
 WRITE  (iout,2040) ufm,nwds
 2040 FORMAT (a23,' 3139, UNABLE TO PROCESS SEQEP DATA IN SUBROUTINE ',  &
     'DPD1 DUE TO INSUFFICIENT CORE.', //5X,  &
     'ADDITIONAL CORE REQUIRED =',i10,7H  words)
 CALL mesage (-61,0,0)
 
!     CHECK FOR MULTIPLE REFERENCES TO EXTRA POINT ID NOS. AND
!     SEQUENCE ID NOS. ON SEQEP CARDS
 
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
     2100 FORMAT (a23,' 3140, MULTIPLE REFERENCES TO EXTRA POINT ID NO.',i9,  &
         ' ON SEQEP CARDS.')
     GO TO 2260
     2110 idseq1 = z(i)/1000
     irmndr = z(i) - 1000*idseq1
     IF (irmndr /= 0 .AND. mult >= 10) GO TO 2140
     WRITE  (iout,2120) ufm,idseq1
     2120 FORMAT (a23,' 3141, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,6X,  &
         'ON SEQEP CARDS.')
     GO TO 2260
     2140 idseq2 = irmndr/100
     irmndr = irmndr - 100*idseq2
     IF (irmndr /= 0 .AND. mult >= 100) GO TO 2180
     WRITE  (iout,2160) ufm,idseq1,idseq2
     2160 FORMAT (a23,' 3141, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,6X,'ON SEQEP CARDS.')
     GO TO 2260
     2180 idseq3 = irmndr/10
     irmndr = irmndr - 10*idseq3
     IF (irmndr /= 0) GO TO 2220
     WRITE  (iout,2200) ufm,idseq1,idseq2,idseq3
     2200 FORMAT (a23,' 3141, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,1H.,i1,4X,'ON SEQEP CARDS.')
     GO TO 2260
     2220 WRITE  (iout,2240) ufm,idseq1,idseq2,idseq3,irmndr
     2240 FORMAT (a23,' 3141, MULTIPLE REFERENCES TO SEQUENCE ID NO.',i6,  &
         1H.,i1,1H.,i1,1H.,i1,'  ON SEQEP CARDS.')
     2260 z(j) = -z(j)
   END DO
   
   2275 IF (jj < kk .OR. mult == 1 .OR. mult == 1000) CYCLE
   l = z(i)
   IF (mult   ==   10) GO TO 2280
   IF (MOD(l,10) /= 0) GO TO 2276
   z(i) = l / 10
   CYCLE
   2276 IF (.NOT.first) CYCLE
   first = .false.
   nogo  = 1
   WRITE  (iout,2277) ufm
   2277 FORMAT (a23,' 2140B, ILLEGAL DATA IN SEQEP CARD, POSSIBLY CAUSED',  &
       ' BY LARGE GRID OR SCALAR POINTS')
   CYCLE
   2280 IF (MOD(l,100) /= 0) GO TO 2276
   z(i) = l / 100
 END DO
 
 IF (k /= n2) GO TO 2290
 jj = kk
 k  = k + 1
 GO TO 2080
 
 2290 DO  i = n2,kk,2
   IF (z(i) < 0) z(i) = -z(i)
 END DO
 IF (nogo == 1) GO TO 2400
 
!     CHECK TO SEE IF ANY SEQUENCE ID NO. ON SEQEP CARDS IS THE SAME
!     AS AN EXTRA POINT ID NO. THAT HAS NOT BEEN RESEQUENCED
 
 loop2390:  DO  i = k,kk,2
   IF (z(i) < 0) CYCLE loop2390
   idseq1 = z(i) / mult
   irmndr = z(i) - mult*idseq1
   IF (irmndr /= 0) CYCLE loop2390
   DO  j = n2,kk,2
     IF (idseq1 == z(j)) CYCLE loop2390
   END DO
   DO  j = 1,n1,3
     IF (idseq1 == z(j)) GO TO 2360
   END DO
   CYCLE loop2390
   2360 nogo = 1
   WRITE  (iout,2380) ufm,idseq1
   2380 FORMAT (a23,' 3142, SEQUENCE ID NO.',i6,  &
       '  ON SEQEP CARDS IS THE SAME AS AN ', /5X,  &
       'EXTRA POINT ID NO. THAT HAS NOT BEEN RESEQUENCED.')
 END DO loop2390
 2400 CONTINUE
 i = -1
 1043 i = i + 2
 IF (i > flag) GO TO 1045
 buf(1) = z(n2+i-1)
 buf(2) = z(n2+i  )
 DO  j = iep,nep,3
   IF (z(j) == buf(1)) GO TO 1042
 END DO
 1044 buf(2) = 0
 CALL mesage (30,64,buf)
 nogo = 1
 GO TO 1043
 1042 IF (z(j+2) /= 0) GO TO 1044
 z(j+1) = buf(2)
 GO TO 1043
 1045 CALL CLOSE (dpool,clsrew)
 1047 IF (luset+noep == 0) GO TO 2004
 
!     IF EXTRA POINTS PRESENT, SORT THE GPL ON SEQ NO.
!     REPLACE SEQ NO WITH INTERNAL GRID NO FOR DYNAMICS.
 
 n = ngpl + 2
 IF (noep /= 0) CALL sort (0,0,3,2,z,n)
 i   = 2
 z(i)= 1
 IF (ngpl == 1) GO TO 1060
 DO  i = 4,ngpl,3
   z(i+1) = z(i-2) + 1
 END DO
 
!     WRITE THE GPLD.
 
 1060 FILE = gpld
 CALL OPEN  (*2001,gpld,z(buf1),wrtrew)
 CALL fname (gpld,buf)
 CALL WRITE (gpld,buf,2,1)
 DO  i = igpl,ngpl,3
   CALL WRITE (gpld,z(i),1,0)
 END DO
 CALL WRITE (gpld,0,0,1)
 CALL CLOSE (gpld,clsrew)
 mcb(1) = gpld
 mcb(2) = n/3
 CALL wrttrl (mcb)
 kn= mcb(2)
 
!     OPEN SILD AND USETD. WRITE HEADER RECORDS.
!     OPEN SIL  AND USET.  SKIP  HEADER RECORD.
!     READ SIL INTO CORE.
 
 FILE = sild
 CALL OPEN  (*2001,sild,z(buf1),wrtrew)
 CALL fname (sild,buf)
 CALL WRITE (sild,buf,2,1)
 IF (luset == 0) GO TO 1082
 FILE = sil
 CALL OPEN (*2001,sil,z(buf2),rdrew)
 CALL fwdrec (*2002,sil)
 isil = ngpl + 3
 CALL READ (*2002,*1081,sil,z(isil),buf3-isil,1,n)
 CALL mesage (-8,0,nam)
 1081 CALL CLOSE  (sil,clsrew)
 nsil = isil + n
 z(nsil)= luset + 1
 1082 FILE = usetd
 CALL OPEN  (*2001,usetd,z(buf3),wrtrew)
 CALL fname (usetd,buf)
 CALL WRITE (usetd,buf,2,1)
 IF (luset == 0) GO TO 1100
 FILE = uset
 CALL OPEN (*2001,uset,z(buf2),rdrew)
 CALL fwdrec (*2002,uset)
 
!     INITIALIZE DISPLACEMENT SET BIT MASKS.
 
 1100 i = igpl
 j = isil - 1
 nbrep = 0
 buf(10) = 1
 DO  k = 2,7
   mcb(k) = 0
 END DO
 mskua  = two(ua)
 mskun  = two(un)
 mskuf  = two(uf)
 mskue  = two(ue)
 mskup  = two(up)
 mskud  = two(ud)
 mskune = two(une)
 mskufe = two(ufe)
 musetd = orf(mskue,orf(mskune,orf(mskufe,orf(mskud,mskup))))
 
!     TEST FOR CURRENT POINT IN G-SET OR IN P-SET (EXTRA POINT).
 
 1110 IF (z(i+2) == 0) GO TO 1130
 
!     POINT IS IN G-SET - READ USET MASKS BELONGING TO POINT.
!     TURN ON APPROPRIATE BITS FOR P-SET. WRITE MASKS ON USETD.
 
 j = j + 1
 m = z(j+1) - z(j)
 CALL READ (*2002,*2003,uset,buf,m,0,flag)
 DO  k = 1,m
   ksw = orf(buf(k),mskup)
   IF (andf(ksw,mskua) /= 0) ksw = orf(ksw,mskud )
   IF (andf(ksw,mskun) /= 0) ksw = orf(ksw,mskune)
   IF (andf(ksw,mskuf) /= 0) ksw = orf(ksw,mskufe)
   mcb(5) = orf(mcb(5),ksw)
   buf(k) = ksw
 END DO
 CALL WRITE (usetd,buf,m,0)
 GO TO 1140
 
!     POINT IS AN EXTRA POINT - WRITE MASK ON USETD.
 
 1130 CALL WRITE (usetd,musetd,1,0)
 mcb(5) = orf(mcb(5),musetd)
 m = 1
 
!     REPLACE INTERNAL DYNAMICS NO. WITH SILD NO. WRITE SILD ENTRY.
!     REPLACE INTERNAL STATICS NO. WITH SIL NO.
 
 1140 z(i+1) = buf(10)
 CALL WRITE (sild,z(i+1),1,0)
 IF (z(i+2) == 0) GO TO 1141
 z(i+2) = z(j)
 GO TO 1150
 1141 nbrep = nbrep + 1
 
!     TEST FOR COMPLETION.
 
 1150 buf(10) = buf(10) + m
 i = i + 3
 IF (i <= ngpl) GO TO 1110
 
!     WRITE SECOND RECORD OF SILD (PAIRS OF SIL NO., SILD NO.)
 
 CALL WRITE (sild,0,0,1)
 CALL WRITE (usetd,0,0,1)
 DO  i = igpl,ngpl,3
   IF (z(i+2) == 0) CYCLE
   buf(1) = z(i+2)
   buf(2) = z(i+1)
   CALL WRITE (sild,buf,2,0)
 END DO
 
!     CLOSE FILES AND WRITE TRAILERS.
 
 CALL CLOSE (sild ,clsrew)
 CALL CLOSE (usetd,clsrew)
 mcb(1) = sild
 lusetd = luset + nbrep
 mcb(2) = lusetd
 mcb(3) = nbrep
 CALL wrttrl (mcb)
 mcb(1) = usetd
 CALL wrttrl (mcb)
 mcb(5) = 0
 CALL CLOSE (uset,clsrew)
 
!     REPLACE SIL NO. IN TABLE WITH CODED SILD NO.
!     THEN SORT TABLE ON EXTERNAL GRID NO.
 
 z(ngpl+4) = lusetd + 1
 DO  i = igpl,ngpl,3
   j = 1
   IF (z(i+4)-z(i+1) /= 1) GO TO 1091
   j = 2
   IF (z(i+2) == 0) j = 3
   1091 z(i+2) = 10*z(i+1) + j
 END DO
 CALL sort (0,0,3,1,z(igpl),ngpl-igpl+3)
 
!     WRITE EQDYN DATA BLOCK. FIRST RECORD IS PAIRS OF EXTERNAL GRID NO,
!     SILD NO. SECOND RECORD IS PAIRS OF EXTERNAL GRID NO., CODED SILD
!     NO.
 
 FILE = eqdyn
 CALL OPEN  (*2001,eqdyn,z(buf1),wrtrew)
 CALL fname (eqdyn,buf)
 CALL WRITE (eqdyn,buf,2,1)
 DO  i = igpl,ngpl,3
   CALL WRITE (eqdyn,z(i),2,0)
 END DO
 CALL WRITE (eqdyn,0,0,1)
 DO  i = igpl,ngpl,3
   buf(1) = z(i  )
   buf(2) = z(i+2)
   CALL WRITE (eqdyn,buf,2,0)
 END DO
 CALL WRITE (eqdyn,0,0,1)
 CALL CLOSE (eqdyn,clsrew)
 mcb(1) = eqdyn
 mcb(2) = kn
 CALL wrttrl (mcb)
 neqdyn = 2*kn - 1
 IF (nbrep == 0) nbrep = -1
 RETURN
 
!     FATAL FILE ERRORS
 
 2001 n = -1
 GO TO 2005
 2002 n = -2
 GO TO 2005
 2003 n = -3
 GO TO 2005
 2004 n = -30
 FILE = 109
 2005 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE dpd1
