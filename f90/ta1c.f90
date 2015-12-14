SUBROUTINE ta1c
     
!     TA1C READS GENERAL ELEMENTS FROM THE ECT AND BUILDS THE GEI.
!     FOR EACH GENERAL ELEMENT, THE UI AND UD LISTS ARE CONVERTED TO
!     SIL NOS. AND SORTED ON SIL NO. THE ELEMENTS OF THE Z AND S
!     MATRICES ARE WRITTEN IN INTERNAL SORT (I.E., ROW AND COL NOS
!     CORRESPOND TO POSITION IN THE SORTED UI AND UD LISTS.
 
 
 INTEGER :: genl  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     est   ,gpect ,gei   ,ecpt  ,gpct  ,scr1  ,scr2  ,  &
     scr3  ,scr4  ,z     ,sysbuf,buf1  ,buf2  ,buf3  ,  &
     FILE  ,flag  ,genel ,rd    ,rdrew ,wrt   ,wrtrew, clsrew,silno ,buf   ,half
 DIMENSION       nam(2),buf(10)      ,genel(2)
 COMMON /BLANK / luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /ta1com/ nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4
 COMMON /tac1ax/ buf1  ,buf2  ,buf3  ,iui   ,nui   ,iud   ,nud   ,  &
     iz    ,nogo  ,idgenl
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,dum38(38)    ,nbpw
 COMMON /setup / nfile(6)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 DATA    genel / 4301  ,43  / ,nam /  4HTA1C,4H    /
 DATA    half  / 65536      /
 
!     ADD MORE BITS TO HALF IF MACHINE WORD IS LARGER THAN 32
 
 IF (nbpw >= 36) half = 4*half
 IF (nbpw > 36) half = 4*half
 
!     SET BUFFER POINTERS, ETC.
 
 buf1 = korsz(z) - sysbuf - 2
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 nogo = 0
 nogenl = 0
 
!     READ THE SIL INTO CORE
 
 FILE = sil
 CALL OPEN (*2001,sil,z(buf1),rdrew)
 CALL fwdrec (*2002,sil)
 CALL READ (*2002,*1011,sil,z,buf2,1,nsil)
 CALL mesage (-8,0,nam)
 1011 CALL CLOSE (sil,clsrew)
 
!     OPEN THE GEI. WRITE HEADER RECORD.
 
 FILE = gei
 CALL OPEN (*2001,gei,z(buf2),wrtrew)
 CALL fname (gei,buf)
 CALL WRITE (gei,buf,2,1)
 
!     OPEN THE ECT. READ ELEMENT ID.
 
 FILE = ect
 CALL preloc (*2001,z(buf1),ect)
 CALL locate (*2006,z(buf1),genel,flag)
 1031 CALL READ (*2002,*1150,ect,buf,1,0,flag)
 idgenl = buf(1)
 nogenl = nogenl + 1
 
!     READ THE UI LIST. STORE POSITION IN UI LIST, SIL NO.,
!     INTERNAL GRID NO., AND COMPONENT CODE.
 
 iui = nsil + 1
 i = iui
 j = 1
 1041 CALL READ (*2002,*2003,ect,z(i+2),2,0,flag)
 IF (z(i+2) == -1) GO TO 1042
 z(i) = j
 k = z(i+2)
 z(i+1) = z(k)
 IF (z(i+3) /= 0) z(i+1) = z(i+1) + z(i+3) - 1
 i = i + 4
 j = j + 1
 GO TO 1041
 1042 nui   = i - 4
 nbrui = j - 1
 nwdui = 4*nbrui
 
!     READ THE UD LIST (IF PRESENT). STORE POSITION IN UD LIST, SIL NO.,
!     INTERNAL GRID NO., AND COMPONENT CODE.
 
 iud = i
 j = 1
 1051 CALL READ (*2002,*2003,ect,z(i+2),2,0,flag)
 IF (z(i+2) == -1) GO TO 1052
 z(i) = j
 k = z(i+2)
 z(i+1) = z(k)
 IF (z(i+3) /= 0) z(i+1) = z(i+1) + z(i+3) - 1
 i = i + 4
 j = j + 1
 GO TO 1051
 1052 nud   = i - 4
 nbrud = j - 1
 nwdud = 4*nbrud
 iz = i
 
!     SORT UI AND UD LISTS ON SIL NO.
!     STORE INTERNAL POSITION IN UI AND UD LISTS.
!     WRITE ELEMENT ID, NO. OF UI-S, NO. OF UD-S.
!     WRITE SIL NOS. FOR UI LIST AND SIL NOS. FOR UD LIST.
 
 CALL sorti (0,0,4,2,z(iui),nwdui)
 buf(2) = nbrui
 buf(3) = nbrud
 CALL WRITE (gei,buf,3,0)
 k = 1
 DO  i = iui,nui,4
   silno  = z(i+1)
   z(i+1) = k
   CALL WRITE (gei,silno,1,0)
   k = k + 1
 END DO
 IF (nbrud == 0) GO TO 1070
 CALL sorti (0,0,4,2,z(iud),nwdud)
 k = 1
 DO  i = iud,nud,4
   silno  = z(i+1)
   z(i+1) = k
   CALL WRITE (gei,silno,1,0)
   k = k + 1
 END DO
 
!     SORT UI LIST ON EXTERNAL POSITION.
 
 1070 CALL sorti (0,0,4,1,z(iui),nwdui)
 
!     DETERMINE IF CORE WILL HOLD THE FULL Z OR K MATRIX
 
 ncore  = buf2 - iz
 nwdz   = nbrui**2
 nocore = 0
 IF (nwdz > ncore) nocore = 1
 
!     READ INDICATOR OF INPUT OF Z OR K MATRIX
 
 CALL READ (*2002,*2003,ect,ijk,1,0,flag)
 CALL WRITE (gei,ijk,1,0)
 koz = 0
 IF (ijk == 2) koz = 1
 
!     READ THE ELEMENTS OF THE Z OR K MATRIX.
!     CONVERT FROM EXTERNAL ROW AND COL NOS. TO INTERNAL ROW AND COL
!     NOS.  IF CORE WILL HOLD Z OR K, STORE THE ELEMENTS IN CORE
!     OTHERWISE, WRITE CODED ROW/COL NOS AND ELEMENTS ON SCRATCH FILE.
 
 IF (nocore /= 0) CALL OPEN (*2001,scr4,z(buf3),wrtrew)
 DO  i = iui,nui,4
   introw = z(i+1)
   krow = iz + (introw-1)*nbrui - 1
   DO  j = i,nui,4
     intcol = z(j+1)
     kcol = iz + (intcol-1)*nbrui - 1
     CALL READ (*2002,*2003,ect,buf(3),1,0,flag)
     IF (nocore /= 0) GO TO 1092
     k    = krow + intcol
     z(k) = buf(3)
     k    = kcol + introw
     z(k) = buf(3)
     GO TO 1093
     1092 m = 3
     buf(1) = intcol
     buf(2) = introw
     IF (introw == intcol) GO TO 1095
     buf(4) = introw
     buf(5) = intcol
     buf(6) =buf(3)
     m = 6
     1095 CALL WRITE (scr4,buf,m,0)
     1093 CONTINUE
   END DO
 END DO
 IF (nocore /= 0) CALL CLOSE (scr4,clsrew)
 
!     IF Z OR K MATRIX IS IN CORE,WRITE IT OUT
!     OTHERWISE,SORT THE MATRIX AND THEN WRITE IT.
 
 IF (nocore == 0) GO TO  1103
 CALL OPEN (*2001,scr4,z(buf3),rdrew)
 nfile(1) = scr1
 nfile(2) = scr2
 nfile(3) = scr3
 CALL sorti (scr4,0,3,2,z(iz),ncore-sysbuf)
 CALL CLOSE (scr4,clsrew)
 IF (nfile(6) == nfile(1)) nfile(1) = scr4
 IF (nfile(6) == nfile(2)) nfile(2) = scr4
 IF (nfile(6) == nfile(3)) nfile(3) = scr4
 jfile = nfile(6)
 CALL OPEN (*2001, jfile, z(buf3), rdrew)
 CALL sorti (jfile, 0, 3, -1, z(iz), ncore-sysbuf)
 CALL CLOSE (jfile, clsrew)
 CALL OPEN  (*2001,nfile(6),z(buf3),rdrew)
 1101 CALL READ  (*2002,*1102,nfile(6),buf,3,0,flag)
 CALL WRITE (gei,buf(3),1,0)
 GO TO 1101
 1102 CALL CLOSE (nfile(6),clsrew)
 GO TO 1110
 1103 CALL WRITE (gei,z(iz),nwdz,0)
 
!     READ FLAG WORD FOR S MATRIX.
!     IF S MATRIX NOT PRESENT, BUT UD IS PRESENT,
!     EXECUTE TA1CA TO COMPUTE AND WRITE S MATRIX.
!     IF S MATRIX AND UD BOTH NOT PRESENT, CLOSE GEI RECORD AND LOOP
!     BACK
 
 1110 CALL READ (*2002,*2003,ect,buf,1,0,flag)
 IF (buf(1) /= 0) GO TO 1120
 IF (nbrud  == 0) GO TO 1111
 CALL sorti (0,0,4,2,z(iui),nwdui)
 CALL ta1ca (koz)
 1111 CALL WRITE (gei,0,0,1)
 GO TO 1031
 
!     S MATRIX IS PRESENT.
!     DETERMINE IF CORE WILL HOLD THE FULL S MATRIX
 
 1120 nwds = nbrud*nbrui
 CALL sorti (0,0,4,1,z(iud),nwdud)
 nocore = 0
 IF (nwds > ncore) nocore = 1
 
!     READ THE ELEMENTS OF THE S MATRIX.
!     CONVERT FROM EXTERNAL ROW AND COL NOS TO INTERNAL ROW AND COL NOS.
!     IF CORE WILL HOLD S, STORE THE ELEMENTS IN CORE.
!     OTHERWISE, WRITE CODED ROW/COL NOS AND ELEMENTS ON SCRATCH FILE.
 
 IF (nocore /= 0) CALL OPEN (*2001,scr4,z(buf3),wrtrew)
 DO  i = iui,nui,4
   introw = z(i+1)
   krow   = iz + (introw-1)*nbrud - 1
   DO  j = iud,nud,4
     intcol = z(j+1)
     k = krow + intcol
     CALL READ (*2002,*2003,ect,buf(3),1,0,flag)
     IF (nocore /= 0) GO TO 1131
     z(k) = buf(3)
     CYCLE
     1131 buf(1) = introw
     buf(2) = intcol
     CALL WRITE (scr4,buf,3,1)
   END DO
 END DO
 IF (nocore /= 0) CALL CLOSE (scr4,clsrew)
 
!     IF S MATRIX IS IN CORE, WRITE IT OUT.
!     OTHERWISE, SORT THE MATRIX AND THEN WRITE IT.
 
 IF (nocore == 0) GO TO 1142
 CALL OPEN (*2001,scr4,z(buf3),rdrew)
 nfile(1) = scr1
 nfile(2) = scr2
 nfile(3) = scr3
 CALL sorti (scr4,0,3,2,z(iz),ncore-sysbuf)
 CALL CLOSE (scr4,clsrew)
 IF (nfile(6) == nfile(1)) nfile(1) = scr4
 IF (nfile(6) == nfile(2)) nfile(2) = scr4
 IF (nfile(6) == nfile(3)) nfile(3) = scr4
 jfile = nfile(6)
 CALL OPEN (*2001, jfile, z(buf3), rdrew)
 CALL sorti (jfile, 0, 3, -1, z(iz), ncore-sysbuf)
 CALL CLOSE (jfile, clsrew)
 CALL OPEN  (*2001,nfile(6),z(buf3),rdrew)
 1141 CALL READ  (*2002,*1143,nfile(6),buf,3,0,FILE)
 CALL WRITE (gei,buf(3),1,0)
 GO TO 1141
 1142 CALL WRITE (gei,z(iz),nwds,0)
 1143 CALL WRITE (gei,0,0,1)
 GO TO 1031
 
!     HERE WHEN NO MORE GENERAL ELEMENTS
 
 1150 CALL CLOSE (ect,clsrew)
 CALL CLOSE (gei,clsrew)
 buf(1) = gei
 buf(2) = nogenl
 CALL wrttrl (buf)
 IF (nogo /= 0) CALL mesage (-61,0,nam)
 RETURN
 
!     FATAL ERRORS
 
 2001 n = -1
 GO TO 2005
 2002 n = -2
 GO TO 2005
 2003 n = -3
 2005 CALL mesage (n,FILE,nam)
 2006 CALL mesage (-30,63,buf)
 RETURN
END SUBROUTINE ta1c
