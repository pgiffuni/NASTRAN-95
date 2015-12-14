SUBROUTINE reduce
     
!     REDUCE BUILDS THE FOLLOWING DATA BLOCKS
 
!     1.  PVX  -  THE REDUCTION PARTITIONING VECTOR
!     2.  USX  -  THE USET EQUIVALENT VECTOR
!     3.  INX  -  THE REDUCTION TRANSFORMATION IDENTITY PARTITION
 
!     THE FOLLOWING BULK DATA CARDS ARE READ
 
!     1.  BDYC
!     2.  BDYS
!     3.  BDYS1
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: inbset,fset,bad,lonly
 REAL :: rz(1)
 DIMENSION       modnam(2),ijk(6),ihd(96),bdys(2),bdys1(2),  &
     bdyc(2),mnem(4),namold(14),namnew(2),aray(6),  &
     isid(100),cset(6),ipset(6),listo(32),listn(32), mcb(7),ibits(6)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /zzzzzz/ z(1)
 COMMON /packx / typin,typout,irow,nrow,incr
 COMMON /system/ sysbuf,outt,x1(6),nlpp,x2(2),line,x3(2),idate(3)
 COMMON /output/ ititl(96),ihead(96)
 COMMON /cmbfnd/ iinam(2),iiierr
 COMMON /two   / tpow(32)
 COMMON /BLANK / step,dry,pora
 EQUIVALENCE     (rz(1),z(1) )
 DATA    ihd   / 4H    , 8*4H****,  &
     4H s u, 4H b s, 4H t r, 4H u c, 4H t u, 4H r e, 4H    ,  &
     4HM o , 4HD u , 4HL e , 4H   r, 4H e d, 4H u c, 4H e *,  &
     9*4H**** , 64*4H      /
 DATA    nheqss, nhbgss,nhcstm,nhplts/4HEQSS,4HBGSS,4HCSTM,4HPLTS/
 DATA    modnam/ 4HREDU,4HCE   /
 DATA    papp  , lods, loap    /4HPAPP,4HLODS,4HLOAP/
!     --------------------
!     CODES TO LOCATE BULK DATA
!     --------------------
 DATA    bdyc/ 910,9 / , bdys/ 1210,12 / , bdys1/ 1310,13 /
!     --------------------
!     CASE CONTROL MNEMONICS
!     --------------------
 DATA    mnem/ 4HNAMA , 4HNAMB , 4HBOUN , 4HOUTP /
!     --------------------
!     GINO FILES FOR DATA BLOCKS AND SCRATCH
!     --------------------
 DATA    casecc/ 101 / , geom4/ 102 /
 DATA    pvx   / 201 / , usx  / 202 / , inx/ 203 /
 DATA    scr1  / 301 / , scr2 / 302 / , i3 / 3   /
 
 
!     I.  COMPUTE OPEN CORE AND DEFINE GINO AND SOF BUFFERS
!     *****************************************************
 
 IF (dry == -2) RETURN
 iba = 128
 ibo = 4
 ibf = 64
 nzwd= korsz(z(1))
 IF (nzwd <= 0 ) CALL mesage (-8,0,modnam)
 
 lonly= .false.
 buf1 = nzwd - sysbuf - 2
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 ib1  = buf3 - sysbuf
 ib2  = ib1  - sysbuf
 ib3  = ib2  - sysbuf
 
!     SCORE IS STARTING ADDRESS OF OPEN CORE AND NZ THE LENGTH
 
 score = 1
 nz = ib3 - 1
 
!     INITIALIZE ACTIVITY ON THE SOF
 
 litm = lods
 IF (pora == papp) litm = loap
 CALL sofopn (z(ib1),z(ib2),z(ib3))
 DO  i = 1,96
   ihead(i) = ihd(i)
 END DO
 
!     II.  PROCESS THE CASE CONTROL DATA BLOCK ( CASECC )
!     ***************************************************
 
 DO  i = 1,14
   namold(i) = 0
 END DO
 ifile = casecc
 CALL OPEN (*2001,casecc,z(buf2),0)
 prtopt = 0
 nrec = step
 IF (nrec == 0) THEN
   GO TO   201
 END IF
 200 DO  i = 1,nrec
   CALL fwdrec (*2002,casecc)
 END DO
 
!     BEGIN READING CASECC
 
 201 inbset = .false.
 CALL READ (*2002,*2003,casecc,z(1),2,0,nnn)
 nwdscc = z(i3-1)
 DO  i = 1,nwdscc,3
   CALL READ (*2002,*2003,casecc,z(1),3,0,nnn)
   
!     CHECK FOR CASE CONTROL MNEMONICS
   
   DO  j = 1,4
     IF (z(1) == mnem(j)) GO TO 205
   END DO
   CYCLE
   205 SELECT CASE ( j )
     CASE (    1)
       GO TO 206
     CASE (    2)
       GO TO 207
     CASE (    3)
       GO TO 208
     CASE (    4)
       GO TO 209
   END SELECT
   206 namold(1) = z(i3-1)
   namold(2) = z(i3  )
   CYCLE
   207 namnew(1) = z(i3-1)
   namnew(2) = z(i3  )
   CYCLE
   208 inbset = .true.
   bset = z(i3)
   CYCLE
   209 prtopt = orf(prtopt,z(i3))
 END DO
 IF (dry == 0) prtopt = 0
 IF (andf(prtopt,1) /= 1) GO TO 2199
 CALL page1
 WRITE  (outt,280) (namold(i),i=1,2),namnew,bset,(namold(i),i=1,2)
 280 FORMAT (//41X,'S U M M A R Y    O F    C U R R E N T    P R O ',  &
     'B L E M', //43X,  &
     'NAME OF PSEUDOSTRUCTURE TO BE REDUCED    - ',2A4, //43X,  &
     'NAME GIVEN TO RESULTANT PSEUDOSTRUCTURE  - ',2A4, //43X,  &
     'BOUNDARY SET IDENTIFICATION NUMBER       - ',i8,  //43X,  &
     'NAMES OF COMPONENT SUBSTRUCTURES CONTAINED IN ',2A4/)
 2199 CONTINUE
 CALL CLOSE (casecc,1)
 
!     CHECK FOR ALLOWABILITY OF INPUT
 
 bad = .false.
 CALL sfetch (namold,nheqss,3,itest)
 IF (itest == 4) GO TO 290
 261 CALL sfetch (namnew,nheqss,3,itest)
 IF (itest /= 4 .AND. dry /= 0) GO TO 291
 IF (itest == 4 .AND. dry == 0) GO TO 297
 262 IF (.NOT.inbset) GO TO 292
 263 IF (.NOT.bad) GO TO 300
 GO TO 2100
 
!     IF NO ERRORS, CONTINUE PROCESSING
 
 
 290 WRITE (outt,293) ufm,(namold(i),i=1,2)
 bad = .true.
 GO TO 261
 291 CALL sfetch (namnew,litm,3,itest)
 IF (itest /= 3) GO TO 296
 lonly = .true.
 GO TO 300
 296 CONTINUE
 WRITE (outt,294) ufm,(namnew(i),i=1,2)
 bad = .true.
 GO TO 262
 292 WRITE (outt,295) ufm
 bad = .true.
 GO TO 263
 297 WRITE  (outt,298) ufm,namnew
 298 FORMAT (a23,' 6613, FOR RUN=GO, THE REDUCED SUBSTRUCTURE ',2A4,  &
     ' MUST ALREADY EXIST.')
 bad = .true.
 GO TO 262
 293 FORMAT (a23,' 6601, REQUEST TO REDUCE PSEUDOSTRUCTURE ',2A4,  &
     ' INVALID. DOES NOT EXIST ON THE SOF.')
 294 FORMAT (a23,' 6602, THE NAME ',2A4,' CAN NOT BE USED FOR THE ',  &
     'REDUCED PSEUDOSTRUCTURE. IT ALREADY EXISTS ON THE SOF.')
 295 FORMAT (a23,' 6603, A BOUNDARY SET MUST BE SPECIFIED FOR A ',  &
     'REDUCE OPERATION.')
 
!     READ FIRST GROUP OF EQSS FOR THE STRUCTURE BEING REDUCED,
!     PLACE THE NAMES OF THE COMPONENT SUBSTRUCTURES INTO THE
!     FIRST NWDS WORDS OF OPEN CORE.
 
 300 ks1 = score
 CALL sfetch (namold,nheqss,1,itest)
 CALL suread (z(ks1),-1,nout,itest )
 
!     NCSUB IS THE NUMBER OF COMPONENT SUBSTRUCTURES
!     NIPOLD IS THE NUMBER OF IP S IN THE STRUCTURE BEING REDUCED
 
 ncsub = z(ks1+2)
 nout  = nout - 4
 DO  i = 1,nout
   ii = i - 1
   z(ks1+ii) = z(ks1+4+ii)
 END DO
 nwds = nout
 score= ks1 + nwds
 kf1  = score - 1
 nz   = nz - nwds
 IF (andf(prtopt,1) /= 1) GO TO 282
 WRITE  (outt,281) (z(jj),jj=ks1,kf1)
 281 FORMAT (48X,2A4,4X,2A4,4X,2A4,4X,2A4)
 282 CONTINUE
 
!     III. READ BOUNDARY SET ( BDYC ) BULK DATA INTO OPEN CORE FOR
!     THE REQUESTED SET ( BSET ) FROM THE GEOM4 INPUT DATA BLOCK.
!     ************************************************************
 
 ks2   = score
 ifile = geom4
 CALL preloc (*2001,z(buf1),geom4)
 CALL locate (*490,z(buf1),bdyc,flag)
 401 CALL READ (*2002,*490,geom4,id,1,0,nnn)
 IF (id == bset) GO TO 402
 403 CALL READ (*2002,*2003,geom4,aray,3,0,nnn)
 IF (aray(3) == -1) GO TO 401
 GO TO 403
 
!     CORRECT BOUNDARY SET HAS BEEN FOUND, STORE DATA IN SECOND NWBS WOR
!     OF OPEN CORE.
 
 402 nwbs = 0
 405 bad  = .false.
 CALL READ (*2002,*2003,geom4,z(ks2+nwbs),3,0,nnn)
 IF (z(ks2+nwbs+2) == -1) GO TO 440
 
!     MUST CHECK THAT THE SUBSTRUCTURE IS A PHASE1 BASIC SUBSTRUCTURE
!     AND THAT IT IS A COMPONENT OF THE STRUCTURE BEING REDUCED.
 
!     CHECK FOR COMPONENT
 
 DO  i = 1,nwds,2
   ii = i - 1
   IF (z(ks1+ii) == z(ks2+nwbs) .AND. z(ks1+ii+1) == z(ks2+nwbs+1)) GO TO 420
 END DO
 
!     NOT A COMPONENT
 
 WRITE (outt,491) ufm,z(ks2+nwbs),z(ks2+nwbs+1)
 bad = .true.
 491 FORMAT (a23,' 6604, A BOUNDARY SET HAS BEEN SPECIFIED FOR ',2A4,  &
     ', BUT IT IS NOT A COMPONENT OF THE', /31X,'PSEUDOSTRUC',  &
     'TURE BEING REDUCED. THE BOUNDARY SET WILL BE IGNORED.')
 
 420 IF (bad) GO TO 405
 nwbs = nwbs + 3
 GO TO 405
 440 score = ks2 + nwbs
 kf2 = score - 1
 nz  = nz - nwbs
 
!     SORT ON SET ID
 
 CALL sort (0,0,3,3,z(ks2),nwbs)
 IF (andf(rshift(prtopt,1),1) /= 1) GO TO 2299
 ii = 0
 2203 CALL page1
 WRITE  (outt,2202) bset
 2202 FORMAT (//44X,'SUMMARY OF COMBINED BOUNDARY SET NUMBER',i9, //55X,  &
     'BASIC',11X,'BOUNDARY', /52X,'SUBSTRUCTURE',8X,'SET ID',  &
     /56X,'NAME',12X,'NUMBER',/)
 line = line + 7
 2206 line = line + 1
 IF (line > nlpp) GO TO 2203
 WRITE  (outt,2205) z(ks2+ii),z(ks2+ii+1),z(ks2+ii+2)
 2205 FORMAT (54X,2A4,9X,i8)
 ii = ii + 3
 IF (ii > nwbs - 3) GO TO 2299
 GO TO 2206
 2299 CONTINUE
 GO TO 500
!WKBR 8/94 ALPHA-VMS  490 WRITE (OUTT,493) IFM,BSET
 490 WRITE (outt,493) ufm,bset
 GO TO 2200
 493 FORMAT (a23,' 6606, BOUNDARY SET ,I8,61H SPECIFIED IN CASE ',  &
     'CONTROL HAS NOT BEEN DEFINED BY BULK DATA.')
 
!     IV. READ BDYS BULK DATA PROCESSING ONLY THE SET ID S REFERENCED ON
!     THE BDYC CARD.  IF DATA DOES NOT EXIST, GO TO BDYS1 PROCESSING SEC
!     ******************************************************************
 
 500 j = 0
 ierr = 0
 CALL locate (*580,z(buf1),bdys,flag)
 502 CALL READ (*2002,*600,geom4,idhid,1,0,nnn)
 
!     CHECK REQUESTED ID
 
 DO  i = ks2,kf2,3
   IF (idhid == z(i+2)) GO TO 503
 END DO
 505 CALL READ (*2002,*2003,geom4,aray,2,0,nnn)
 IF (aray(1) /= -1 .AND. aray(2) /= -1) GO TO 505
 GO TO 502
 503 CALL READ (*2002,*2003,geom4,aray,2,0,nnn)
 IF (aray(1) == -1 .AND. aray(2) == -1) GO TO 502
 z(score+j  ) = idhid
 z(score+j+1) = aray(1)
 z(score+j+2) = aray(2)
 j = j + 3
 GO TO 503
 580 ierr = ierr + 1
 
!     V. READ BDYS1 BULK DATA AND MERGE WITH BDYS IN OPEN CORE.
!     *********************************************************
 
 600 CALL locate (*620,z(buf1),bdys1,flag)
 606 CALL READ (*2002,*602,geom4,aray(1),2,0,nnn)
 
!     CHECK ID
 
 DO  i = ks2,kf2,3
   IF (aray(1) == z(i+2)) GO TO 604
 END DO
 605 CALL READ (*2002,*2003,geom4,aray(3),1,0,nnn)
 IF (aray(3) /= -1) GO TO 605
 GO TO 606
 604 CALL READ (*2002,*2003,geom4,aray(3),1,0,nnn)
 IF (aray(3) == -1) GO TO 606
 z(score+j  ) = aray(1)
 z(score+j+1) = aray(3)
 z(score+j+2) = aray(2)
 j = j + 3
 GO TO 604
 620 ierr = ierr + 1
 602 CALL CLOSE (geom4,1)
 IF (ierr /= 2) GO TO 650
 WRITE (outt,691) ufm,bset
 GO TO 2200
 691 FORMAT (a23,' 6607, NO BDYS OR BDYS1 BULK DATA HAS BEEN INPUT TO',  &
     ' DEFINE BOUNDARY SET',i8)
 
!     SORT COMPLETE BOUNDARY SET DATA ON SET ID IN OPEN CORE
 
 650 CALL sort (0,0,3,1,z(score),j)
 
!     TRANSLATE COMPONENT NUMBER TO BIT PATTERN
 
 it = score + j - 1
 DO  i = score,it,3
   CALL encode (z(i+2))
 END DO
 IF (andf(rshift(prtopt,2),1) /= 1) GO TO 2399
 iinc = 0
 2303 CALL page1
 WRITE  (outt,2302)
 2302 FORMAT (1H0,46X,44HTABLE of grid points composing boundary sets, /  &
     /52X,8HBOUNDARY ,/52X , 34H set id      grid point       dof ,  &
     /52X,34H NUMBER      id  NUMBER       code ,/ )
 line = line + 7
 2305 line = line + 1
 IF (line > nlpp) GO TO 2303
 icode = z(score+iinc+2)
 CALL bitpat (icode, ibits)
 WRITE (outt,2304) z(score+iinc),z(score+iinc+1), ibits(1),ibits(2)
 2304 FORMAT (52X,i8,6X,i8,7X,a4,a2)
 iinc = iinc + 3
 IF (iinc > j-3) GO TO 2399
 GO TO 2305
 2399 CONTINUE
 
!     WRITE BOUNDARY SET DATA ON TO FILE SCR1, ONE LOGICAL RECORD FOR EA
!     SET ID.
 
 CALL OPEN (*2001,scr1,z(buf2),1)
 ist  = score + 3
 ifin = score + j - 1
 n    = 1
 nsid = 1
 isid(1) = z(score)
 CALL WRITE (scr1,z(score+1),2,0)
 DO  i = ist,ifin,3
   IF (z(i) == isid(n)) GO TO 661
   n    = n + 1
   nsid = nsid + 1
   isid(n) = z(i)
   CALL WRITE (scr1,aray,0,1)
   661 CALL WRITE (scr1,z(i+1),2,0)
 END DO
 CALL WRITE (scr1,aray,0,1)
 CALL CLOSE (scr1,1)
 
 
!     SCR1 NOW CONTAINS BOUNDARY SET DATA FOR ALL GRID POINTS
 
!     CHECK THAT ALL REQUESTED SID S HAVE BEEN FOUND
 
 nrsid = nwbs/3
 j = 0
 DO  i = ks2,kf2,3
   z(score+j) = z(i+2)
   j = j + 1
 END DO
 DO  i = 1,nrsid
   ii = i - 1
   DO  j = 1,nsid
     IF (isid(j) == z(score+ii)) GO TO 677
   END DO
   CYCLE
   677 z(score+ii) = 0
 END DO
 ibad = 0
 DO  i = 1,nrsid
   ii = i - 1
   IF (z(score+ii) == 0) CYCLE
   INDEX = (i-1)*3
   WRITE (outt,692) ufm,z(ks2+INDEX+2),z(ks2+INDEX),z(ks2+INDEX+1)
   ibad = 1
 END DO
 IF (ibad == 1) GO TO 2300
 692 FORMAT (a23,' 6608, THE REQUEST FOR BOUNDARY SET ',i8,  &
     ' SUBSTRUCTURE ',2A4,' WAS NOT DEFINED.')
 
!     VI. PROCESS THE EQSS FROM THE SOF FOR EACH COMPONENT SUBSTRUCTURE.
!     ******************************************************************
 
 CALL OPEN (*2001,scr1,z(buf3),0)
 CALL OPEN (*2001,scr2,z(buf2),1)
 CALL sfetch (namold,nheqss,1,itest)
 ngrp = 1
 CALL sjump (ngrp)
 
!     READ AND PROCESS EQSS
 
 bad = .false.
 DO  i = 1,ncsub
   ii = 2*(i-1)
   CALL suread (z(score),-1,nout,itest)
   IF (andf(rshift(prtopt,3),1) /= 1) GO TO 2499
   CALL cmiwrt (1,namold,z(ks1+ii),score,nout,z,z)
   2499 CONTINUE
   
!     FIND A BOUNDARY SET FOR THE COMPONENT
   
   inxt = 1
   fset = .false.
   737 DO  j = inxt,nwbs,3
     jj = j - 1
     IF (z(ks2+jj) == z(ks1+ii) .AND. z(ks2+jj+1) == z(ks1+ii+1)) GO TO 704
   END DO
   IF (fset) GO TO 735
   
!     NO BOUNDARY SET FOR COMPONENT - IMPLIES ENTIRE SUBSTRUCTURE WILL B
!     REDUCED - POSSSIBLE ERROR.
   
   IF (nout /= 0) WRITE(outt,791) uim,z(ks1+ii),z(ks1+ii+1),  &
       (namold(j),j=1,2)
   791 FORMAT (a29,' 6609, NO BOUNDARY SET HAS BEEN SPECIFIED FOR ',  &
       'COMPONENT ',2A4,' OF PSEUDOSTRUCTURE ',2A4, /35X,  &
       'ALL DEGREES OF FREEDOM WILL BE REDUCED.')
   CALL WRITE (scr2,aray(1),0,1)
   CYCLE
   
!     COMPONENT HAS A BOUNDARY SET, CALL EQSCOD TO ACCOUNT FOR POSSIBLE
!     MULTIPLE IP NUMBERS.
   
   704 IF (fset) GO TO 736
   CALL eqscod (score,nout,z)
   
!     DEFINE ARRAY TO CB - DEGREES OF FREEDOM RETAINED AS BOUNDARY POINT
   
   ist  = score + nout
   ifin = ist + nout/3 - 1
   DO  j = ist,ifin
     z(j) = 0
   END DO
   
!     LOCATE BOUNDARY SET ON SCR1
   
   736 inxt = jj + 4
   fset = .true.
   nset = z(ks2+jj+2)
   DO  j = 1,nsid
     IF (nset == isid(j)) EXIT
   END DO
   766 nrec = j - 1
   IF (nrec == 0) GO TO 716
   DO  jj = 1,nrec
     CALL fwdrec (*2002,scr1)
   END DO
   
!     READ BOUNDARY DATA AND UPDATE CB
   
   716 CALL READ (*2002,*730,scr1,aray,2,0,nnn)
   
!     LOCATE GRID ID IN EQSS AND SETS OF VALUES IF THE GRID IS MULTIPLY
   
   IF (nout == 0) GO TO 717
   CALL gridip (aray(1),score,nout,ipset,cset,no,z,loc)
   IF (iiierr /= 1) GO TO 718
   717 bad = .true.
   WRITE  (outt,714) ufm,aray(1),nset,z(ks1+ii),z(ks1+ii+1)
   714 FORMAT (a23,' 6611, GRID POINT',i9,' SPECIFIED IN BOUNDARY SET',  &
       i9,' FOR SUBSTRUCTURE ',2A4,' DOES NOT EXIST.')
   718 iadd = loc
   IF (no > 1) GO TO 710
   icomp = z(iadd+2) - lshift(rshift(z(iadd+2),26),26)
   GO TO 711
   710 icomp = 0
   DO  j = 1,no
     cset(j) = cset(j) - lshift(rshift(cset(j),26),26)
     icomp = orf(icomp,cset(j))
   END DO
   
!     CHECK THAT THE RETAINED DOF ARE A SUBSET OF THE ORIGINAL.
   
   711 IF (andf( aray(2),icomp ) == aray(2).OR.iiierr == 1) GO TO 715
   WRITE  (outt,792) uwm,aray(1),z(ks1+ii),z(ks1+ii+1)
   792 FORMAT (a25,' 6610, DEGREES OF FREEDOM AT GRID POINT',i9,  &
       ' COMPONENT SUBSTRUCTURE ',2A4, /31X,'INCLUDED IN A ',  &
   'BOUNDARY SET DO NOT EXIST. REQUEST WILL BE IGNORED.')
     aray(2) = aray(2) - (orf(aray(2),icomp)-icomp)
     
!     UPDATE CB ARRAY
     
     715 IF (no > 1) GO TO 757
     nent = (iadd-score)/3
     z(ist+nent) = orf(z(ist+nent),aray(2))
     GO TO 716
     757 nent = (iadd-score)/3
     DO  j = 1,no
       z(ist+nent+j-1) = orf(z(ist+nent+j-1),aray(2))
     END DO
     GO TO 716
     
!     BOUNDARY SET COMPLETE, IS THERE ANOTHER
     
     730 CALL REWIND (scr1)
     GO TO 737
     
!     WRITE IP AND CB ON SCR2
     
     735 i1 = score
     i2 = i1 + nout - 1
     ii = -1
     DO  j = i1,i2,3
       ii = ii + 1
       aray(1) = andf(z(j+2),z(ist+ii))
       IF (aray(1) /= 0) CALL WRITE (scr2,z(j+1),1,0)
       IF (aray(1) /= 0) CALL WRITE (scr2,aray(1),1,0)
     END DO
     CALL WRITE (scr2,aray(1),0,1)
   END DO
   CALL CLOSE (scr1,1)
   CALL CLOSE (scr2,1)
   IF (bad) GO TO 2300
   
!     VII. PROCESS MASTER SIL LIST AND ALLOCATE SPACE FOR CNEW
!     ********************************************************
   
   j = 0
   800 CALL suread (z(score+j),2,nout,itest)
   IF (itest == 3) GO TO 810
   j = j + 3
   GO TO 800
   810 nw = j - 3
   DO  i = 1,nw,3
     jj = i - 1
     z(score+jj+2) = 0
   END DO
   CALL OPEN (*2001,scr2,z(buf2),0)
   840 CALL READ (*860,*850,scr2,aray,2,0,nnn)
   iloc = 3*aray(1) - 3
   z( score+iloc+2 ) = orf(z(score+iloc+2),aray(2))
   GO TO 840
   
!     READ NEXT COMPONENT
   
   850 GO TO 840
   
!     PROCESSING COMPLETE
   
   860 CALL CLOSE (scr2,1)
   ks3   = score
   score = score + nw
   kf3   = score - 1
   
!     VIII. DEFINE PARTITIONING VECTORS PVX AND USX
!     *********************************************
   
   CALL gopen (pvx,z(buf2),1)
   
!     GENERATE PVX DATA BLOCK IN CORE
   
   jjj = 0
   DO  i = 1,nw,3
     icode = z(ks3+i)
     CALL decode (icode,listo,nrow)
     DO  j = 1,nrow
       rz(score+jjj+j-1) = 0.0
     END DO
     icode = z(ks3+i+1)
     CALL decode (icode,listn,nnew)
     DO  j = 1,nrow
       listo(j) = listo(j) + 1
     END DO
     IF (nnew == 0) GO TO 960
     DO  j = 1,nnew
       listn(j) = listn(j) + 1
     END DO
     
!     FIND DOF THAT REMAIN AT GIVEN IP
     
     DO  j  = 1,nnew
       DO  jj = 1,nrow
         IF (listn(j) == listo(jj)) GO TO 943
       END DO
       CYCLE
       943 ijk(j) = jj
     END DO
     DO  j = 1,nnew
       ik = ijk(j)
       rz(score+jjj+ik-1) = 1.0
     END DO
     960 jjj = jjj + nrow
   END DO
   
!     SET PARAMETERS AND CALL PACK
   
   mcb(1) = pvx
   mcb(2) = 0
   mcb(3) = jjj
   mcb(4) = 2
   mcb(5) = 1
   mcb(6) = 0
   mcb(7) = 0
   typin  = 1
   typout = 1
   incr   = 1
   irow   = 1
   nrow   = jjj
   CALL pack (rz(score),pvx,mcb)
   CALL wrttrl (mcb)
   CALL CLOSE (pvx,1)
   IF (lonly) GO TO 1070
   
!     PROCESS USX USET EQUIVALENT
   
   CALL OPEN  (*2001,usx,z(buf2),1)
   CALL fname (usx,aray )
   CALL WRITE (usx,aray,2,0)
   CALL WRITE (usx,0.0 ,1,0)
   CALL WRITE (usx,0.0 ,1,1)
   mcb(1) = usx
   mcb(2) = 0
   mcb(3) = jjj
   mcb(4) = 0
   mcb(5) = iba + ibo + ibf
   mcb(6) = 0
   mcb(7) = 0
   DO  j = 1,jjj
     jj = j - 1
!WKBDB 8/94 ALPHA-VMS
!      IF (RZ(SCORE+JJ) .EQ. 0.0) Z(SCORE+JJ) = IBF + IBO
!      IF (RZ(SCORE+JJ) .EQ. 1.0) Z(SCORE+JJ) = IBF + IBA
!WKBDE 8/94 ALPHA-VMS
!WKBNB 8/94 ALPHA-VMS
     IF (rz(score+jj) /= 0.0) GO TO 976
     z(score+jj) = ibf + ibo
     GO TO 977
     976   IF (rz(score+jj) == 1.0) z(score+jj) = ibf + iba
     977   CONTINUE
!WKBNE 8/94 ALPHA-VMS
   END DO
   CALL WRITE (usx,z(score),jjj,1)
   CALL wrttrl (mcb)
   CALL CLOSE (usx,1)
   
!     IX. PROCESS THE SOF FOR THE REDUCED STRUCTURE
!     *********************************************
   
   
!     PROCESS THE EQSS FOR EACH COMPONENT SUBSTRUCTURE
   
   CALL OPEN (*2001,scr1,z(buf1),1)
   CALL sfetch (namold,nheqss,1,itest)
   
!     UPDATE (SIL,C) REPLACING SIL WITH IPNEW
   
   ipnew = 1
   DO  i = ks3,kf3,3
     IF (z(i+2) == 0.0) THEN
       GO TO  1004
     ELSE
       GO TO  1003
     END IF
     1004 z(i) = 0
     CYCLE
     1003 z(i)  = ipnew
     ipnew = ipnew + 1
   END DO
   nipnew = ipnew - 1
   ngrp   = 1
   CALL sjump (ngrp)
   DO  j = 1,ncsub
     CALL suread (z(score),-1,nout,itest)
     
!     WRITE EQSS ENTRY ON SCR1 IF THE OLD IP NUMBER STILL EXISTS IN THE
!     REDUCED STRUCTURE, ALSO UPDATE DOF CODE.
     
     IF (nout == 0) GO TO 1015
     DO  i = 1,nout,3
       ii  = i - 1
       ipo = z(score+ii+1)
       iadd= ks3 + (ipo-1)*3
       IF (z(iadd) == 0) CYCLE
       aray(1) = z(score+ii)
       aray(2) = z(iadd  )
       aray(3) = z(iadd+2)
       CALL WRITE (scr1,aray,3,0)
     END DO
     1015 CALL WRITE (scr1,0,0,1)
   END DO
   
!     GENERATE NEW MASTER (SIL,C) LIST
   
   isil = 1
   DO  i = ks3,kf3,3
     IF (z(i) == 0) CYCLE
     icode = z(i+2)
     CALL decode (icode,listn,ndof)
     aray(1) = isil
     aray(2) = z(i+2)
     CALL WRITE (scr1,aray,2,0)
     isil = isil + ndof
   END DO
   CALL WRITE (scr1,aray,0,1)
   CALL CLOSE (scr1,1)
   IF (dry == 0) GO TO 8612
   
!     WRITE FIRST GROUP OF EQSS
   
   CALL OPEN (*2001,scr1,z(buf1),0)
   CALL setlvl (namnew,1,namold,itest,28)
   IF (itest == 8) GO TO 6518
   itest = 3
   CALL sfetch (namnew,nheqss,2,itest)
   itest = 1
   CALL suwrt (namnew,2,itest)
   itest = 1
   CALL suwrt (ncsub,1,itest)
   itest = 1
   CALL suwrt (nipnew,1,itest)
   DO  i = ks1,kf1,2
     itest = 1
     CALL suwrt (z(i),2,itest)
   END DO
   itest = 2
   CALL suwrt (z(i),0,itest)
   1043 CALL READ (*1041,*1042,scr1,z(score),nz,0,nnn)
   GO TO 2004
   1042 CALL suwrt (z(score),nnn,2)
   GO TO 1043
   1041 itest = 3
   CALL suwrt (aray,0,itest)
   CALL CLOSE (scr1,1)
   
!     WRITE BGSS FILE
   
   CALL sfetch (namold,nhbgss,1,itest)
   ngrp = 1
   CALL sjump (ngrp)
   CALL suread (z(score),-1,nout,itest)
   j = 0
   
!     THE CID S THAT BELONG TO POINTS THAT ARE COMPLETELY REDUCED
!     WILL BE ACCUMULATED IN BUF3.
   
   jjj1 = 2
   DO  i = 1,nout,4
     ii = i - 1
     IF (z(ks3+jjj1) == 0.0) THEN
       GO TO  1051
     END IF
     1052 IF (z(score+ii) == 0) GO TO 1053
     z(buf3+j) = z(score+ii)
     j = j + 1
     GO TO 1053
     1051 z(score+ii) = -1*tpow(2)
     1053 jjj1 = jjj1 + 3
   END DO
   ncsred = j
   itest = 3
   CALL sfetch (namnew,nhbgss,2,itest)
   itest = 1
   CALL suwrt (namnew,2,itest)
   itest = 2
   CALL suwrt (nipnew,1,itest)
   DO  i = 1,nout,4
     ii = i - 1
     IF (z(score+ii) == -tpow(2)) CYCLE
     itest = 1
     CALL suwrt (z(score+ii),4,itest)
   END DO
   itest = 2
   CALL suwrt (aray,0,itest)
   itest = 3
   CALL suwrt (aray,0,itest)
   
!     PROCESS THE CSTM FILES
   
   IF (ncsred /= 0) GO TO 1063
   CALL sfetch (namold,nhcstm,1,itest)
   IF (itest == 3) GO TO 1070
   CALL suread (z(score),-2,nout,itest)
   z(score  ) = namnew(1)
   z(score+1) = namnew(2)
   itest = 3
   CALL sfetch (namnew,nhcstm,2,itest)
   itest = 3
   CALL suwrt (z(score),nout,itest)
   GO TO 1070
   1063 CALL sfetch (namold,nhcstm,1,itest)
   IF (itest == 3) GO TO 1070
   ngrp = 1
   CALL sjump (ngrp)
   
!     SORT THE DELETED CID S
   
   CALL sort (0,0,1,1,z(buf3),ncsred)
   
!     READ ALL RETAINED CSTM DATA INTO OPEN CORE
   
   j = 0
   1065 CALL suread (z(score+j),14,nout,itest)
   IF (itest == 2) GO TO 1066
   IF (z(score+j) == 0) GO TO 1065
   kid = z(score+j)
   CALL bisloc (*1065,kid,z(buf3),1,ncsred,jp)
   j = j + 14
   GO TO 1065
   1066 itest = 3
   CALL sfetch (namnew,nhcstm,2,itest)
   itest = 2
   CALL suwrt (namnew,2,itest)
   itest = 2
   CALL suwrt (z(score),j,itest)
   itest = 3
   CALL suwrt (aray,0,itest)
   1070 CONTINUE
   
!     PROCESS LODS ITEM
   
   CALL sfetch (namold,litm,1,itest)
   IF (itest == 3) GO TO 1080
   CALL suread (z(score),-2,nout,itest)
   z(score  ) = namnew(1)
   z(score+1) = namnew(2)
   itest = 3
   CALL sfetch (namnew,litm,2,itest)
   itest = 3
   CALL suwrt (z(score),nout,itest)
   1080 CONTINUE
   IF (lonly) GO TO 8511
   
!     PROCESS PLTS ITEM
   
   CALL sfetch (namold,nhplts,1,itest)
   IF (itest == 3) GO TO 1090
   CALL suread (z(score),-1,nout,itest)
   z(score  ) = namnew(1)
   z(score+1) = namnew(2)
   itest = 3
   CALL sfetch (namnew,nhplts,2,itest)
   itest = 2
   CALL suwrt (z(score),nout,itest)
   itest = 3
   CALL suwrt (z(score),0,itest)
   1090 CONTINUE
   
!     PROCESS OUTPUT REQUESTS
   
   IF (andf(rshift(prtopt,4),1) /= 1) GO TO 8211
   
!     WRITE EQSS FOR NEW STRUCTURE
   
   CALL sfetch (namnew,nheqss,1,itest)
   CALL suread (z(score),4,nout,itest)
   CALL suread (z(score),-1,nout,itest)
   ist = score + nout
   DO  i = 1,ncsub
     CALL suread (z(ist),-1,nout,itest)
     iadd = score + 2*(i-1)
     CALL cmiwrt (1,namnew,z(iadd),ist,nout,z,z)
   END DO
   CALL suread (z(ist),-1,nout,itest)
   CALL cmiwrt (8,namnew,0,ist,nout,z,z)
   8211 IF (andf(rshift(prtopt,5),1) /= 1) GO TO 8311
   
!     WRITE NEW BGSS
   
   CALL sfetch (namnew,nhbgss,1,itest)
   ngrp = 1
   CALL sjump (ngrp)
   ist = score
   CALL suread (z(ist),-1,nout,itest)
   CALL cmiwrt (2,namnew,namnew,ist,nout,z,z)
   8311 IF (andf(rshift(prtopt,6),1) /= 1) GO TO 8411
   
!     WRITE CSTM ITEM
   
   CALL sfetch (namnew,nhcstm,1,itest)
   IF (itest == 3) GO TO 8411
   ngrp = 1
   CALL sjump (ngrp)
   ist = score
   CALL suread (z(ist),-1,nout,itest)
   CALL cmiwrt (3,namnew,namnew,ist,nout,z,z)
   8411 IF (andf(rshift(prtopt,7),1) /= 1) GO TO 8511
   
!     WRITE PLTS ITEM
   
   CALL sfetch (namnew,nhplts,1,itest)
   IF (itest == 3) GO TO 8511
   ist = score
   CALL suread (z(ist),3,nout,itest)
   CALL suread (z(ist),-1,nout,itest)
   CALL cmiwrt (4,namnew,namnew,ist,nout,z,z)
   8511 IF (andf(rshift(prtopt,8),1) /= 1) GO TO 8611
   
!     WRITE LODS ITEM
   
   CALL sfetch (namnew,lods,1,itest)
   IF (itest == 3) GO TO 8611
   CALL suread (z(score),4,nout,itest)
   CALL suread (z(score),-1,nout,itest)
   ist   = score + nout
   itype = 5
   IF (litm == loap) itype = 7
   DO  i = 1,ncsub
     iadd = score+2*(i-1)
     CALL suread (z(ist),-1,nout,itest)
     CALL cmiwrt (itype,namnew,z(iadd),ist,nout,z,z)
     itype = 6
   END DO
   8611 CONTINUE
   IF (lonly) GO TO 1105
   
!     X. GENERATE THE INX OUTPUT DATA BLOCK
!     *************************************
   
   8612 CALL gopen (inx,z(buf2),1)
   mcb(1) = inx
   mcb(2) = 0
   mcb(3) = isil - 1
   mcb(4) = 1
   mcb(5) = 1
   mcb(6) = 0
   mcb(7) = 0
   typin  = 1
   typout = 1
   incr   = 1
   isilm1 = isil - 1
   DO  i = 1,isilm1
     irow = i
     nrow = i
     CALL pack (1.0,inx,mcb)
   END DO
   CALL wrttrl (mcb)
   CALL CLOSE (inx,1)
   1105 CALL sofcls
   RETURN
   
   2100 WRITE  (outt,2101) ufm
   2101 FORMAT (a23,' 6535, MODULE REDUCE TERMINATING DUE TO ABOVE ',  &
       'SUBSTRUCTURE CONTROL ERRORS.')
   GO TO 2400
   
   2200 WRITE  (outt,2201) ufm
   2201 FORMAT (a23,' 6536, MODULE REDUCE TERMINATING DUE TO ABOVE ',  &
       'ERRORS IN BULK DATA.')
   CALL CLOSE (geom4,1)
   GO TO 2400
   
   2300 WRITE  (outt,2301) ufm
   2301 FORMAT (a23,' 6537, MODULE REDUCE TERMINATING DUE TO ABOVE ',  &
       'ERRORS.')
   2400 dry = -2
   CALL sofcls
   RETURN
   
   6518 WRITE  (outt,6519) ufm
   6519 FORMAT (a23,' 6518, ONE OF THE COMPONENT SUBSTRUCTURES HAS BEEN ',  &
       'USED IN A PREVIOUS COMBINE OR REDUCE.')
   GO TO 2300
   2001 imsg = -1
   GO TO 2998
   2002 imsg = -2
   GO TO 2998
   2003 imsg = -3
   GO TO 2998
   2004 imsg = -8
   2998 CALL mesage (imsg,ifile,modnam)
   RETURN
 END SUBROUTINE reduce
