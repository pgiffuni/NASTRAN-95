SUBROUTINE apd
     
 EXTERNAL        andf,orf
 LOGICAL :: lmkaer ,lskip,lset,lsplin
 INTEGER :: iz(1),flag,silgp
 INTEGER :: eid,pid,cp,cidbx,acsid,silb,scr1,scr2,scr3,scr4,  &
     scr5,ecta,bgpa,gpla,useta,sila,cstma,acpt,buf10, buf11,buf12
 INTEGER :: caero2(3),caero3(3),caero4(3),caero5(3)
 INTEGER :: paero2(3),paero3(3),paero4(3),paero5(3)
 INTEGER :: kspl(3),andf,FILE
 INTEGER :: splin3(3)
 INTEGER :: caero1(3),paero1(3),aero(3),splin1(3),splin2(3),  &
     set1(3),set2(3),mkaer2(3),mkaer1(3),fluttr(3),  &
     aefact(3),flfact(3),buf(7),msg(7),  &
     buf1,buf2,buf3,buf4,buf5,buf6,buf7,buf8,buf9,  &
     edt,eqdyn,ect,bgpdt,sild,usetd,cstm,  &
     eqaero,spline,aeror,flist,gpld,nam(2),  &
     aerx(3),symxz,symxy,corwds,rdrew,clsrew, orf,sysbuf,wtrew,pspa,nbca(3)
 INTEGER :: ca2s,ca2e,ca3s,ca3e,ca4s,ca4e,ca5s,ca5e
 INTEGER :: pa2s,pa2e,pa3s,pa3e,pa4s,pa4e,pa5s,pa5e
 INTEGER :: msg1(9),msg2(5),msg3(6),msg4(10)
 COMMON /BLANK / nk,nj,luseta,bov
 COMMON /system/ sysbuf,NOT
 COMMON /apd1c / eid,pid,cp,nspan,nchord,lspan,lchord,igid,  &
     x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
     ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,  &
     scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,  &
     sila,cstma,acpt,buf10,buf11,buf12,next,left,isiln,  &
     ncam,naef1,naef2,nca1,nca2,ca2s,ca2e,ca3s,ca3e,  &
     ca4s,ca4e,npa1,npa2,pa2s,pa2e,pa3s,pa3e,pa4s,pa4e, ca5s,ca5e,pa5s,pa5e
 COMMON /two   / itwo(32)
 COMMON /bitpos/ ibit(64)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1), iz(1))
 EQUIVALENCE     (aerx(1),symxz),(aerx(2),symxy),(aerx(3),bref)
 DATA    rdrew , clsrew,wtrew /0,1,1 /
 DATA    msg1  / 4HSETI,4H AND,4H/OR ,4HSPLI,4HNEI ,4HCARD,4HS re,  &
     4HQUIR,4HED  /
 DATA    msg2  / 4HNO a,4HERO ,4HCARD,4H fou,4HND    /
 DATA    msg3  / 4HNO c,4HAERO,4H  ca,4HARDS,4HFOUN,4HD    /
 DATA    msg4  / 4HNEIT,4HHER ,4HMKAE,4HRO1 ,4HOR  ,4HMKAE,4HRO2 ,  &
     4HCARD,4HS fo,4HUND /
 DATA    caero2/ 4301,43, 0/ , caero3 /4401,44, 0 /
 DATA    caero4/ 4501,45, 0/ , paero2 /4601,46, 0 /
 DATA    paero3/ 4701,47, 0/ , paero4 /4801,48, 0 /
 DATA    caero5/ 5001,50, 0/ , paero5 /5101,51, 0 /
 DATA    splin3/ 4901,49, 0/
 DATA    kspl  / 200 , 2, 0/
 DATA    caero1/ 3002,30,16 /, paero1 /3102,31,0  /,  &
     aero  / 3202,32,0  /, splin1 /3302,33,0  /,  &
     splin2/ 3402,34,0  /, set1   /3502,35,0  /,  &
     set2  / 3602,36,0  /, mkaer2 /3702,37,0  /,  &
     mkaer1/ 3802,38,0  /, fluttr /3902,39,0  /,  &
     aefact/ 4002,40,0  /, flfact /4102,41,0  /, nbca  / 3002,46,0/
 DATA    edt   , eqdyn,  ect, bgpdt, sild, usetd, cstm, gpld /  &
     101   , 102,    103, 104,   105,  106,   107,  108  /
 DATA    eqaero, spline, aeror, flist  / 201   , 206,    207,   209    /       
 DATA    msg   / 7*0 /,  nam / 4HAPD ,4H    /
 
 lca  = caero1(3)
 nogo = 0
 buf1 = korsz(iz) - sysbuf
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 buf5 = buf4 - sysbuf
 buf6 = buf5 - sysbuf
 buf7 = buf6 - sysbuf
 buf8 = buf7 - sysbuf
 buf9 = buf8 - sysbuf
 buf10= buf9 - sysbuf
 buf11= buf10- sysbuf
 buf12= buf11- sysbuf
 scr1 = 301
 scr2 = 302
 scr3 = 303
 scr4 = 304
 scr5 = 305
 ecta = 202
 bgpa = 203
 sila = 204
 useta= 205
 acpt = 208
 cstma= 210
 gpla = 211
 silgp= 212
 last = buf12 - 1
 IF (last <= 0) GO TO 995
 nj   = 0
 nk   = 0
 i17  = ibit(17)
 i20  = ibit(20)
 pspa = orf(itwo(i17),itwo(i20))
 
!     READ AERO CARDS
 
 left = last
 FILE = edt
 CALL preloc (*940,z(buf1),edt)
 CALL locate (*800,z(buf1),aero,flag)
 CALL READ (*960,*970,edt,z(1),6,1,flag)
 acsid= iz(1)
 izx  = 2
 vsound = z(izx)
 izx  = 3
 bref = z(izx)
 bov  = 0.0
 IF (vsound /= 0.0) bov = bref/(2.0*vsound)
 izx  = 5
 symxz= iz(izx)
 izx  = 6
 symxy= iz(izx)
 
!     READ AEFACT CARDS
 
 naef2 = 0
 CALL apdr (edt,z,left,naef1,naef2,flag,buf1,aefact)
 
!     READ CSTM TABLE
 
 FILE  = cstm
 ncst2 = naef2
 ncst1 = 0
 mcstm = 0
 buf(1)= cstm
 CALL rdtrl(buf)
 IF (buf(1) /= cstm) GO TO 100
 CALL gopen (cstm,z(buf2),rdrew)
 ncst1 = ncst2 + 1
 CALL READ (*980,*80,cstm,z(ncst1),left,0,ncst2)
 GO TO 970
 80 CALL CLOSE (cstm,clsrew)
 left  = left  - ncst2
 ncst2 = ncst1 + ncst2 - 1
 
!     FIND LARGEST CID OF CSTM
 
 DO  j = ncst1,ncst2,14
   IF (iz(j) < mcstm) CYCLE
   mcstm = iz(j)
 END DO
 
!     FIND AC TRANS
 
 100 IF (acsid == 0) GO TO 120
 IF (ncst1 == 0) GO TO 880
 DO  iacs = ncst1,ncst2,14
   IF (iz(iacs) == acsid) GO TO 120
 END DO
 GO TO 880
 
!     WRITE CSTM TO CSTMA
 
 120 CALL gopen (cstma,z(buf2),wtrew)
 IF (mcstm /= 0) CALL WRITE (cstma,z(ncst1),ncst2-ncst1+1,0)
 ncsa = mcstm
 
!     READ EQDYN INTO CORE
 
 next = ncst2 + 1
 FILE = eqdyn
 CALL gopen (eqdyn,z(buf3),rdrew)
 CALL skprec (eqdyn,1)
 CALL READ (*980,*140,eqdyn,z(next),left,0,nx)
 GO TO 970
 140 CALL CLOSE (eqdyn,clsrew)
 buf(1) = eqdyn
 CALL rdtrl (buf)
 nextra = buf(3)
 
!     CIDBX = LARGEST ID
!     NCRD  = NUMBER OF GRID AND SCALAR POINTS
 
 ncrd  = buf(2) - nextra
 ncrdo = ncrd
 cidbx = 1000000
 
!     WRITE SECOND RECORD OF EQDYN ONTO SCR1
 
 CALL gopen (scr1,z(buf3),wtrew)
 CALL WRITE (scr1,z(next),nx,0)
 
!     READ BGPDT
 
 FILE = bgpdt
 CALL gopen (bgpdt,z(buf4),rdrew)
 CALL READ (*980,*150,bgpdt,z(next),left,0,nx)
 GO TO 970
 150 CALL CLOSE (bgpdt,clsrew)
 
!     WRITE BGPDT TO BGPA
 
 CALL gopen (bgpa,z(buf4),wtrew)
 CALL WRITE (bgpa,z(next),nx,0)
 
!     READ USETD
 
 FILE = usetd
 CALL gopen (usetd,z(buf5),rdrew)
 CALL READ (*980,*160,usetd,z(next),left,0,nx)
 GO TO 970
 160 CALL CLOSE (usetd,clsrew)
 
!     MASK IN PS AND PA BITS
 
 n2 = next + nx - 1
 DO  i = next,n2
   iz(i) = orf(iz(i),pspa)
 END DO
 
!     WRITE USETD TO USETA
 
 FILE = useta
 CALL gopen (useta,z(buf5),wtrew)
 CALL WRITE (useta,z(next),nx,0)
 
!     READ ECT
 
 FILE   = ect
 buf(1) = ect
 CALL gopen (ecta,z(buf6),wtrew)
 CALL rdtrl (buf)
 IF (buf(1) /= ect) GO TO 210
 CALL gopen (ect,z(buf7),rdrew)
 180 CALL READ (*200,*190,ect,z(next),left,0,nx)
 GO TO 970
 190 CALL WRITE (ecta,z(next),nx,1)
 GO TO 180
 200 CALL CLOSE (ect,clsrew)
 210 CALL WRITE (ecta,nbca,3,0)
 
!     READ FIRST RECORD OF SILD INTO CORE
 
 FILE = sild
 CALL gopen (sild,z(buf8),rdrew)
 CALL READ (*980,*220,sild,z(next),left,0,nx)
 GO TO 970
 
!     WRITE FIRST RECORD OF SILD ONTO SILA
 
 220 buf(1) = sild
 CALL rdtrl (buf)
 
!     SILB  + 6 = NEXT DOF IN PROBLEM
!     ISILN + 6 = NEXT DOF WITHOUT EXTRA POINTS
 
 silb  = buf(2) - 5
 isiln = silb - nextra
 CALL gopen (sila,z(buf7),wtrew)
 CALL WRITE (sila,z(next),nx,0)
 
!     READ SECOND RECORD OF SILD INTO CORE
 
 CALL READ (*980,*230,sild,z(next),left,0,nx)
 GO TO 970
 230 CALL CLOSE (sild,clsrew)
 
!     WRITE SECOND RECORD OF SILD ONTO SCR2
 
 CALL gopen (scr2,z(buf8),wtrew)
 CALL WRITE (scr2,z(next),nx,0)
 
!     COPY GPLD TO GPLA
 
 FILE = gpld
 CALL gopen (gpld,z(buf9),rdrew)
 CALL READ (*980,*235,gpld,z(next),left,0,nx)
 GO TO 970
 235 CALL CLOSE (gpld,clsrew)
 CALL gopen (gpla,z(buf9),wtrew)
 CALL WRITE (gpla,z(next),nx,0)
 
!     READ CAERO CARDS INTO CORE
 
 nca2 = ncst2
 lcas = nca2 + 1
 CALL apdr (edt,z,left,nca1,nca2,flag,buf1,caero1)
 ca2e = nca2
 CALL apdr (edt,z,left,ca2s,ca2e,flag,buf1,caero2)
 ca3e = ca2e
 CALL apdr (edt,z,left,ca3s,ca3e,flag,buf1,caero3)
 ca4e = ca3e
 CALL apdr (edt,z,left,ca4s,ca4e,flag,buf1,caero4)
 ca5e = ca4e
 CALL apdr (edt,z,left,ca5s,ca5e,flag,buf1,caero5)
 lcae = MAX0(nca2,ca2e,ca3e,ca4e,ca5e)
 
!     READ PAERO CARDS INTO CORE
 
 npa2 = ca5e
 CALL apdr (edt,z,left,npa1,npa2,flag,buf1,paero1)
 pa2e = npa2
 CALL apdr (edt,z,left,pa2s,pa2e,flag,buf1,paero2)
 pa3e = pa2e
 CALL apdr (edt,z,left,pa3s,pa3e,flag,buf1,paero3)
 pa4e = pa3e
 CALL apdr (edt,z,left,pa4s,pa4e,flag,buf1,paero4)
 pa5e = pa4e
 CALL apdr (edt,z,left,pa5s,pa5e,flag,buf1,paero5)
 next = pa5e + 1
 CALL CLOSE (edt,clsrew)
 IF (nca1 == 0 .AND. ca2s == 0 .AND. ca3s == 0 .AND. ca4s == 0  &
     .AND. ca5s == 0) GO TO 820
 
!     OPEN ACPT
 
 CALL gopen (acpt,z(buf1),wtrew)
 
!     CALL CAERO TYPE
 
 IF (nca1 /= 0 .OR. ca2s /= 0) CALL apd12
 IF (ca3s /= 0) CALL apd3
 IF (ca4s /= 0) CALL apd4
 IF (ca5s /= 0) CALL apd5
 luseta = luseta + 5
 CALL WRITE (cstma,0,0,1)
 CALL CLOSE (cstma,clsrew)
 CALL CLOSE (acpt,clsrew)
 CALL WRITE (ecta,0,0,1)
 CALL CLOSE (ecta,clsrew)
 CALL WRITE (bgpa,0,0,1)
 CALL CLOSE (bgpa,clsrew)
 CALL WRITE (gpla,0,0,1)
 CALL CLOSE (gpla,clsrew)
 CALL WRITE (useta,0,0,1)
 CALL CLOSE (useta,clsrew)
 CALL WRITE (sila,0,0,1)
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,clsrew)
 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr2,clsrew)
 
!     READ SECOND RECORD OF EQAERO TABLE OFF SCR1
 
 FILE = scr1
 CALL gopen (scr1,z(buf3),rdrew)
 i = next
 411 CONTINUE
 IF (i+2 > next+left) GO TO 970
 CALL READ (*980,*420,scr1,z(i),2,0,nx)
 iz(i+2) = 0
 i = i + 3
 GO TO 411
 420 CALL CLOSE (scr1,clsrew)
 nx = i - next
 
!     SORT TABLE ON SILD VALUE
 
 CALL sort (0,0,3,2,z(next),nx)
 ny = next + nx - 1
 
!     REPLACE THIRD ENTRIES WITH INTERNAL GRID ID WITH OUT EXTRA
 
 k = 0
 DO  i = next,ny,3
   IF (iz(i+1)-(iz(i+1)/10)*10 == 3) CYCLE
   k = k + 1
   iz(i+2) = k
 END DO
 
!     SORT EQAERO TABLE
 
 CALL sort (0,0,3,1,z(next),nx)
 
!     CHECK FOR DUPLICATE EXT ID
 
 n1 = next + 3
 DO  i = n1,ny,3
   IF (iz(i-3) /= iz(i)) CYCLE
   CALL emsg (0,2329,1,2,0)
   nogo = 1
   WRITE  (NOT,421) iz(i)
   421 FORMAT (10X,26HDUPLICATE EXTERNAL id no. ,i8,11H generated.)
 END DO
 
!     WRITE FIRST RECORD OF EQAERO TABLE
 
 CALL gopen (eqaero,z(buf3),wtrew)
 DO  i = next,ny,3
   buf(1) = iz(i  )
   buf(2) = iz(i+2)
   CALL WRITE (eqaero,buf,2,0)
 END DO
 CALL WRITE (eqaero,0,0,1)
 
!     WRITE SECOND RECORD OF EQAERO TABLE
 
 DO  i = next,ny,3
   CALL WRITE (eqaero,iz(i),2,0)
 END DO
 CALL WRITE (eqaero,0,0,1)
 CALL CLOSE (eqaero,clsrew)
 
!     PUT ON SPLINE A RECORD OF K POINTS WITH
!     EXTERNAL ID , BGPA POINTERS, AND K COLUMN NUMBER
 
 FILE = useta
 n1   = next + nx
 CALL gopen (useta,z(buf3),rdrew)
 CALL READ (*980,*442,useta,z(n1),left-nx,0,n2)
 GO TO 970
 442 CALL CLOSE (useta,clsrew)
 CALL gopen (spline,z(buf3),wtrew)
 CALL WRITE (spline,kspl,3,0)
 mask = ibit(19)
 mask = itwo(mask)
 ko   = 1
 n3   = (ncrdo+nextra)*3 + next
 DO  i = next,ny,3
   IF (MOD(iz(i+1),10) /= 1) CYCLE
   k  = 0
   n4 = iz(i+1)/10 - 2
   DO  j = 1,6
     IF (andf(iz(n1+n4+j),mask) /= 0) k = k + 1
   END DO
   IF (k == 0) CYCLE
   buf(1) = iz(i  )
   buf(2) = iz(i+2)
   buf(3) = ko
   CALL WRITE (spline,buf,3,0)
   ko = ko + k
 END DO
 CALL WRITE (spline,0,0,1)
 CALL CLOSE (spline,2)
 
!     READ SECOND RECORD OF SILA TABLE
 
 CALL gopen (scr2,z(buf8),rdrew)
 CALL READ  (*980,*450,scr2,z(next),left,0,nx)
 GO TO 970
 450 CALL CLOSE (scr2,clsrew)
 CALL WRITE (sila,z(next),nx,1)
 CALL CLOSE (sila,clsrew)
 
!     BUILD SILGP TABLE
 
 CALL gopen (silgp,z(buf8),w_trew)
 ny = next + nx - 1
 k  = 0
 DO  i = next,ny,2
   iz(next+k) = iz(i)
   k  = k + 1
 END DO
 CALL WRITE (silgp,iz(next),k,1)
 CALL CLOSE (silgp,clsrew)
 
!     WRITE RECORD
 
 CALL gopen (aeror,z(buf2),wtrew)
 CALL WRITE (aeror,aerx,3,1)
 
!     READ IN MKAERO1 CARDS
 
 FILE = edt
 CALL preloc (*940,z(buf1),edt)
 lmkaer = .false.
 CALL locate (*510,z(buf1),mkaer1,flag)
 CALL READ (*980,*460,edt,z(next),left,0,nx)
 GO TO 970
 460 n1 = next
 lmkaer = .true.
 470 n2 = n1 + 7
 loop490:  DO  i = n1,n2
   IF (iz(i) == -1) EXIT loop490
   buf(1) = iz(i)
   n3 = n2 + 1
   n4 = n3 + 7
   DO  j = n3,n4
     IF (iz(j) == -1) CYCLE loop490
     buf(2) = iz(j)
     CALL WRITE (aeror,buf,2,0)
   END DO
 END DO loop490
 500 IF (n4-next+1 >= nx) GO TO 510
 n1 = n1 + 16
 GO TO 470
 
!     READ IN MKAER2 CARDS
 
 510 CALL locate (*530,z(buf1),mkaer2,flag)
 CALL READ (*980,*520,edt,z(next),left,0,nx)
 GO TO 970
 520 CALL WRITE (aeror,z(next),nx,0)
 lmkaer =.true.
 530 CALL WRITE (aeror,0,0,1)
 CALL CLOSE (aeror,clsrew)
 IF (lmkaer) GO TO 540
 GO TO 870
 
!     PROCESS SET1 CARDS
 
 540 CALL OPEN (*940,spline,z(buf2),3)
 lset =.false.
 CALL locate (*610,z(buf1),set1,flag)
 lset =.true.
 CALL READ (*980,*550,edt,z(next),left,0,nx)
 GO TO 970
 550 n3 = next + nx
 CALL gopen (eqaero,z(buf3),rdrew)
 left = corwds(iz(n3),iz(last))
 FILE = eqaero
 CALL READ (*980,*560,eqaero,z(n3),left,0,n4)
 GO TO 970
 560 n1 = next
 FILE = edt
 n2 = n1 + nx - 1
 CALL CLOSE (eqaero,clsrew)
 left = corwds(iz(next),iz(last))
 
!     CONVERT SET1 TO INTERNAL COOR NO
 
 lskip = .true.
 DO  i = n1,n2
   IF (iz(i) == -1) GO TO 590
   IF (lskip) GO TO 580
   kid = iz(i)
   CALL bisloc (*930,kid,iz(n3),2,n4/2,jp)
   k = n3 + jp
   IF (iz(k) > ncrdo) GO TO 930
   iz(i) = iz(k)
   CYCLE
   580 lskip =.false.
   CYCLE
   930 CALL emsg (0,2330,1,2,0)
   nogo = 1
   WRITE  (NOT,931) iz(n1),iz(i)
   931 FORMAT (10X,24HSET1 OR splin3 card no. ,i8,28H references EXTERNAL  &
       id no. ,i8,22H which does NOT EXIST.)
   CYCLE
   590 lskip = .true.
 END DO
 
!     WRITE OUT SET1 CARD ON SPLINE
 
 CALL WRITE (spline,set1,3,0)
 CALL WRITE (spline,z(next),nx,1)
 
!     PROCESS SET2 CARDS
 
 610 CALL locate (*660,z(buf1),set2,flag)
 lset =.true.
 CALL WRITE (spline,set2,3,0)
 620 CALL READ (*980,*650,edt,z(next),8,0,nx)
 CALL WRITE (spline,z(next),10,0)
 nx = iz(next+1)
 DO  i = lcas,lcae,lca
   IF (iz(i) == nx) GO TO 640
 END DO
 GO TO 830
 640 CALL WRITE (spline,z(i),lca,0)
 GO TO 620
 650 CALL WRITE (spline,0,0,1)
 
!     PROCESS SPLINE1 CARDS
 
 660 lsplin =.false.
 CALL locate (*710,z(buf1),splin1,flag)
 lsplin =.true.
 CALL WRITE (spline,splin1,3,0)
 ASSIGN 670 TO iret
 670 CALL READ (*980,*700,edt,z(next),6,0,nx)
 GO TO 671
 
!     INTERNAL ROUTINE TO ATTACH CAERO DATA TO SPLINE
 
 671 CONTINUE
 CALL WRITE (spline,z(next),10,0)
 nx = iz(next+1)
 DO  i = lcas,lcae,lca
   IF (iz(i) == nx) GO TO 690
 END DO
 GO TO 810
 690 CALL WRITE (spline,z(i),lca,0)
 IF (iz(next+2) > iz(next+3)) GO TO 691
 j1 = iz(i+4)*iz(i+3) + iz(i) - 1
 IF (iz(next+2) < iz(i) .OR. iz(next+3) > j1) GO TO 691
 GO TO 693
 691 nogo = 1
 CALL emsg (0,2331,1,2,0)
 WRITE  (NOT,692) iz(next),iz(i)
 692 FORMAT (10X,30HBOX picked on spline card no. ,i8,32HNOT generated  &
     by caero card no. ,i8,1H.)
 693 GO TO iret, (670,720,7651)
 700 CALL WRITE (spline,0,0,1)
 
!     PROCESS SPLINE2 CARDS
 
 710 CALL locate (*760,z(buf1),splin2,flag)
 lsplin =.true.
 CALL WRITE (spline,splin2,3,0)
 ASSIGN 720 TO iret
 720 CALL READ (*980,*750,edt,z(next),10,0,nx)
 GO TO 671
 750 CALL WRITE (spline,0,0,1)
 
!     PROCESS  SPLINE3 CARDS
 
 760 nsplie = next - 1
 CALL apdr (edt,z,left,nsplis,nsplie,flag,buf1,splin3)
 IF (nsplis == 0) GO TO 769
 FILE = eqaero
 n3   = nsplie + 1
 CALL gopen (eqaero,z(buf3),rdrew)
 CALL READ (*980,*761,eqaero,z(n3),left,0,n4)
 GO TO 970
 761 FILE = edt
 CALL CLOSE (eqaero,clsrew)
 n4   = n4/2
 lset = .true.
 left = left + flag
 lsplin = .true.
 ASSIGN 7651 TO iret
 CALL WRITE (spline,splin3,3,0)
 isp = nsplis
 nls = 0
 
!     PICK UP NEXT SPLIN3 AND ATTACHED CAERO
 
 765 isp = isp + nls
 IF (isp >= nsplie) GO TO 768
 CALL apdoe (iz(isp),z,isp,nsplie,flag,nls)
 nls = nls + 1
 nx  = iz(isp+1)
 DO  i = lcas,lcae,lca
   j   = i
   IF (iz(i) == nx) GO TO 767
 END DO
 GO TO 810
 767 j1  = iz(i+3)*iz(i+4) + nx - 1
 iz(next) = iz(isp)
 IF (iz(isp+2) < iz(isp+1)) GO TO 691
 IF (iz(isp+2) > j1) GO TO 691
 
!     CONVERT TO INTERNAL ID
 
 n2 = nls - 4
 DO  i = 1,n2,3
   n1 = iz(isp+i+3)
   CALL bisloc (*7653,n1,iz(n3),2,n4,jp)
   IF (iz(n3+jp) > ncrdo) GO TO 7653
   iz(isp+i+3) = iz(n3+jp)
 END DO
 7651 CALL WRITE (spline,nls+lca,1,0)
 CALL WRITE (spline,iz(isp),nls,0)
 nls = nls + 1
 CALL WRITE (spline,iz(j),lca,0)
 GO TO 765
 7653 nogo = 1
 CALL emsg (0,2330,1,2,0)
 WRITE (NOT,931) iz(isp),n1
 GO TO 7651
 768 CALL WRITE (spline,0,0,1)
 769 CALL CLOSE (spline,clsrew)
 
!     CREATE FLIST TABLE
 
 CALL gopen (flist,z(buf2),wtrew)
 CALL locate (*800,z(buf1),aero,flag)
 CALL READ (*980,*770,edt,z(next),left,0,nx)
 GO TO 970
 770 CALL WRITE (flist,aero,3,0)
 CALL WRITE (flist,z(next),nx,1)
 CALL locate (*785,z(buf1),flfact,flag)
 CALL READ (*980,*780,edt,z(next),left,1,nx)
 GO TO 970
 780 CALL WRITE (flist,flfact,3,0)
 CALL WRITE (flist,z(next),nx,1)
 785 CALL locate (*900,z(buf1),fluttr,flag)
 CALL READ (*980,*790,edt,z(next),left,0,nx)
 GO TO 970
 790 CALL WRITE (flist,fluttr,3,0)
 CALL WRITE (flist,z(next),nx,1)
 900 CALL CLOSE (flist,clsrew)
 CALL CLOSE (edt,clsrew)
 msg(1) = aeror
 msg(2) = 1
 CALL wrttrl (msg)
 msg(1) = eqdyn
 CALL rdtrl (msg)
 msg(1) = eqaero
 msg(2) = ncrd + nextra
 CALL wrttrl (msg)
 msg(1) = bgpdt
 CALL rdtrl (msg(1))
 msg(3) = ncrd - msg(2)
 msg(1) = bgpa
 msg(2) = ncrd
 CALL wrttrl (msg)
 msg(1) = sila
 msg(2) = luseta
 msg(3) = nextra
 CALL wrttrl (msg)
 msg(1) = acpt
 msg(2) = 1
 CALL wrttrl (msg)
 msg(1) = gpla
 msg(2) = ncrd + nextra
 CALL wrttrl (msg)
 msg(1) = cstm
 CALL rdtrl (msg)
 IF (msg(1) < 0) msg(3) = 0
 msg(1) = cstma
 msg(3) = msg(3) + mcstm - ncsa
 CALL wrttrl (msg)
 msg(1) = useta
 msg(2) = luseta
 msg(3) = nextra
 msg(4) = pspa
 CALL wrttrl (msg)
 msg(1) = edt
 CALL rdtrl (msg)
 msg(1) = flist
 CALL wrttrl (msg)
 msg(1) = edt
 CALL rdtrl (msg)
 msg(1) = spline
 msg(2) = orf(msg(2),itwo(18))
 CALL wrttrl (msg)
 msg(1) = ect
 CALL rdtrl (msg)
 n1     = (nbca(2)-1)/16 + 2
 n2     = nbca(2) - (n1-2)*16 + 16
 msg(n1)= orf(msg(n1),itwo(n2))
 msg(1) = ecta
 CALL wrttrl (msg)
 
!     PUT OUT SILGP TRAILER
 
 msg(1) = silgp
 msg(2) = ncrd
 msg(3) = luseta - nextra
 msg(4) = 0
 msg(5) = 0
 msg(6) = 0
 msg(7) = 0
 CALL wrttrl (msg)
 IF (nogo == 1) CALL mesage (-37,0,nam)
 IF (lset .AND. lsplin) RETURN
 
!     ERROR MESSAGES
 
 CALL emsg (35,-2328,1,2,msg1)
 800 CALL emsg (18,-2318,1,3,msg2)
 810 CALL emsg (0,-2324,1,2,0)
 WRITE  (NOT,811) nx
 811 FORMAT (10X,19HCAERO  element no. ,i8,  &
     45H referenced on a splinei card does NOT EXIST.)
 812 CALL mesage (-61,0,nam)
 820 CALL emsg (21,-2319,1,2,msg3)
 830 CALL emsg (0,-2325,1,2,0)
 WRITE  (NOT,831) nx
 831 FORMAT (10X,19HCAERO  element no. ,i8,  &
     42H referenced on a set2 card does NOT EXIST.)
 GO TO 812
 870 CALL emsg (38,-2322,1,2,msg4)
 880 CALL mesage (-30,25,acsid)
 940 ip1 = -1
 950 CALL mesage (ip1,FILE,nam)
 995 ip1 = -8
 GO TO 950
 960 ip1 = -2
 GO TO 950
 970 ip1 = 3
 GO TO 950
 980 GO TO 960
END SUBROUTINE apd
