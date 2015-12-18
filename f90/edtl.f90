SUBROUTINE edtl (nedt,ilist,pg)
     
!     THIS SUBROUTINE COMPUTES THE ELEMENT TEMPERATURE AND ENFORCED
!     DEFORMATION LOADS
 
 
 INTEGER, INTENT(IN)                      :: nedt
 INTEGER, INTENT(IN)                      :: ilist(1)
 INTEGER, INTENT(IN OUT)                  :: pg(7)
 LOGICAL :: eorflg,endid,bufflg,record
 INTEGER :: pcomp(2),pcomp1(2),pcomp2(2), iparm(2),tlist(1080),outpt
 REAL :: core,ti
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /BLANK / nrowsp,iparam,comps
 COMMON /system/ ksystm(64)
 COMMON /packx / itya,ityb,ii,jj,incur
 COMMON /zzzzzz/ core(1)
 COMMON /xcstm / tgb(3,3)
 COMMON /tranx / idum1(14)
 COMMON /fpt   / TO,nsil,ngptt,nstart,lcore
 COMMON /loadx / lcare,n(3),cstm,sil,nnn,ecpt,mpt,gptt,edt,impt,  &
     igptt,iec,nn(3),dit,icm
 COMMON /trimex/ mecpt(200)
 COMMON /sgtmpd/ ti(33)
 COMMON /matin / matid,inflag,temp,stress,sinth,costh
 COMMON /matout/ e1,g,nu,rho,alpha,to1,GE,sigmat,sigmac,sigmas, SPACE(10)
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /ssgett/ eltype,oldel,eorflg,endid,bufflg,itemp,ideft, idefm,record
 COMMON /ssgwrk/ dum(300)
 COMMON /compst/ ipcmp,npcmp,ipcmp1,npcmp1,ipcmp2,npcmp2
 EQUIVALENCE     (ksystm( 1),sysbuf),(ksystm( 2),outpt ),  &
     (ksystm(55),iprec ),(ksystm(56),ithrml),  &
     (ti(7)     ,icheck),(ti(6)     ,iflag )
 DATA    iparm , ipgtt/ 4HEDTL,4H    ,4HGPTT   /
 DATA    crod  , ctube, conrod, cbar, pcomps   /  &
     1     , 3    , 10    , 34  , 112      /
 DATA    pcomp ,        pcomp1,       pcomp2   /  &
     5502  , 55,    5602, 56,     5702, 57 /
 
 igptt = ipgtt
 
!     CHECK IF HEAT FORMULATION
 
 IF (ithrml /= 0) RETURN
 
 itemp = 0
 ideft = nedt
 GO TO 10
 
 
 ENTRY templ (ntemp,ilist,pg)
!     ============================
 
 IF (ithrml /= 0) RETURN
 ideft = 0
 itemp = ntemp
 
!     START SEARCH POINTERS AT ZERO
 
 10 itya  = 1
 CALL delset
 ityb  = 1
 ipr   = iprec
 IF (ipr /= 1) ipr = 0
 ii    = 1
 jj    = nrowsp
 incur = 1
 nnn   = 0
 nogptt= 0
 idum1(1) = 0
 icm   = 1
 noedt = 0
 CALL delset
 lpcomp = 0
 
!     SET CORE SIZE AND BUFFERS
 
 lcore= korsz(core) - nrowsp
 buf1 = lcore - sysbuf - 2
 buf2 = buf1  - sysbuf - 2
 buf3 = buf2  - sysbuf - 2
 buf4 = buf3  - sysbuf - 2
 buf5 = buf4  - sysbuf - 2
 
!     OPEN FILES--
 
!     READ FILE PCOMPS INTO CORE ONLY IF PARAM COMPS = -1,
!     INDICATING THE PRESENCE OF LAMINATED COMPOSITE ELEMENTS
 
 IF (comps /= -1) GO TO 25
 
 ipm = pcomps
 CALL preloc (*750,core(buf2),pcomps)
 
 ipcmp  = nrowsp + 1
 ipcmp1 = ipcmp
 npcmp  = 0
 npcmp1 = 0
 npcmp2 = 0
 
 lcore = buf5 - nrowsp - 1
 
!     LOCATE PCOMP DATA AND READ INTO CORE
 
 CALL locate (*14,core(buf2),pcomp,flag)
 
 CALL READ (*860,*12,pcomps,core(ipcmp),lcore,0,npcmp)
 GO TO 820
 12 ipcmp1 = ipcmp + npcmp
 lcore  = lcore - npcmp
 IF (ipcmp1 >= lcore) GO TO 820
 
!     LOCATE PCOMP1 DATA AND READ INTO CORE
 
 14 CALL locate (*18,core(buf2),pcomp1,flag)
 
 ipcmp1 = ipcmp + npcmp
 CALL READ (*20,*16,pcomps,core(ipcmp1),lcore,0,npcmp1)
 GO TO 820
 16 ipcmp2 = ipcmp1 + npcmp1
 lcore  = lcore  - npcmp1
 IF (ipcmp2 >= lcore) GO TO 820
 
!     LOCATE PCOMP2 DATA AND READ INTO CORE
 
 18 CALL locate (*20,core(buf2),pcomp2,flag)
 
 ipcmp2 = ipcmp1 + npcmp1
 CALL READ (*20,*20,pcomps,core(ipcmp2),lcore,0,npcmp2)
 GO TO 820
 
 20 lpcomp = npcmp + npcmp1 + npcmp2
 
 lcore = lcore - npcmp2
 IF (lcore <= 0) GO TO 820
 
 CALL CLOSE (pcomps,1)
 
 
 25 CALL gopen (ecpt,core(buf2),0)
 IF (itemp == 0) THEN
   GO TO    40
 END IF
 30 ipm = gptt
 CALL OPEN (*750,gptt,core(buf3),0)
 
!     BRING IN MAT ETC
 
 CALL READ (*860,*810,gptt,tlist(1),  -2,0,ntlist)
 CALL READ (*860,*40 ,gptt,tlist(1),1080,1,ntlist)
 WRITE  (outpt,35) ufm
 35 FORMAT (a23,' 4013, PROBLEM LIMITATION OF 360 TEMPERATURE SETS ',  &
     ' HAS BEEN EXCEEDED.')
 n1 = -37
 GO TO 760
 40 IF (ideft /= 0) CALL gopen (edt,core(buf4),0)
 nloop = ideft + itemp
 IF (ideft /= 0) ldefm = 0
 
!     INITIALIZE MATERIAL ROUTINE
 
 imat  = nrowsp + lpcomp
 lcore = buf5   - imat
 CALL premat (core(imat+1),core(imat+1),core(buf5),lcore,nmat, mpt,dit)
 nstart = imat  + nmat
 lcore  = lcore - nstart
 IF (lcore <= 0) GO TO 820
 IF (ideft /= 0) ldefm = 0
 
 DO  illop = 1,nloop
   
   idefm = ilist(illop)
   IF (itemp > 0) THEN
     GO TO    70
   ELSE
     GO TO    75
   END IF
   70 CALL REWIND (gptt)
   
   75 IF (nnn == 1) GO TO 95
   
!     BRING SIL INTO CORE
   
   IF (lcore < 0) GO TO 820
   CALL gopen (sil,core(buf5),0)
   ipm = sil
   CALL READ (*860,*80,sil,core(nstart+1),lcore,1,nsil)
   GO TO 820
   80 CALL CLOSE (sil,1)
   lcore  = lcore  - nsil
   nstart = nstart + nsil
   
!     READ CSTM INTO OPEN CORE AND MAKE INITIAL CALLS TO PRETRD/PRETRS
   
   IF (lcore < 0) GO TO 820
   CALL OPEN (*90,cstm,core(buf5),0)
   icm = 0
   CALL skprec (cstm,1)
   ipm = cstm
   CALL READ (*860,*85,cstm,core(nstart+1),lcore,1,ncstm)
   GO TO 820
   85 CONTINUE
   
!     FOR THOSE SUBROUTINES WHICH USE BASGLB INSTEAD OF TRANSS/TRANSD,
!     WE NEED TO REPOSITION THE CSTM FILE AND LEAVE THE GINO BUFFER
!     AVAILABLE FOR LATER CALLS TO READ BY SUBROUTINE BASGLB.
   
   CALL REWIND (cstm)
   CALL skprec (cstm,1)
   
   CALL pretrd (core(nstart+1),ncstm)
   CALL pretrs (core(nstart+1),ncstm)
   
   lcore  = lcore  - ncstm
   nstart = nstart + ncstm
   IF (lcore <= 0) GO TO 820
   
   90 nnn = 1
   95 IF (itemp > 0) THEN
     GO TO    99
   ELSE
     GO TO   150
   END IF
   
   99 DO  i = 1,ntlist,3
     IF (idefm == tlist(i)) GO TO 110
   END DO
   
!     THERMAL LOAD NOT FOUND IN GPTT
   
   iparm(2) = iparm(1)
   iparm(1) = igptt
   CALL mesage (-32,idefm,iparm(1))
   110 TO = tlist(i+1)
   IF (tlist(i+2) == 0) GO TO 140
   i = tlist(i+2)
   DO  j = 1,i
     CALL fwdrec (*800,gptt)
   END DO
   
!     READ SETID AND VERIFY CORRECT RECORD.  FAILSAFE
   
   CALL READ (*121,*121,gptt,iddd,1,0,dummy)
   IF (iddd == idefm) GO TO 125
   121 WRITE  (outpt,122) idefm
   122 FORMAT (98H0*** system fatal error 4014, routine edtl detects bad &
       data on temperature data block for set id =,i9)
   n1 = -61
   GO TO 760
   125 record = .true.
   GO TO 150
   
!     THE GPTT (ELEMENT TEMPERATURE TABLE) IS NOW POSITIONED TO THE
!     TEMPERATURE DATA FOR THE SET REQUESTED.  SUBROUTINE SSGETD WILL
!     READ THE DATA.
   
   140 CONTINUE
   record = .false.
   
   150 CONTINUE
   CALL CLOSE (cstm,1)
   CALL OPEN (*151,cstm,core(buf5),0)
   CALL skprec (cstm,1)
   icm = 0
   151 DO  i = 1,nrowsp
     core(i) =  0.0
   END DO
   
!     INITIALIZE /SSGETT/ VARIABLES
   
   oldel  = 0
   eorflg = .false.
   endid  = .true.
   bufflg = .false.
   
!     ELEMENT CALL PROCESSING
   
   
!     READ THE ELEMENT TYPE
   
   170 CALL READ (*710,*830,ecpt,eltype,1,0,flag)
   IF (eltype >= 1 .AND. eltype <= nelems) GO TO 174
   CALL mesage (-7,0,NAME)
   172 WRITE  (outpt,173) swm,eltype
   173 FORMAT (a27,' 4015, ELEMENT THERMAL AND DEFORMATION LOADING NOT ',  &
       'COMPUTED FOR ILLEGAL ELEMENT TYPE',i9, /34X, 'IN MODULE SSG1.')
   GO TO 610
   174 idx    = (eltype-1)*incr
   jltype = 2*eltype - ipr
   nwords = NE(idx+12)
   
!     READ AN ENTRY FOR ONE ELEMENT FROM ECPT
   
   175 CALL READ (*840,*170,ecpt,mecpt(1),nwords,0,flag)
   IF (itemp /= 0) GO TO 176
   
!     ELEMENT DEFORMATION LOAD
   
   IF (idefm /= ldefm) CALL fedtst (idefm)
   ldefm = idefm
   IF (eltype == crod  ) GO TO 180
   IF (eltype == ctube ) GO TO 200
   IF (eltype == conrod) GO TO 190
   IF (eltype == cbar  ) GO TO 210
   GO TO 610
   
!     THERMAL LOAD
   
   176 CONTINUE
   
!     BRANCH TO THE DESIRED ELEMENT TYPE
   
   local = jltype - 100
   IF (local > 0) THEN
     GO TO   178
   END IF
   
   
!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION
   
!             1 CROD      2 CBEAM     3 CTUBE     4 CSHEAR    5 CTWIST
   177 GO TO( 180,  180,  172,  172,  200,  200,  360,  360,  610,  610 &   
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
   ,      270,  270,  240,  240,  250,  250,  220,  220,  190,  190 &   
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
   ,      610,  610,  610,  610,  610,  610,  610,  610,  260,  260 &   
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
   ,      230,  230,  280,  280,  300,  300,  290,  290,  610,  610 &   
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
   ,      610,  610,  610,  610,  610,  610,  610,  610,  610,  610 &   
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
   ,      610,  610,  610,  610,  610,  610,  610,  610,  610,  610 &   
!            31 PLOTEL   32 CREACT   33 CQUAD3   34 CBAR     35 CCONE  &
   ,      610,  610,  172,  172,  172,  172,  210,  210,  350,  350 &   
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
   ,      320,  320,  330,  330,  340,  340,  390,  390,  400,  400 &   
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
   ,      410,  410,  420,  420,  610,  610,  610,  610,  610,  610 &   
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
   ,      610,  610,  610,  610,  610,  610,  610,  610,  610,  610 &
    ), jltype
   
!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
   178 GO TO( 610,  610,  610,  610,  553,  553,  554,  554,  555,  555 &   
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
   ,      556,  556,  557,  557,  558,  558,  559,  559,  560,  560 &   
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQUAD4   65 CIHEX1  &
   ,      561,  561,  562,  562,  563,  563,  564,  564,  425,  425 &   
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
   ,      425,  425,  425,  425,  172,  172,  172,  172,  428,  428 &   
!             71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
   ,       429,  429,  172,  172,  430,  430,  431,  431,  432,  432 &
!             76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
   ,       172,  172,  172,  172,  172,  172,  172,  172,  433,  433 &   
!             81 CELBOW   82 FTUBE    83 TRIA3  &
   ,       172,  172,  610,  610,  566,  566 &
    ), local
   
!     ROD
   
   180 CALL rod
   GO TO 175
   
!     CONROD
   
   190 GO TO 180
   
!     TUBE
   
   200 GO TO 180
   
!     BAR
   
   210 CALL bar (core(1),idefm,itemp,ideft)
   GO TO 175
   
!     TRMEM
   
   220 CALL ssgetd (mecpt(1),ti(1),0)
   CALL trimem (0,ti,core(1))
   GO TO 175
   
!     QDMEM
   
   230 CALL ssgetd (mecpt(1),ti(1),0)
   CALL qdmem  (ti,core(1))
   GO TO 175
   
!     TRBSC
   
   240 CALL ssgetd (mecpt(1),ti(1),0)
   CALL trbsc  (0,ti)
   GO TO 175
   
!     TRPLT
   
   250 CALL ssgetd (mecpt(1),ti(1),0)
   CALL trplt  (ti)
   GO TO 175
   
!     QDPLT
   
   260 CALL ssgetd (mecpt(1),ti(1),0)
   CALL qdplt  (ti)
   GO TO 175
   
!     TRIA1
   
   270 kk = 1
   GO TO 301
   
!     TRIA2
   
   280 kk = 2
   GO TO 301
   
!     QUAD1
   
   290 kk = 3
   GO TO 301
   
!     QUAD2
   
   300 kk = 4
   301 CALL ssgetd (mecpt(1),ti(1),0)
   CALL triqd  (kk,ti(1))
   GO TO 175
   
!     TRIARG
   
   320 CALL ssgetd (mecpt(1),ti(1),3)
   CALL ttrirg (ti(2),core(1))
   GO TO 175
   
!     TRAPRG
   
   330 CALL ssgetd (mecpt(1),ti(1),4)
   CALL ttrapr (ti(2),core(1))
   GO TO 175
   
!     TORDRG
   
   340 CALL ssgetd (mecpt(1),ti(1),2)
   CALL ttordr (ti(2),core(1) )
   GO TO 175
   
!     CONE
   
   350 CALL ssgetd (mecpt(1),ti(1),2)
   CALL cone   (ti(2),core(1))
   GO TO 175
   
!     SHEAR PANEL
   
   360 CALL tshear
   GO TO 175
   
!     TETRA
   
   390 CALL ssgetd (mecpt(1),ti(1),4)
   CALL tetra  (ti(2),core(1),0)
   GO TO 175
   
!     WEDGE
   
   400 iijj = 1
   npts = 6
   GO TO 421
   
!     HEXA1
   
   410 iijj = 2
   npts = 8
   GO TO 421
   
!     HEXA2
   
   420 iijj = 3
   npts = 8
   421 CALL ssgetd (mecpt(1),ti(1),npts)
   CALL solid  (ti(2),core(1),iijj)
   GO TO 175
   
!     IHEX1, IHEX2, IHEX3
   
   425 npts=12*(eltype-64)-4
   CALL ssgetd (mecpt(1),ti(1),npts)
   CALL ihex   (ti(1),core(1),eltype-64)
   GO TO 175
   
!     TRIAAX
   
   428 CALL ssgetd (mecpt,ti,3)
   CALL trttem (ti(2),core)
   GO TO 175
   
!     TRAPAX
   
   429 CALL ssgetd (mecpt,ti,4)
   CALL tpztem (ti(2),core)
   GO TO 175
   
!     TRIM6
   
   430 CALL ssgetd (mecpt(1),ti,6)
   CALL tlodm6 (ti(1))
   GO TO 175
   
!     TRPLT1
   
   431 CALL ssgetd (mecpt(1),ti,0)
   CALL tlodt1 (ti(1),ti(1))
   GO TO 175
   
!     TRSHL
   
   432 CALL ssgetd (mecpt(1),ti(1),0)
   CALL tlodsl (ti(1),ti(1))
   GO TO 175
   
!     IS2D8
   
   433 CALL ssgetd (mecpt(1),ti(1),8)
   CALL tis2d8 (ti(2),core)
   GO TO 175
   
!     DUMMY ELEMENTS
   
   553 CALL dum1 (core(1))
   GO TO 175
   554 CALL dum2 (core(1))
   GO TO 175
   555 CALL dum3 (core(1))
   GO TO 175
   556 CALL dum4 (core(1))
   GO TO 175
   557 CALL dum5 (core(1))
   GO TO 175
   558 CALL dum6 (core(1))
   GO TO 175
   559 CALL dum7 (core(1))
   GO TO 175
   560 CALL dum8 (core(1))
   GO TO 175
   561 CALL dum9 (core(1))
   GO TO 175
   
!     QDMEM1
   
   562 CALL ssgetd (mecpt(1),ti(1),0)
   CALL qdmm1  (ti,core(1))
   GO TO 175
   
!     QDMEM2
   
   563 CALL ssgetd (mecpt(1),ti(1),0)
   CALL qdmm2  (ti,core(1))
   GO TO 175
   
!     QUAD4
   
   564 DO  iti = 1,7
     ti(iti) = 0.0
   END DO
   CALL ssgetd (mecpt(1),ti,4)
   IF (ipr /= 0) CALL tlqd4s
   IF (ipr == 0) CALL tlqd4d
   GO TO 175
   
!     TRIA3
   
   566 DO  iti = 1,7
     ti(iti) = 0.0
   END DO
   CALL ssgetd (mecpt(1),ti,3)
   IF (ipr /= 0) CALL tltr3s
   IF (ipr == 0) CALL tltr3d
   GO TO 175
   
!     NO LOAD, SKIP THE ECPT ENTRY ONLY
   
   610 CALL fwdrec (*840,ecpt)
   GO TO 170
   
!     PACK THE LOAD VECTOR FROM CORE TO OUTPUT DATA BLOCK -PG-
   
   710 CALL pack (core,pg(1),pg)
   CALL REWIND (ecpt)
   CALL fwdrec (*840,ecpt)
   IF (ideft /= 0 .AND. idefm /= 0) CALL fedted (idefm)
   
 END DO
 
 IF (noedt  == 0) CALL CLOSE (edt ,1)
 IF (nogptt == 0) CALL CLOSE (gptt,1)
 IF (icm    == 0) CALL CLOSE (cstm,1)
 CALL CLOSE (ecpt,1)
 RETURN
 
 750 n1 = -1
 760 CALL mesage (n1,ipm,iparm)
 800 ipm = gptt
 GO TO 860
 810 n1 = -3
 GO TO 760
 820 n1 = -8
 GO TO 760
 830 ipm = ecpt
 GO TO 810
 840 ipm = ecpt
 860 n1 = -2
 GO TO 760
END SUBROUTINE edtl
