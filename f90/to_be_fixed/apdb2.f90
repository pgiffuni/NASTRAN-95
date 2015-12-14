SUBROUTINE apdb2 (ibuf1,ibuf2,next,left,nstns,nlines,xsign,  &
        lcstm,acstm,nodex,nodei,isilc,xyzb)
     
!     GENERATE GTKA TRANSFORMATION MATRIX FOR SWEPT TURBOPROP
!     BLADES (AERODYNAMIC THEORY NUMBER 7).
 
 
 INTEGER, INTENT(IN OUT)                  :: ibuf1
 INTEGER, INTENT(IN OUT)                  :: ibuf2
 INTEGER, INTENT(IN OUT)                  :: next
 INTEGER, INTENT(IN OUT)                  :: left
 INTEGER, INTENT(IN)                      :: nstns
 INTEGER, INTENT(IN)                      :: nlines
 REAL, INTENT(IN OUT)                     :: xsign
 INTEGER, INTENT(IN OUT)                  :: lcstm
 REAL, INTENT(IN OUT)                     :: acstm(1)
 INTEGER, INTENT(OUT)                     :: nodex(1)
 INTEGER, INTENT(OUT)                     :: nodei(1)
 INTEGER, INTENT(OUT)                     :: isilc(1)
 REAL, INTENT(OUT)                        :: xyzb(4,nstns)
 EXTERNAL        andf
 LOGICAL :: multi,omit,single,debug
 INTEGER :: gm,GO,gtka,scr1,scr2,core,idata(7),  &
     um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,uset1,  &
     gtkg,gknb,gkm,gkab,gkf,gks,gko,gkn,gsize,  &
     andf,rd,rdrew,wrt,wrtrew,clsrew,tgkg(7)
 DIMENSION       itrl(7), iz(1),z(1),rdata(7),  &
     ta(3,3),tbl(3),tbla(3),tblt(3),tblr(3),
 COMMON /system/ ksystm(54),iprec
 COMMON /two   / itwo(32)
 COMMON /zzzzzz/ core(1)
 COMMON /zblpkx/ ap(4),ii
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /patx  / lc,n,no,ny,uset1,ibc(7)
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew
 COMMON /apdbug/ debug
 EQUIVALENCE     (z(1),core(1))
 EQUIVALENCE     (z(1),iz(1)), (idata(1),rdata(1))
 DATA    single, multi,omit /.true.,.true.,.true./
 
 uset = 102
 gm   = 106
 GO   = 107
 gtka = 204
 scr1 = 301
 scr2 = 302
 gknb = 303
 gkm  = 304
 gkab = 305
 itrl(1) = uset
 CALL rdtrl (itrl)
 gsize = itrl(3)
 IF (andf(itrl(5),itwo(um)) == 0) multi = .false.
 IF (andf(itrl(5),itwo(us)) == 0) single= .false.
 IF (andf(itrl(5),itwo(uo)) == 0) omit  = .false.
 IF (.NOT.(multi .OR. single .OR. omit)) scr2 = gtka
 gtkg = scr2
 
!     OPEN SCR1 TO READ BLADE NODE DATA
 
!                         T
!     OPEN SCR2 TO WRITE G   MATRIX OF ORDER (GSIZE X KSIZE)
!                         KG
 
 CALL gopen (scr1,z(ibuf1),rdrew)
 CALL gopen (gtkg,z(ibuf2),wrtrew)
 tgkg(1) = gtkg
 tgkg(2) = 0
 tgkg(3) = gsize
 tgkg(4) = 2
 tgkg(5) = 1
 tgkg(6) = 0
 tgkg(7) = 0
 
!     SET-UP CALL TO TRANSS VIA PRETRS
 
 IF (lcstm > 0) CALL pretrs (acstm,lcstm)
 
!     LOOP ON STREAMLINES
 
 DO  nline = 1,nlines
   
!     READ STREAMLINE NODE DATA FROM SCR1
   
   DO  nst = 1,nstns
     CALL fread (scr1,idata,7,0)
     IF (debug) CALL bug1 ('SCR1 IDATA',10,idata,7)
     nodex(nst) = idata(1)
     nodei(nst) = idata(2)
     isilc(nst) = idata(3)
     xyzb(1,nst)= rdata(4)
     xyzb(2,nst)= rdata(5)
     xyzb(3,nst)= rdata(6)
     xyzb(4,nst)= rdata(7)
   END DO
   
!     GENERATE BASIC TO LOCAL TRANSFORMATION MATRIX FOR THIS STREAMLINE
   
   CALL apdb2a (nlines,nline,scr1,nstns,xsign,xyzb(2,1),  &
       xyzb(2,nstns),tblt,tblr)
   
!     SET TRANSFORMATION TO TRANSLATION FIRST
   
   DO  nn = 1,3
     tbl(nn) = tblt(nn)
   END DO
   
!     LOOP FOR TRANSLATION THEN ROTATION
   
   ndeg = 0
   DO  nloop = 1,2
     IF (debug) CALL bug1 ('MAT-TBL   ' ,18,tbl,3)
     
!     LOOP ON COMPUTING STATIONS
     
     DO  ncs = 1,nstns
       
!     LOCATE GLOBAL TO BASIC TRANSFORMATION MATRIX
       
       rdata(1) = xyzb(1,ncs)
       IF (lcstm == 0 .OR. idata(1) == 0) GO TO 20
       CALL transs (xyzb(1,ncs),ta)
       CALL gmmats (tbl,1,3,0,ta,3,3,0,tbla)
       GO TO 25
       20 tbla(1) = tbl(1)
       tbla(2) = tbl(2)
       tbla(3) = tbl(3)
       25 CONTINUE
       IF (debug) CALL bug1 ('MAT-TBLA  ',25,tbla,3)
       
!     COMPUTE LOCATION IN G-SET USING SIL
!     KODE = 1 FOR GRID POINT
!     KODE = 2 FOR SCALAR POINT (NOT ALLOWED, CHECK WAS MADE BY APDB)
       
       isil = isilc(ncs)/10
       CALL bldpk (1,1,gtkg,0,0)
       
!     OUTPUT GKG(TRANSPOSE) = GTKG
!     II IS ROW POSITION
       
       DO  icol = 1,3
         ii = isil + ndeg
         ap(1) = tbla(icol)
         IF (debug) CALL bug1 ('ISIL      ',28,isil,1)
         IF (debug) CALL bug1 ('MAT-AP    ',29,ap,1)
         CALL zblpki
         isil = isil + 1
       END DO
       CALL bldpkn (gtkg,0,tgkg)
     END DO
     
!     CHANGE BASIC TO LOCAL TRANSFORMATION TO ROTATION
     
     DO  nn = 1,3
       tbl(nn) = tblr(nn)
     END DO
     ndeg = 3
   END DO
 END DO
 CALL CLOSE (scr1,clsrew)
 CALL CLOSE (gtkg,clsrew)
 CALL wrttrl (tgkg)
 
!     CREATE GTKA MATRIX
 
 IF (multi .OR. single .OR. omit) GO TO 60
 GO TO 100
 60 CONTINUE
 lc  = korsz(core)
 gkf = gknb
 gks = gkm
 gko = gks
 uset1 = uset
 
!     REDUCE TO N-SET IF MULTI POINT CONSTRAINTS
 
 gkn = gtkg
 IF (.NOT.multi) GO TO 70
 IF (.NOT.single .AND. .NOT.omit) gkn = gtka
 CALL calcv (scr1,ug,un,um,core)
 CALL ssg2a (gtkg,gknb,gkm,scr1)
 CALL ssg2b (gm,gkm,gknb,gkn,1,iprec,1,scr1)
 
!     PARTITION INTO F-SET IF SINGLE POINT CONSTRAINTS
 
 70 IF (.NOT.single) GO TO 80
 IF (.NOT.omit  ) gkf = gtka
 CALL calcv (scr1,un,uf,us,core)
 CALL ssg2a (gkn,gkf,0,scr1)
 GO TO 90
 
!     REDUCE TO A-SET IF OMITS
 
 80 gkf = gkn
 90 IF (.NOT.omit) GO TO 100
 CALL calcv (scr1,uf,ua,uo,core)
 CALL ssg2a (gkf,gkab,gko,scr1)
 CALL ssg2b (GO,gko,gkab,gtka,1,iprec,1,scr1)
 100 RETURN
END SUBROUTINE apdb2
