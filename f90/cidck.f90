SUBROUTINE cidck (z,buf,nopen)
     
!     BULK DATA CARD COORDINATE CHECK
!     THIS ROUTINE IS CALLED ONLY BY IFP, IN LINK1
 
!     WRITTEN BY G.CHAN/UNISYS   9/1989
 
!     LIST OF NASTRAN BULK DATA CARDS THAT REFERENCE COORDINATE CID -
 
!     BULK DATA      CID      NO. OF     GINO         LOCATE
!     CARD           FIELD    WORDS      FILE         INDEX
!     ----------   -------   -------   ---------   ------------
!     AXIF            1          1       AXIC         8815,88
!     BFIELD          1          2       GEOM1        3101,31
!     CEMLOOP*        13        13       GEOM3        3109,31
!     CONM2           3         13       GEOM2        1501,15
!     CORD1C          1          6       GEOM1        1701,17
!     CORD1R          1          6       GEOM1        1801,18
!     CORD1S          1          6       GEOM1        1901,19
!     CORD2C          1,4       13       GEOM1        2001,20
!     CORD2R          1,4       13       GEOM1        2101,21
!     CORD2S          1,4       13       GEOM1        2201,22
!     FORCE           3          7       GEOM3        4201,42
!     GEMLOOP*        3          -       GEOM3        3309,33
!     GRAV            2          6       GEOM3        4401,44
!     GRID/GRDSET     2,6        8       GEOM1        4501,45
!     GRIDB           3          5       GEOM1        8115,81
!     MDIPOLE*        2         10       GEOM3        3509,35
!     MOMENT          3          7       GEOM3        4801,48
!     PIHEX           3          7       EPT          7002,70
!     PLOAD4          9         12       GEOM3        6709,67
!     REMFLUX*        2          -       GEOM3        3409,34
!WKBR 2/95 SPR94015 RFORCE    3     8       GEOM3        5509,55
!     RFORCE          3          7       GEOM3        5509,55
!     SPCFLD*         2          -       GEOM3        3209,32
 
!     * THE CID'S ON THESE CARDS CURRENTLY MUST BE ZERO OR BLANK, AND
!       WERE CHECKED ALREADY IN IFS4P. THEREFORE THEY ARE NOT CHECKED
!       HERE.
 
 
 INTEGER, INTENT(IN OUT)                  :: z(1)
 INTEGER, INTENT(IN OUT)                  :: buf(1)
 INTEGER, INTENT(IN)                      :: nopen
 EXTERNAL andf
 LOGICAL :: abort
 INTEGER :: axif(2),   bfield(2), conm2(2),  cord(2),   force(2),  &
     grav(2),   grid(2),   gridb(2),  moment(2), pihex(2),  &
     pload4(2), rforce(2),  trl(7),    NAME(2)
 CHARACTER (LEN=7) :: cc,    pcc,       caxif,     cbfiel,    cconm2,  &
     ccord2,    cforce,    cgrav,     cgrid,     cgridb,  &
     cmment,    cpihex,    cplod4,    crforc
 CHARACTER (LEN=23) :: ufm
 COMMON   /xmssg /   ufm
 COMMON   /system/   ibuf,      nout,      abort
 COMMON   /two   /   two(1)
 DATA     geom1,     geom2,     geom3,     ept        /  &
     201,       208,       209,       202        /,  &
     axic,      NAME,      pcd,       pcc        /  &
     215,       4HCIDC,    2HK , 0,   'XXXX   '  /
 DATA     axif,      caxif  /   8815,88,   'AXIF   '  /  &
     bfield,    cbfiel /   3101,31,   'BFIELD '  /,  &
     conm2,     cconm2 /   1501,15,   'CONM2  '  /  &
     cord,      ccord2 /   1601,16,   'CORD2  '  /,  &
     force,     cforce /   4201,42,   'FORCE  '  /  &
     grav,      cgrav  /   4401,44,   'GRAV   '  /,  &
     grid,      cgrid  /   4501,45,   'GRID   '  /  &
     gridb,     cgridb /   8115,81,   'GRIDB  '  /,  &
     moment,    cmment /   4801,48,   'MOMENT '  /  &
     pihex,     cpihex /   7002,70,   'PIHEX  '  /,  &
     pload4,    cplod4 /   6709,67,   'PLOAD4 '  /  &
     rforce,    crforc /   5509,55,   'RFORCE '  /
 
 
!     OPEN GEOM1 AND SAVE ALL COORDINATE IDS IN Z(1) THRU Z(NCORD)
!     AND REFERENCED COORD ID IN Z(NRID) THRU Z(NOPEN). NOPEN IS
!     LENGTH OF THE AVAILABLE OPEN CORE.
!     SORT AND CHECK ID UNIQUENESS
 
 ncord= 1
 nrid = nopen
 FILE = geom1
 CALL preloc (*960,buf,geom1)
 k = 6
 DO  i = 1,6
   cord(1) = cord(1)+100
   cord(2) = cord(2)+1
   IF (i == 4) k = 13
   CALL locate (*20,buf,cord(1),m)
   10   CALL READ (*910,*20,geom1,z(ncord),k,0,m)
   ncord = ncord+1
   IF (i < 4 .OR. z(ncord+2) == 0) GO TO 10
   z(nrid) = z(ncord+2)
   nrid  = nrid-1
   GO TO 10
   20   CONTINUE
 END DO
 ncord = ncord-1
 nrid  = nrid +1
 IF (ncord <= 1) GO TO 60
 CALL sort (0,0,1,1,z(1),ncord)
 j = 1
 DO  i = 2,ncord
   IF (z(i) /= z(i-1)) GO TO 40
   CALL page2 (-2)
   WRITE  (nout,30) ufm,z(i)
   30   FORMAT (a23,' 328, DUPLICATE COORDINATE ID',i9)
   CYCLE
   40   j = j+1
   z(j) = z(i)
 END DO
 ncord = j
 
!     IF CORD2C/R/S CARDS ARE PRESENT, CHECK REFERENCE COORDINATE ID
 
 60   IF (nrid > nopen) GO TO 100
 cc = ccord2
 loop90:  DO  i = nrid,nopen
   cid = z(i)
   IF (ncord <= 0) GO TO 80
   DO  j = 1,ncord
     IF (cid == z(j)) CYCLE loop90
   END DO
   80   CALL page2 (-2)
   WRITE (nout,830) ufm,cid,cc
   abort = .true.
 END DO loop90
 
!     DOUBLE THE COORDINATE ID ARRAY FOR 'CIRCULAR' SEARCH, AND MOVE
!     THE ARRAY TO HIGH END OF OPEN CORE, Z(II) THRU Z(NOPEN-1)
 
 100  ii = nopen-2*ncord-1
 IF (ncord == 0) GO TO 120
 DO  i = 1,ncord
   z(ii+i      ) = z(i)
   z(ii+i+ncord) = z(i)
 END DO
 120  nz = ii
 im = ii+ncord
 ii = ii+1
 z(nopen) = -999
 
!     CHECK CID ON GRID CARDS
 
 cc = cgrid
 CALL locate (*190,buf,grid(1),m)
 nzx = (nz/8)*8
 130  CALL READ (*910,*140,geom1,z(1),nzx,0,m)
 m = nzx
 IF (m <= 0) GO TO 190
 140  pvcid = 0
 ASSIGN 150 TO irtn
 i = -6
 150  i = i+8
 IF (i > m) GO TO 160
 cid = z(i)
 IF (cid /= 0 .AND. cid /= pvcid) GO TO 790
 GO TO 150
 160  pvcid = 0
 ASSIGN 170 TO irtn
 i = -2
 170  i = i+8
 IF (i > m) GO TO 180
 cid = z(i)
 IF (cid /= 0 .AND. cid /= pvcid) GO TO 790
 GO TO 170
 180  IF (m == nzx) GO TO 130
 
!     CHECK GRIDB CARDS
 
 190  cc = cgridb
 CALL locate (*240,buf,gridb(1),m)
 nzx = (nz/5)*5
 200  CALL READ (*910,*210,geom1,z(1),nzx,0,m)
 m = nzx
 IF (m <= 0) GO TO 240
 210  pvcid = 0
 ASSIGN 220 TO irtn
 i = -2
 220  i = i+5
 IF (i > m) GO TO 230
 cid = z(i)
 IF (cid /= 0 .AND. cid /= pvcid) GO TO 790
 GO TO 220
 230  IF (m == nzx) GO TO 200
 
!     CHECK BFIELD CARDS
 
 240  cc = cbfiel
 CALL locate (*270,buf,bfield(1),m)
 CALL READ (*910,*250,geom1,z(1),nz,1,m)
 GO TO 930
 250  pvcid = 0
 ASSIGN 260 TO irtn
 i = -1
 260  i = i+2
 IF (i > m) GO TO 270
 cid = z(i)
 IF (cid /= 0 .AND. cid /= pvcid) GO TO 790
 GO TO 260
 
!     END OF GEOM1 PROCESSING
 
 270  CALL CLOSE (geom1,1)
 
 
!     CHECK THE PRESENCE OF CONM2 CARDS IN GEOM2
 
 FILE = geom2
 k = conm2(2)
 ASSIGN 300 TO jrtn
 GO TO 860
 300  IF (k == 0) GO TO 400
 
!     OPEN GEOM2, AND CHECK CONM2 CARDS
 
 cc = cconm2
 CALL preloc (*400,buf,geom2)
 CALL locate (*350,buf,conm2(1),m)
 nzx = (nz/13)*13
 310  CALL READ (*910,*320,geom2,z(1),nzx,0,m)
 m = nzx
 IF (m <= 0) GO TO 350
 320  pvcid = 0
 ASSIGN 330 TO irtn
 i = -10
 330  i = i+13
 IF (i > m) GO TO 340
 cid = z(i)
 IF (cid /= 0 .AND. cid /= pvcid) GO TO 790
 GO TO 330
 340  IF (m == nzx) GO TO 310
 
 350  CONTINUE
 CALL CLOSE (geom2,1)
 
 
!     CHECK THE PRESENCE OF BULK DATA CARDS IN GEOM3
!     (FORCE, MOMENT, RFORCE, GRAV AND PLOAD4)
 
 400  FILE = geom3
 k = force(2)
 ASSIGN 410 TO jrtn
 GO TO 860
 410  IF (k /= 0) GO TO 500
 k = moment(2)
 ASSIGN 420 TO jrtn
 GO TO 870
 420  IF (k /= 0) GO TO 500
 k = rforce(2)
 ASSIGN 430 TO jrtn
 GO TO 870
 430  IF (k /= 0) GO TO 500
 k = grav(2)
 ASSIGN 440 TO jrtn
 GO TO 870
 440  IF (k /= 0) GO TO 500
 k = pload4(2)
 ASSIGN 450 TO jrtn
 GO TO 870
 450  IF (k /= 0) GO TO 500
 GO TO 650
 
!     OPEN GEOM3, AND CHECK CID ON BULK DATA CARDS
 
 500  CALL preloc (*650,buf,geom3)
 CALL locate (*510,buf,force(1),m)
 cc = cforce
 ib = 3
 ic = 7
 ASSIGN 510 TO krtn
 GO TO 600
 510  CALL locate (*520,buf,moment(1),m)
 cc = cmment
 ib = 3
 ic = 7
 ASSIGN 520 TO krtn
 GO TO 600
 520  CALL locate (*530,buf,rforce(1),m)
 cc = crforc
 ib = 3
!WKBR 2/95 SPR94015      IC = 8
 ic = 7
 ASSIGN 530 TO krtn
 GO TO 600
 530  CALL locate (*540,buf,grav(1),m)
 cc = cgrav
 ib = 2
 ic = 6
 ASSIGN 540 TO krtn
 GO TO 600
 540  CALL locate (*550,buf,pload4(1),m)
 cc = cplod4
 ib = 9
 ic = 12
 ASSIGN 550 TO krtn
 GO TO 600
 550  CONTINUE
 GO TO 630
 
 600  CALL READ (*910,*610,geom3,z(1),nz,1,m)
 GO TO 930
 610  ASSIGN 620 TO irtn
 i = ib-ic
 620  i = i +ic
 IF (i > m) GO TO krtn, (510,520,530,540,550)
 cid = z(i)
 IF (cid /= 0) GO TO 800
 GO TO 620
 
 630  CALL CLOSE (geom3,1)
 
 
!     CHECK THE PRESENCE OF PIHEX CARD IN EPT. IF PRESENT, OPEN EPT,
!     AND CHECK CID ON PIHEX CARDS
 
 650  FILE = ept
 k = pihex(2)
 ASSIGN 660 TO jrtn
 GO TO 860
 660  IF (k == 0) GO TO 700
 CALL preloc (*700,buf,ept)
 CALL locate (*690,buf,pihex(1),m)
 CALL READ (*910,*670,ept,z(1),nz,1,m)
 GO TO 930
 670  cc = cpihex
 ASSIGN 680 TO irtn
 i = -4
 680  i = i+7
 IF (i > m) GO TO 690
 cid = z(i)
 IF (cid /= 0) GO TO 800
 GO TO 680
 690  CALL CLOSE (ept,1)
 
 
!     CHECK THE PRESENCE OF AXIF CARD IN AXIC. IF PRESENT, OPEN AXIC,
!     AND CHECK CID ON AXIF CARD. ONLY ONE AXIF CARD EXISTS
 
 700  FILE = axic
 k = axif(2)
 ASSIGN 710 TO jrtn
 GO TO 860
 710  IF (k == 0) GO TO 750
 CALL preloc (*750,buf,axic)
 CALL locate (*730,buf,axif(1),m)
 CALL READ (*910,*720,axic,cid,1,1,m)
 720  cc = caxif
 ASSIGN 730 TO irtn
 IF (cid /= 0) GO TO 800
 730  CALL CLOSE (axic,1)
 
 750  RETURN
 
 
!     INTERNAL ROUTINE TO LOOK FOR CID MATCH
!     CID ARRAY (DOUBLE) IS AT HIGH END OF CORE, Z(II) THRU Z(NOPEN)
 
 790  pvcid = cid
 800  IF (cid == z(ii)) GO TO 850
 IF (ncord  <=  1) GO TO 820
 ie = ii+ncord-1
 DO  j = ii,ie
   IF (cid == z(j)) GO TO 840
 END DO
 820  IF (cc == pcc .AND. cid == pcd) GO TO 850
 CALL page2 (-2)
 WRITE  (nout,830) ufm,cid,cc
 830  FORMAT (a23,' 328, UNDEFINED COORDINATE',i9,' IS REFERENCED BY A '  &
     ,       a7,' CARD')
 pcc = cc
 pcd = cid
 abort = .true.
 GO TO 850
 840  ii = j
 IF (ii > im) ii = ii-ncord
 850  GO TO irtn, (150,170,220,260,330,620,680,730)
 
 
!     INTERNAL ROUTINE TO CHECK THE PRESENCE OF A PARTICULAR BULK DATA
!     CARD
 
 860  trl(1) = FILE
 CALL rdtrl (trl(1))
 870  IF (trl(1) < 0) GO TO 880
 j = (k-1)/16
 l = k-16*j
 IF (andf(trl(j+2),two(l+16)) /= 0) GO TO 890
 880  k = 0
 890  GO TO jrtn, (300,410,420,430,440,450,660,710)
 
!     ERRORS
 
 910  j = -2
 GO TO 950
 930  j = -8
 950  CALL mesage (j,FILE,NAME)
 960  RETURN
END SUBROUTINE cidck
