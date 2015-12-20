SUBROUTINE flbelm
     
!     READS CFLSTR AND CFREE BULK DATA AND BUILDS INCORE TABLES TO
!     DESCRIBE THE CONNECTIVITY BETWEEN THE STRUCTURE AND FLUID
 
 LOGICAL :: error
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict   ,  &
     z        ,FILE     ,NAME(2)  ,mcb(7)   ,cflstr(2),  &
     card(10) ,id(3)    ,grid(4)  ,cfree(2) ,elm2d(7,3) ,               elmfl(4,3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /flbfil/ geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,nbgpdt   ,  &
     isil     ,nsil     ,igrav    ,ngrav    ,igrid    ,  &
     ngrid    ,ibuf1    ,ibuf2    ,ibuf3    ,ibuf4    , ibuf5
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf   ,nout
 COMMON /BLANK / nograv   ,nofree
 DATA    cflstr/ 7610,76/ ,cfree / 4810,48 /  ,mcb / 7*0 /
 DATA    NAME  / 4HFLBE , 4HLM   /
 
!     TWO DIMENSIONAL STRUCTURAL ELEMENTS DESCRIPTIONS
 
 DATA    n2d   / 7 /
 DATA    elm2d / & 
!                     TRIA1  TRIA2   TRMEM  QUAD1  QUAD2  QDMEM  SHEAR
!     1  IFP CARD NUMBERS
!     2  NUMBER OF GRIDS
!     3  NUMBER OF WORDS IN ECT RECORD
     52    ,53     ,56    ,57    ,58    ,60    ,61   ,  &
     3     ,3      ,3      ,4     ,4     ,4     ,4   ,  &
     6     ,6      ,6      ,7      ,7    ,7      ,6  /
 
!     FLUID ELEMENT DESCRIPTIONS
 
 DATA    nfl   / 4 /
 DATA    elmfl / & 
!                     FHEX1     FHEX2     FTETRA    FWEDGE
!    1  IFP CARD NUMBERS
!    2  NUMBER OF GRIDS
!    3  NUMBER OF WORDS IN ECT RECORD
     333      ,334      ,335      ,336     ,  &
     8        ,8        ,4        ,6       ,  &
     10       ,10       ,6        ,8       /
 
 
!     READ BGPDT INTO OPEN CORE
 
 ibgpdt = 1
 FILE   = bgpdt
 CALL gopen (bgpdt,z(ibuf1),0)
 nz = ibuf3 - 1
 CALL READ (*1002,*10,bgpdt,z(ibgpdt),nz,1,nbgpdt)
 GO TO 1008
 10 icore = ibgpdt + nbgpdt
 ngrdt = nbgpdt/4
 CALL CLOSE (bgpdt,1)
 
!     LOCATE CFLSTR CARDS ON GEOM2 AND READ THEM INTO ELEMENT TABLE
!     IN CORE.   ONE ELEMENT TABLE RECORD WILL LOOK AS FOLLOWS -
 
!                  WORD      DESCRIPTION
 
!                  1         STRUCTURE ELEMENT ID
!                  2         FLUID ELEMENT ID
!                  3-6       ZERO
!                  7         GRAV LOAD ID
 
 FILE = geom2
 CALL preloc (*1001,z(ibuf1),geom2)
 CALL locate (*1200,z(ibuf1),cflstr,id)
 ielmt = icore
 20 CALL READ (*1002,*40,geom2,id,2,0,n)
 30 CALL READ (*1002,*1003,geom2,ids,1,0,n)
 IF (ids < 0) GO TO 20
 IF (icore+7 >= ibuf3) GO TO 1008
 z(icore  ) = ids
 z(icore+1) = id(1)
 z(icore+2) = 0
 z(icore+3) = 0
 z(icore+4) = 0
 z(icore+5) = 0
 z(icore+6) = id(2)
 icore = icore + 7
 GO TO 30
 
 40 nelmt = icore - ielmt
 nelm  = nelmt/7
 
!     SORT ELEMENT TABLE BY STRUCTUREAL ELEMENT ID
 
 CALL sort (0,0,7,1,z(ielmt),nelmt)
 
!     READ ECT AND PROCESS 2D STRUCTURAL ELEMENTS
 
 FILE = ect
 CALL gopen (ect,z(ibuf2),0)
 50 CALL READ (*100,*1002,ect,card,3,0,n)
 DO  i = 1,n2d
   IF (card(3) == elm2d(i,1)) GO TO 70
 END DO
 
!     SKIP RECORD BECAUSE NOT ACCEPTABLE 2D ELEMENT TYPE
 
 CALL fwdrec (*1001,ect)
 GO TO 50
 
!     PROCESS THE 2D ELEMENT
 
 70 ngrds = elm2d(i,2)
 nwds  = elm2d(i,3)
 
!     READ DATA FOR ONE 2D ELEMENT
 
 80 CALL READ (*1001,*50,ect,card,nwds,0,n)
 
!     CHECK IF STRUCTURAL ELEMENT IS CONNECTED TO ANY FLUID ELEMENT
!     MAKE SURE BISLOC FINDS FIRST OF SEVERAL POSSIBLE ENTRIES
 
 CALL bisloc (*80,card(1),z(ielmt),7,nelm,jloc)
 82 IF (jloc == 1 .OR. z(ielmt+jloc-8) /= card(1)) GO TO 84
 jloc = jloc - 7
 GO TO 82
 
!     INSERT ELEMENT GRID POINTS INTO ELEMENT TABLE WORDS 3-6
 
 84 DO  i = 1,ngrds
   z(ielmt+jloc+i) = card(i+2)
 END DO
 IF (ngrds == 3) z(ielmt+jloc+4) = -1
 
!     CHECK IF NEXT ENTRY IS FOR THE SAME STRUCTURAL ELEMENT
 
 IF (jloc+7 >= nelmt .OR. z(ielmt+jloc+6) /= card(1)) GO TO 80
 jloc = jloc + 7
 GO TO 84
 
 100 CONTINUE
 
!     PASS THROUGH ELEMENT TABLE AND CHECK THAT EACH ENTRY HAS GRIDS.
!     ALSO SWITCH THE STRUCTURE AND FLUID ELEMENTS IN THE TABLE FOR
!     FUTURE WORD WITH FLUID ID.
 
 lelmt = ielmt + nelmt - 1
 DO  i = ielmt,lelmt,7
   ids  = z(i  )
   z(i) = z(i+1)
   IF (z(i+2) /= 0) GO TO 110
   error = .true.
   WRITE (nout,8002) ufm,ids
   ids  = 0
   110 z(i+1) = ids
 END DO
 
!     ALLOCATE AND ZERO THE GRID POINT CONNECTIVE TABLE AT THE BOTTOM
!     OF CORE
 
!     TABLE ENTRIES WILL BE AS FOLLOWS
 
!     POSITIVE LESS THEN 1,000,000  - NUMBER OF STRUCTURAL POINTS
!                                     CONNECTED TO THIS FLUID POINT
!     MULTIPLES OF 1,000,000        - NUMBER OF FREE SURFACE POINTS
!                                     CONNECTED TO THIS FLUID POINT
!     NEGATIVE                      - NUMBER OF STRUCTURAL POINTS
!                                     CONNECTED TO THIS STRUCTURAL
!                                     POINT
 
 igrid = ibuf3 - ngrdt - 1
 IF (igrid < icore) GO TO 1008
 ngrid = ngrdt
 lgrid = ibuf3 - 1
 DO  i = igrid,lgrid
   z(i) = 0
 END DO
 
!     LOCATE CFREE CARDS ON GEOM2 AND ADD THEM TO THE ELEMENT TABLE.
!     THESE ELEMENT RECORDS WILL APPEAR AS FOLLOWS
 
!                  WORD      DESCRIPTION
 
!                  1         FLUID ELEMENT ID
!                  2         -1
!                  3         FACE ID
!                  4-6       ZERO
!                  7         GRAV ID
 
 FILE = geom2
 CALL locate (*124,z(ibuf1),cfree,id)
 nofree = 1
 120 CALL READ (*1002,*130,geom2,id,3,0,n)
 IF (icore+7 >= igrid) GO TO 1008
 z(icore  ) = id(1)
 z(icore+1) = -1
 z(icore+2) = id(3)
 z(icore+3) = 0
 z(icore+4) = 0
 z(icore+5) = 0
 z(icore+6) = id(2)
 icore = icore + 7
 GO TO 120
 
!     NO CFREE CARDS - THIS IMPLIES THAT THERE WILL BE NO FREE SURFACE
 
 124 nofree = -1
 
!     COMPLETE CORE ALLOCATION FOR THIS PHASE
 
 130 nelmt = icore - ielmt
 nelm  = nelmt/7
 CALL CLOSE (geom2,1)
 
!     SORT ELEMENT TABLE BY FLUID ID
 
 CALL sort (0,0,7,1,z(ielmt),nelmt)
 
!     OPEN FBELM AND FRELM SCRATCH FILES
 
 CALL gopen (fbelm,z(ibuf1),1)
 CALL gopen (frelm,z(ibuf3),1)
 
!     READ ECT AND PROCESS FLUID ELEMENTS
 
 FILE = ect
 CALL REWIND (ect)
 CALL fwdrec (*1002,ect)
 140 CALL READ (*220,*1003,ect,card,3,0,n)
 DO  i = 1,nfl
   IF (card(3) == elmfl(i,1)) GO TO 160
 END DO
 
!     SKIP RECORD BECAUSE NOT FLUID ELEMENT TYPE
 
 CALL fwdrec (*1001,ect)
 GO TO 140
 
!     PRECESS FLUID ELEMENT
 
 160 ntype = elmfl(i,1)
 nwds  = elmfl(i,3)
 
!     READ DATA FOR ONE FLUID ELEMENT
 
 170 CALL READ (*1001,*140,ect,card,nwds,0,n)
 
!     FIND IF FLUID ELEMENT IS ON FREE SURFACE OR STRUCTURAL BOUNDARY.
!     MAKE SURE BISLOC FINDS THE FIRST OF SEVERAL POSSIBLE ENTRIES.
 
 CALL bisloc (*170,card(1),z(ielmt),7,nelm,jloc)
 175 IF (jloc == 1 .OR. z(ielmt+jloc-8) /= card(1)) GO TO 180
 jloc = jloc - 7
 GO TO 175
 
!     DETERMINE IF ENTRY IS EITHER A BOUNDARY OR FREE SURFACE
!     DESCRIPTION - IGNORE ENTRY IF IT WAS IN ERROR DURING STRUCTURAL
!     ELEMENT PROCESSING
 
 180 IF (z(ielmt+jloc) >  0) GO TO 190
 IF (z(ielmt+jloc) == -1) GO TO 200
 GO TO 210
 
!     THIS ENTRY DESCRIBES THE FLUID / STRUCTURE BOUNDARY - FIND THE
!     FLUID GRID POINTS WHICH COINCIDE WITH THE STRUCTURAL POINTS
 
 190 CALL flface (ntype,card,z(ielmt+jloc-1),grid)
 IF (error) GO TO 210
 
!     INCLUDE CONNECTIONS IN GRID POINT CONNECTIVITY TABLE
!        1) NUMBER OF STRUCTURE GRID POINTS CONNECTED TO EACH FLUID
!        2) NUMBER OF STRUCTURAL GRID POINTS CONNECTED TO EACH
!           STRUCTURE POINT
 
 ngrdf = 4
 IF (grid(4) < 0) ngrdf = 3
 ngrds = 4
 IF (z(ielmt+jloc+4) < 0) ngrds = 3
 DO  i = 1,ngrdf
   j = grid(i) - 1
   z(igrid+j) = z(igrid+j) + ngrds
 END DO
 DO  i = 1,ngrds
   j = z(ielmt+jloc+i) - 1
   z(igrid+j) = z(igrid+j) - ngrds
 END DO
 
!     WRITE 12 WORD RECORD FOR THIS ENTRY ON FBELM
 
!                  WORD      DESCRIPTION
 
!                  1         FLUID ELEMENT ID
!                  2         STRUCTURAL ELEMENT ID
!                  3-6       STRUCTURE GRID POINTS
!                  7         GRAVITY LOAD ID
!                  8         MATERIAL ID
!                  9-12      FLUID GRID POINTS
 
 CALL WRITE (fbelm,z(ielmt+jloc-1),7,0)
 CALL WRITE (fbelm,card(2),1,0)
 CALL WRITE (fbelm,grid,4,0)
 GO TO 210
 
!     THIS ENTRY DESCRIBES THE FREE SURFACE - FIND THE FLUIDS GRID
!     POINTS WHICH DEFINE THE FACE ID GIVEN
 
 200 CALL flface (ntype,card,z(ielmt+jloc-1),grid)
 IF (error) GO TO 210
 
!     INCLUDE CONNECTIONS IN GRID POINT CONNECTIVITY TABLE
!        1) NUMBER OF FREE SURFACE POINTS CONNECTED TO THIS FREE
!           SURFACE POINT
 
 ngrdf = 4
 IF (grid(4) < 0) ngrdf = 3
 DO  i = 1,ngrdf
   j = grid(i) - 1
   z(igrid+j) = z(igrid+j) + ngrdf*1000000
 END DO
 
!     WRITE 7 WORD RECORD ON FRELM FILE
 
!                  WORD      DESCRIPTION
 
!                  1         FLUID ELEMENT ID
!                  2         MATERIAL FLAG
!                  3-6       FLUID GRID POINTS
!                  7         GRAVITY LOAD ID
 
 z(ielmt+jloc) = card(2)
 CALL WRITE (frelm,z(ielmt+jloc-1),2,0)
 CALL WRITE (frelm,grid,4,0)
 CALL WRITE (frelm,z(ielmt+jloc+5),1,0)
 
!     FLAG THE ELEMENT TABLE ENTRY AS BEEN PROCESSED AND CHECK IF
!     THE NEXT ENTRY IS FOR THE SAME FLUID ELEMENT
 
 210 z(ielmt+jloc) = -2
 IF (jloc+7 >= nelmt .OR. z(ielmt+jloc+6) /= card(1)) GO TO 170
 jloc = jloc + 7
 GO TO 180
 
 220 CALL CLOSE (ect,1)
 CALL CLOSE (fbelm,1)
 CALL CLOSE (frelm,1)
 mcb(1) = fbelm
 mcb(2) = ngrdt
 mcb(3) = nelm
 CALL wrttrl (mcb)
 mcb(1) = frelm
 CALL wrttrl (mcb)
 
!     MAKE ONE FINAL PASS THROUGH ELEMENT TABLE AND VERIFY THAT
!     EVERY FLUID ELEMENT WAS PROCESSED
 
 lelmt = ielmt + nelmt - 1
 DO  i = ielmt,lelmt,7
   IF (z(i+1) == -2) CYCLE
   IF (z(i+1) == -1) GO TO 230
   error = .true.
   WRITE (nout,8003) ufm,z(i)
   CYCLE
   
   230 error = .true.
   WRITE (nout,8004) ufm,z(i)
   
 END DO
 
!     ELEMENT TABLE IS NO LONGER NEEDED SO DELETE IT AND RETURN
 
 icore = ielmt
 RETURN
 
!     ERROR CONDITIONS
 
 1001 n = -1
 GO TO 1100
 1002 n = -2
 GO TO 1100
 1003 n = -3
 GO TO 1100
 1008 n = -8
 1100 CALL mesage (n,FILE,NAME)
 
!     NO FLUID / STRUCTURE BOUNDARY DEFINED.  FATAL ERROR BECAUSE DMAP
!     CANNOT HANDLE THIS CONDITION
 
 1200 error = .true.
 WRITE (nout,8001) ufm
 
 RETURN
 
!     ERROR FORMATS
 
 8001 FORMAT (a23,' 8001. THERE MUST BE A FLUID/STRUCTURE BOUNDARY IN ',  &
     'HYDROELASTIC ANALYSIS.')
 8002 FORMAT (a23,' 8002, ELEMENT ID',i9,' ON A CFLSTR CARD DOES NOT ',  &
     'REFERENCE A VALID 2D STRUCTURAL ELEMENT.')
 8003 FORMAT (a23,' 8003. ELEMENT ID',i9,' ON A CFLSTR CARD DOES NOT ',  &
     'REFERENCE A VALID FLUID ELEMENT.')
 8004 FORMAT (a23,' 8004. ELEMENT ID',i9,' ON A CFFREE CARD DOES NOT ',  &
     'REFERENCE A VALID FLUID ELEMENT.')
     
END SUBROUTINE flbelm
