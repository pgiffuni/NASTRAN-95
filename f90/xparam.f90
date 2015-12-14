SUBROUTINE xparam
     
!     THE PURPOSE OF XPARAM IS TO GENERATE THE PARAMETER SECTION OF AN
!     OSCAR ENTRY,AND TO GENERATE THE VPS TABLE.
 
!          ... DESCRIPTION OF PROGRAM VARIABLES ...
!     ITMP   = TEMPORARY STORAGE FOR PARAMETER NAME AND VALUE.
!     IPVAL  = HIGHEST PRIORITY NOMINAL VALUE IN ITMP.
!     IPRVOP = PREVIOUS OPERATOR OR OPERAND RECEIVED FROM DMAP.
!     INDEX  = TABLE CONTAINING ROW INDEXES FOR ISYNTX TABLE.
!     ISYNTX = SYNTAX TABLE USED TO PROCESS DMAP PARAMETER LIST.
!     NVSTBL = NOMINAL VALUE SOURCE TABLE.
!     NOSPNT = POINTER TO PARAMETER COUNT IN PARAMETER SECTION OF OSCAR.
!     IOSPNT = POINTER TO NEXT AVAILABLE WORD IN OSCAR.
!     ENDCRD = END OF CARD FLAG
!     MPLLN  = LENGTH(IN WORDS) OF MPL PARAMETER VALUE
!     ITYPE  = TABLE FOR TRANSLATING NUMBER TYPE CODES TO WORD LENGTH.
!     ENDCRD = FLAG INDICATING END OF CARD SENSED.
 
!     RETURN CODES FROM XSCNDM
 
!        1  DELIMITOR
!        2  BCD
!        3  VALUE
!        4  END OF CARD
!        5  ERROR ENCOUNTERED
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 DIMENSION       itmp(7),INDEX(2,2),isyntx(4,5),nvstbl(4,4),  &
     itype(6),oscar(1),os(5)
 COMMON /system/ bufsz,optape,nogo
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     ndiag,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt, nchkpt,npurge,nequiv,  &
     ncpw,nbpc,nwpc, maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /xgpid / xxgpid(8),modflg
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xgpi3 / pvt(2)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xvps  / vps(2)
 COMMON /autosm/ nwords,savnam(100)
 EQUIVALENCE     (core(1),os(1),loscar),(os(2),osprc),  &
     (os(3),osbot),(os(4),ospnt),(os(5),oscar(1))
 
 DATA   INDEX  / 1,3,2,4/, isyntx / 3*1,8,3*2,7,3*3,5,4*4,4*6/,  &
     nvstbl / 1,1,3,3,1,1,4,4,1,1,4,4,1,2,4,2/, itype  / 1,1,2,2,2,4/,  &
     ic  /4HC   /,  iv/4HV   /,  iy   /4HY   /,  in/4HN   /,  &
     nvps/4HVPS /,  is/4HS   /,  iastk/4H*   /,  &
     NAME/1 /, ival/2/, NONE/1/, impl/2/, idmap/3/, ipvt/4/
 
!     INITIALIZE
 
 OR (i,j) = orf(i,j)
 AND(i,j) = andf(i,j)
 endcrd = 0
 iprvop = islsh
 nospnt = oscar(ospnt) + ospnt
 iospnt = nospnt + 1
 oscar(nospnt) = 0
 mplbot = mpl(mplpnt-7) + mplpnt - 7
 
!     GET FIRST/NEXT TYPE AND MODIFY CODES FROM DMAP,CHECK FOR $
 
 10 newtyp = 0
 isave  = 0
 15 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 600
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 601
   CASE (    4)
     GO TO 410
   CASE (    5)
     GO TO 570
 END SELECT
 20 IF (dmap(dmppnt) == nblank) GO TO 15
 oscar(nospnt) = 1 + oscar(nospnt)
 j = dmap(dmppnt)
 IF (j /= ic .AND. j /= iv .AND. j /= is) GO TO 602
 IF (j == is) isave = 1
 i = 1
 IF (j == ic) i = 2
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 470
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 470
   CASE (    4)
     GO TO 470
   CASE (    5)
     GO TO 570
 END SELECT
 30 k = dmap(dmppnt)
 IF (k /= iy .AND. k /= in) GO TO 470
 j = 1
 IF (k == in .OR. k == is) j = 2
 
!     USE I AND J TO OBTAIN ROW INDEX FOR SYNTAX TABLE.
 
 75 i = INDEX(i,j)
 
!     INITIALIZE IPVAL,AND ITMP WITH MPL DATA
 
 IF (mplpnt >= mplbot) GO TO 580
 DO  k = 1,7
   itmp(k) = 0
 END DO
 itmp(3) = IABS(mpl(mplpnt))
 
!     CONVERT PARAMETER TYPE CODE TO WORD LENGTH
 
 k = itmp(3)
 mplln = itype(k)
 ipval = NONE
 IF (mpl(mplpnt) < 0) GO TO 60
 DO  k = 1,mplln
   mplpnt = mplpnt + 1
   itmp(k+3) = mpl(mplpnt)
 END DO
 ipval  = impl
 60 mplpnt = mplpnt + 1
 IF (newtyp == 1) THEN
    SELECT CASE ( irturn )
     CASE (    1)
       GO TO 620
     CASE (    2)
       GO TO 100
     CASE (    3)
       GO TO 110
     CASE (    4)
       GO TO 120
     CASE (    5)
       GO TO 570
   END SELECT
 END IF
 
!     SCAN DMAP FOR PARAMETER NAME AND VALUE IF ANY, AND CODE DMAP ENTRY
!     FOR USE AS COLUMN INDEX IN SYNTAX TABLE.
 
 70 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 110
   CASE (    4)
     GO TO 120
   CASE (    5)
     GO TO 570
 END SELECT
 90 IF (dmap(dmppnt+1) /= iequl .AND. dmap(dmppnt+1) /= islsh .AND.  &
     dmap(dmppnt+1) /= iastk) GO TO 470
 IF (dmap(dmppnt+1) == iastk) GO TO 70
 j = 2
 IF (dmap(dmppnt+1) == islsh) j = 4
 GO TO 130
 100 j = 1
 
!     CHECK FOR BLANK
 
 IF (dmap(dmppnt) == nblank) GO TO 70
 GO TO 130
 110 j = 3
 GO TO 130
 120 j = 5
 
!     BRANCH ON SYNTAX TABLE VALUE
 
 130 k = isyntx(i,j)
 SELECT CASE ( k )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 180
   CASE (    3)
     GO TO 210
   CASE (    4)
     GO TO 280
   CASE (    5)
     GO TO 200
   CASE (    6)
     GO TO 290
   CASE (    7)
     GO TO 470
   CASE (    8)
     GO TO 190
 END SELECT
 
!     NAME FOUND. NAME TO TEMP,UPDATE PREVOP AND SEARCH PVT FOR VALUE.
 
 140 IF (iprvop == iequl) GO TO 190
 IF (iprvop /= islsh) GO TO 470
 itmp(1) = dmap(dmppnt  )
 itmp(2) = dmap(dmppnt+1)
 iprvop  = NAME
 
!     SCAN PVT
 k = 3
 150 l = andf(pvt(k+2),nosgn)
 l = itype(l)
 IF (dmap(dmppnt) == pvt(k) .AND. dmap(dmppnt+1) == pvt(k+1)) GO TO 160
 k = k + 3 + l
 IF (k-pvt(2) < 0) THEN
   GO TO   150
 ELSE
   GO TO    70
 END IF
 
!     CHECK LENGTH OF PVT VALUE
 
 160 ipval = ipvt
 pvt(k+2) = orf(pvt(k+2),isgnon)
 IF (andf(pvt(k+2),nosgn) /= itmp(3)) GO TO 490
 
!     TRANSFER VALUE TO ITMP
 
 DO  m = 1,l
   j = k + m + 2
   itmp(m+3) = pvt(j)
 END DO
 GO TO 70
 
!     DMAP ENTRY IS = OPERATOR
 
 180 IF (iprvop /= NAME) GO TO 470
 iprvop = iequl
 GO TO 70
 
!     BCD PARAMETER VALUE FOUND
 
 190 IF (itmp(3) /= 3) GO TO 500
 length = 2
 dmppnt = dmppnt - 1
 dmap(dmppnt) = itmp(3)
 GO TO 220
 
!     DMAP ENTRY IS BINARY VALUE
 
 200 IF (iprvop == islsh) GO TO 220
 210 IF (iprvop /= iequl) GO TO 470
 220 iprvop = ival
 IF (ipval == ipvt) GO TO 70
 
!     DMAP VALUE IS HIGHEST PRIORITY
 
 ipval = idmap
 IF (andf(dmap(dmppnt),nosgn) /= itmp(3)) GO TO 500
 
! TRANSFER DMAP VALUE TO ITMP
 
 DO  m = 1,length
   j = dmppnt + m
   itmp(m+3) = dmap(j)
 END DO
 GO TO 70
 
!     DMAP ENTRY IS / OPERATOR
 
 280 IF (iprvop == iequl) GO TO 470
 iprvop = islsh
 GO TO 300
 
!     END OF DMAP INSTRUCTION
 
 290 IF (iprvop == iequl) GO TO 470
 
!     PARAMETER SCANNED,CHECK CORRECTNESS OF NAME AND VALUE AND
!     PROCESS ITMP ACCORDING TO NVSTBL
 
 300 IF (i < 4 .AND. itmp(1) == 0) GO TO 510
 k = nvstbl(i,ipval)
 
 SELECT CASE ( k )
   CASE (    1)
     GO TO 310
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 390
 END SELECT
 
!     VARIABLE PARAMETER,VALUE TO VPS,POINTER TO OSCAR
 
 310 k = 3
 320 IF (itmp(1) == vps(k) .AND. itmp(2) == vps(k+1)) GO TO 330
 k = k + AND(vps(k+2),maskhi) + 3
 IF (k-vps(2) < 0) THEN
   GO TO   320
 ELSE
   GO TO   350
 END IF
 
!     PARAMETER IS ALREADY IN VPS - MAKE SURE TYPES AGREE.
 
 330 l = andf(rshift(vps(k+2),16),15)
 IF (l == 0) GO TO 335
 IF (l /= andf(itmp(3),15)) GO TO 555
 
!     CHECK VALUE MODIFIED FLAG
 
 335 IF (andf(modflg,vps(k+2)) == 0) GO TO 340
 
!     VALUE HAS BEEN MODIFIED FOR RESTART - DO NOT CHANGE.
 
 GO TO 380
 
!     CHECK IF PREVIOUSLY DEFINED
 
 340 IF (vps(k+2) < 0) GO TO 540
 GO TO 360
 
!     NAME NOT IN VPS,MAKE NEW ENTRY
 
 350 k = vps(2) + 1
 vps(2) = k + 2 + mplln
 IF (vps(2)-vps(1) > 0.0) THEN
   GO TO   560
 END IF
 
!     ITMP NAME,LENGTH,FLAG,VALUE TO VPS
 
 360 l = mplln + 3
 DO  m = 1,l
   j = k + m - 1
   vps(j  ) = itmp(m)
 END DO
 vps(k+2) = OR(mplln,lshift(itmp(3),16))
 IF (ipval == idmap) vps(k+2) = OR(vps(k+2),isgnon)
 
!     LOCATION OF VALUE IN VPS TO OSCAR
 
 380 oscar(iospnt) = k + 3
 IF (isave /= 1) GO TO 385
 nwords = nwords + 1
 savnam(nwords) = k+3
 385 CONTINUE
 oscar(iospnt) = OR(oscar(iospnt),isgnon)
 iospnt = iospnt + 1
 GO TO 10
 
!     CONSTANT PARAMETER,VALUE TO OSCAR
 
 390 oscar(iospnt) = mplln
 DO  m = 1,mplln
   j = iospnt + m
   oscar(j) = itmp(m+3)
 END DO
 iospnt = iospnt + mplln + 1
 GO TO 10
 
!     PROCESS ANY INTEGER, REAL, OR COMPLEX CONSTANTS
 
 601 i = 2
 j = 2
 newtyp = 1
 oscar(nospnt) = oscar(nospnt) + 1
 GO TO 75
 
!     PROCESS POSSIBLE DELIMITERS - SLASH AND ASTERISK
 
 600 IF (dmap(dmppnt+1) /= iastk) GO TO 610
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 470
   CASE (    2)
     GO TO 605
   CASE (    3)
     GO TO 470
   CASE (    4)
     GO TO 410
   CASE (    5)
     GO TO 570
 END SELECT
 605 i = 2
 j = 2
 newtyp = 1
 oscar(nospnt) = oscar(nospnt) + 1
 GO TO 75
 
!     PROCESS MPL DEFAULTS IF // IS ENCOUNTERED
 
 610 IF (dmap(dmppnt+1) /= islsh) GO TO 470
 i = 2
 j = 2
 newtyp = 1
 oscar(nospnt) = oscar(nospnt) + 1
 GO TO 75
 
!     USE DEFAULT MPL VALUE FOR PARAMETER
 
 620 IF (ipval == NONE) GO TO 470
 oscar(iospnt) = mplln
 DO  m = 1,mplln
   j = iospnt + m
   oscar(j) = itmp(m+3)
 END DO
 iospnt = iospnt + mplln + 1
 GO TO 10
 
!     PROCESS V,N,NAME PARAMETER TYPES AS /NAME/
 
 602 i = 1
 j = 2
 newtyp = 1
 GO TO 75
 
!     ALL PARAMETERS ON DMAP CARD PROCESSED,PROCESS ANY REMAINING ON
!     MPL
 
 410 IF (mplpnt >= mplbot) GO TO 450
 endcrd = 1
 length = IABS(mpl(mplpnt))
 length = itype(length)
 oscar(nospnt) = 1 + oscar(nospnt)
 IF (mpl(mplpnt) < 0) THEN
   GO TO   530
 ELSE IF (mpl(mplpnt) == 0) THEN
   GO TO   480
 END IF
 420 oscar(iospnt) = length
 DO  m = 1,length
   j = iospnt + m
   mplpnt = mplpnt + 1
   oscar(j) = mpl(mplpnt)
 END DO
 440 mplpnt = mplpnt + 1
 iospnt = iospnt + length + 1
 GO TO 410
 
!     RETURN TO XOSGEN
 
 450 oscar(ospnt) = iospnt - ospnt
 irturn = 1
 460 RETURN
 
!     ERROR MESSAGES -
 
!     DMAP CARD FORMAT ERROR
 
 470 CALL xgpidg (3,ospnt,oscar(nospnt),0)
 GO TO 450
 
!     MPL PARAMETER ERROR
 
 480 CALL xgpidg (4,ospnt,oscar(nospnt),0)
 GO TO 450
 
!     PARA CARD ERROR
 
 490 CALL xgpidg (5,0,itmp(1),itmp(2))
 GO TO 70
 
!     ILLEGAL DMAP PARAMETER VALUE
 
 500 CALL xgpidg (6,ospnt,oscar(nospnt),0)
 GO TO 70
 
!     DMAP PARAMETER NAME MISSING
 
 510 CALL xgpidg (7,ospnt,oscar(nospnt),0)
 GO TO 390
 
!     ILLEGAL PARA CARD
 
 520 CALL xgpidg (8,0,itmp(1),itmp(2))
 IF (i-2 > 0) THEN
   GO TO   390
 ELSE
   GO TO   310
 END IF
 
!     CONSTANT PARAMETER NOT DEFINED
 
 530 CALL xgpidg (9,ospnt,oscar(nospnt),0)
 IF (endcrd == 1) GO TO 440
 GO TO 390
 
!     WARNING - PARAMETER ALREADY HAD VALUE ASSIGNED PREVIOUSLY
 
 540 IF (ipval /= idmap) GO TO 550
 CALL xgpidg (-42,ospnt,itmp(1),itmp(2))
 550 IF (AND(rshift(vps(k+2),16),15) == AND(itmp(3),15)) GO TO 380
 
!     INCONSISTENT LENGTH USED FOR VARIABLE PARAMETER.
 
 555 CALL xgpidg (15,ospnt,itmp(1),itmp(2))
 GO TO 380
 
!     VPS TABLE OVERFLOW
 
 560 CALL xgpidg (14,nvps,nblank,dmpcnt)
 570 nogo   = 2
 irturn = 2
 GO TO 460
 
!     TOO MANY PARAMETERS IN DMAP PARAMETER LIST.
 
 580 CALL xgpidg (18,ospnt,0,0)
 GO TO 450
END SUBROUTINE xparam
