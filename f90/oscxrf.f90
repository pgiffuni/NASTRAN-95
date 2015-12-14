SUBROUTINE oscxrf (iop,avail)
 
 INTEGER, INTENT(IN OUT)                  :: iop
 INTEGER, INTENT(IN)                      :: avail
 EXTERNAL        lshift,rshift,andf,orf,complf
 INTEGER :: dbent(3),BLOCK(6),lab(6),iout(32),ihd1(32),  &
            ihd2(32),ihd3(32),ihd4(32),ihd5(32), op
 COMMON /output/ ititl(96),ihead(96)
 COMMON /moddmp/ iflg(5)
 COMMON /system/ ksys(65)
 COMMON /lnklst/ i,nvail,iseqn,kind,itype,mask3,mask4,mask5
 COMMON /zzzzzz/ z(1)
 COMMON /xvps  / vps(3)
 EQUIVALENCE     (ksys(2),op), (ksys(9),nlpp), (ksys(12),nline)
 DATA    ihd1  / 7*4H    ,4HCOSM,4HIC /,4H nas,4HTRAN,4H dma,  &
                 4HP co  ,4HMPIL,4HER -,4H dma,4HP cr,4HOSS ,  &
                 4HREFE  ,4HRENC,4HE li,4HSTIN,4HG   ,9*4H    /
 DATA    ihd2  / 32*4H     /
 DATA    ihd3  / 4HMODU,4HLE n,4HAME ,4HDMAP,4H sta,4HTEME,  &
                 4HNT n,4HUMBE,4HRS  ,23*4H     /
 DATA    ihd4  / 4HDATA,4H blo,4HCK  ,4HDMAP,4H sta,4HTEME,  &
                 4HNT n,4HUMBE,4HRS  ,23*4H     /
 DATA    ihd5  / 4HPARA,4HMETE,4HR   ,4HTYPE,4H    ,4HDMAP,  &
                 4H sta,4HTEME,4HNT n,4HUMBE,4HRS  ,21*4H     /
 DATA    lab   / 4HI   ,4HR   ,4HBCD ,4HRDP ,4HCSP ,4HCDP     /
 DATA    pool  / 4HPOOL  /
 DATA    nblank/ 4H      /    ,iout  /32*4H     /
 DATA    nastk / 4H*     /    ,notapp/4HN.a.    /
 
!     RESTRICT OPEN CORE DUE TO LIMITED FIELD SIZE FOR POINTERS
 
 nvail = avail
 IF (nvail > 16350) nvail = 16350
 
!     PROCESS VARAIABLE PARAMETER LIST
 
 mask2 = lshift(1,16) - 1
 mask1 = andf(lshift(1,20)-1,complf(mask2))
 mask3 = lshift(1,14) - 1
 mask4 = lshift(mask3,14)
 mask5 = complf(orf(mask3,mask4))
 nosgn = complf(lshift(1,ksys(40)-1))
 
 DO  i = 1,1600
   z(i) = 0
 END DO
 k = 3
 i = 1
 kind  =-5
 nparam= 1
 20 itype = andf(vps(k+2),mask1)
 itype = rshift(itype,16)
 LEN   = andf(vps(k+2),mask2)
 CALL linkup (*999,vps(k))
 k = k + LEN + 3
 IF (k > vps(2)) GO TO 30
 nparam = nparam + 1
 GO TO 20
 
!     PROCESS NAMES OF MODULES AND DATA BLOCKS
 
 30 pseq = 0
 40 CALL READ (*200,*998,pool,BLOCK,6,0,q)
 iauto = 0
 mi    = rshift(BLOCK(3),16)
 itype = andf(mask2,BLOCK(3))
 iseqn = andf(nosgn,BLOCK(6))
 kind  = 1
 IF (pseq == iseqn .AND. (mi == 3 .OR. mi == 8)) iauto = 1
 IF (iauto == 1) kind = 2
 pseq = iseqn
 CALL linkup (*999,BLOCK(4))
 kind = 3
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 50
   CASE (    3)
     GO TO 70
   CASE (    4)
     GO TO 90
 END SELECT
 
!     PROCESS FUNCTIONAL MODULE IO SECTIONS
 
 50 irlh = 0
 51 irlh = irlh + 1
 CALL READ (*998,*998,pool,ndb,1,0,q)
 DO  j = 1,ndb
   CALL READ (*998,*998,pool,dbent,3,0,q)
   IF (dbent(1) == 0) CYCLE
   CALL linkup (*999,dbent)
 END DO
 kind = 4
 IF (itype == 1 .AND. irlh == 1) GO TO 51
 kind = 5
 CALL READ (*998,*998,pool,ndb, -1,0,q)
 CALL READ (*998,*998,pool,nparm,1,0,q)
 IF (nparm == 0) GO TO 80
 DO  j = 1,nparm
   CALL READ (*998,*998,pool,il,1,0,q)
   IF (il < 0) THEN
     GO TO    60
   END IF
   55 CALL READ (*998,*998,pool,dbent,-il,0,q)
   CYCLE
   60 il = andf(nosgn,il)
   CALL linkup (*999,vps(il-3))
 END DO
 GO TO 80
 70 IF (mi /= 7) GO TO 80
 kind = 5
 CALL READ (*998,*998,pool,il,1,0,q)
 il = andf(mask2,il)
 CALL linkup (*999,vps(il-3))
 80 CALL fwdrec (*998,pool)
 GO TO 40
 90 mi = mi - 7
 IF (mi < 0) mi = 4
 SELECT CASE ( mi )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 170
   CASE (    4)
     GO TO 120
 END SELECT
 100 CALL READ (*998,*998,pool,ndb,1,0,q)
 kind = 5
 IF (iauto == 1) kind = 6
 DO  j = 1,ndb
   CALL READ (*998,*998,pool,dbent,2,0,q)
   il = dbent(1)
   CALL linkup (*999,vps(il-3))
 END DO
 GO TO 80
 120 CALL READ (*998,*40,pool,ndb,1,0,q)
 kind = 3
 130 DO  j = 1,ndb
   CALL READ (*998,*998,pool,dbent,2,0,q)
   IF (dbent(1) == 0) CYCLE
   CALL linkup (*999,dbent)
 END DO
 IF (mi == 4) GO TO 80
 CALL READ (*998,*998,pool,il,1,0,q)
 IF (il > 0) THEN
   GO TO   150
 ELSE
   GO TO   160
 END IF
 150 kind = 5
 CALL linkup (*999,vps(il-3))
 160 IF (mi == 2) GO TO 120
 170 CALL READ (*998,*40,pool,ndb,1,0,q)
 kind = 3
 CALL READ (*998,*998,pool,dbent,3,0,q)
 IF (dbent(1) == 0) GO TO 180
 CALL linkup (*999,dbent)
 180 ndb = ndb - 1
 GO TO 130
 
!     SORT PARAMETER AND MODULE NAMES, 8-BCD WORD SORT
 
 200 nwds = 4*nparam
 CALL sorta8 (0,0,4,1,z(1),nwds)
 ist = nwds + 1
 j   = i - 1 - nwds
 CALL sorta8 (0,0,4,1,z(ist),j)
 nwds = i - 1
 
!     TRAVERSE LINKED LISTS AND GENERATE OUTPUT
 
 k   = 1
 kdh = 0
 DO  j = 1,32
   ihead(j   ) = ihd1(j)
   ihead(j+32) = ihd2(j)
   ihead(j+64) = ihd5(j)
 END DO
 CALL page
 WRITE (op,900)
 nline = nline + 1
 
!     PROCESS PARAMETER NAMES
 
 270 iout(2) = z(k  )
 iout(3) = z(k+1)
 ntype   = rshift(z(k+2),28)
 iout(4) = nblank
 iout(5) = lab(ntype)
 IF (ntype == 0 .OR. ntype > 6) iout(5) = notapp
 iout(6) = nblank
 
!     TRACE THROUGH LINKED LIST
 
 ii = 7
 280 link = andf(mask3,z(k+2))
 310 isn  = andf(mask3,z(link))
 IF (kdh == 0) isn = -isn
 CALL outpak (ii,iout,isn)
 itemp = rshift(z(link),28)
 IF (itemp == 2 .OR. itemp == 4 .OR. itemp == 6) iout(ii+1) = nastk
 link = rshift(andf(z(link),mask4),14)
 IF (link == 0) GO TO 320
 ii = ii + 2
 GO TO 310
 
!     PRINT OUTPUT
 
 320 nline = nline + 1
 IF (nline <= nlpp) GO TO 321
 CALL page
 nline = nline + 1
 WRITE (op,900)
 nline = nline + 1
 321 WRITE (op,902) (iout(ll),ll=2,32)
 DO  ll = 2,32
   iout(ll) = nblank
 END DO
 328 k = k + 4
 IF (k >= ist) GO TO 330
 IF (kdh == 1) GO TO 337
 IF (kdh == 2) GO TO 425
 GO TO 270
 
!     PROCESS MODULE NAMES
 
 330 IF (kdh > 0) GO TO 340
 kdh = 1
 DO  j = 1,32
   ihead(j+64) = ihd3(j)
 END DO
 WRITE (op,910)
 CALL page
 nline = nline + 1
 WRITE (op,900)
 k   = ist
 ist = nwds
 337 IF (rshift(z(k+3),28) >= 3) GO TO 328
 339 iout(2) = z(k  )
 iout(3) = z(k+1)
 iout(4) = nblank
 ii = 5
 GO TO 280
 
!     PROCESS DATA BLOCKS
 
 340 IF (kdh > 1) GO TO 430
 kdh = 2
 DO  j = 1,32
   ihead(j+64) = ihd4(j)
 END DO
 WRITE (op,905)
 CALL page
 nline = nline + 1
 WRITE (op,900)
 k   = 4*nparam + 1
 ist = nwds
 425 IF (rshift(z(k+3),28) >= 3) GO TO 339
 GO TO 328
 430 WRITE (op,906)
 CALL REWIND (pool)
 CALL skpfil (pool,iop)
 CALL fwdrec (*998,pool)
 GO TO 1000
 998 CALL xgpidg (59,0,0,0)
 GO TO 1000
 999 CALL xgpidg (60,0,0,0)
 1000 RETURN
 
 900 FORMAT (1H )
 902 FORMAT (5X,31A4)
 905 FORMAT (//6X,'* DENOTES AUTOMATICALLY GENERATED INSTRUCTIONS',  &
     /8X,'STATEMENT NUMBER REFERS TO DMAP SEQUENCE NUMBER OF ',  &
     'PREVIOUS INSTRUCTION')
 906 FORMAT (//6X,'* DENOTES STATEMENTS IN WHICH THE DATA BLOCK ',  &
     'APPEARSRS AS OUTPUT.')
 910 FORMAT (//6X,'* DENOTES APPEARANCE OF PARAMETER IN AUTOMATICALLY',  &
     ' GENERATED SAVE INSTRUCTION')
     
END SUBROUTINE oscxrf
