SUBROUTINE itmprt (NAME,item,nz,iopt)
     
!     WILL PRINT SOF ITEM - USING  E15.7,I10, OR ALPHA FORMAT
 
 
 REAL, INTENT(IN)                         :: NAME(2)
 REAL, INTENT(IN)                         :: item
 INTEGER, INTENT(IN)                      :: nz
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER :: sysbuf,otpe,two1,rc
 REAL :: subs(3),itm, lods,loap
 DIMENSION       icore(4)
 CHARACTER (LEN=1) :: ccore(2000)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /machin/ machx
 COMMON /xmssg / ufm,uwm
 COMMON /system/ sysbuf,otpe,inx(6),nlpp,inx1(2),line,inx2(26)
 COMMON /zzzzzz/ core(1)
 COMMON /two   / two1(32)
 COMMON /output/ head1(96),head2(96)
 EQUIVALENCE     (ccore,core)
 EQUIVALENCE     (icore(1),core(1))
 DATA    oparen, cparen,ec,ec1,ec2,intgc,alphc,alphc1,cont,uned,d/  &
     4H(1X , 4H)   ,4H,1P,,4HE13.,4H6   ,4H,i13,4H,9X,,4HA4  ,  &
     4HCONT, 4HINUE,4HD   /
 DATA    BLANK , subs,itm/4H    ,4HSUBS,4HTRUC,4HTURE,4HITEM/
 DATA    eqss  / 4HEQSS/, bgss/4HBGSS/, cstm/4HCSTM /,  &
     plts  / 4HPLTS/, lods/4HLODS/, loap/4HLOAP /
 
 
!     TEST FOR FORMATED TABLE PRINT
 
 IF (iopt /= 2) GO TO 5
 IF (item == eqss) GO TO 2000
 IF (item == bgss) GO TO 2100
 IF (item == cstm) GO TO 2200
 IF (item == plts) GO TO 2300
 IF (item == lods) GO TO 2400
 IF (item == loap) GO TO 2500
 5 CONTINUE
 
!     PERFORM UNFORMATED DUMP OF TABLE
 
 CALL sfetch (NAME,item,1,rc)
 IF (rc /= 1) GO TO 190
 DO  i = 1,96
   head2(i) = BLANK
 END DO
 DO  i = 1,3
   head2( i) = subs(i)
 END DO
 head2( 5) = NAME(1)
 head2( 6) = NAME(2)
 head2( 8) = itm
 head2(10) = item
 CALL page
 head2(12) = cont
 head2(13) = uned
 head2(14) = d
 inum = nz/2 - 1
 ns   = inum + 1
 llen = 0
 core(1) = oparen
 irec = 0
 20 WRITE (otpe,30)irec
 irec = irec + 1
 30 FORMAT ('0GROUP NO.',i4)
 line = line + 2
 IF (line >= nlpp) CALL page
 ix   = inum
 nred = 0
 np   = inum - 1
 iv   = 4
 40 ix   = ix + 1
 iout = 4
 nred = nred + 1
 np   = np + 1
 CALL suread (core(ix),1,flag,rc)
 IF (rc-2 < 0.0) THEN
   GO TO    45
 ELSE IF (rc-2 == 0.0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 45 i  = numtyp(core(ix)) + 1
 IF (i == 1 .AND. iv /= 4) i = iv
 iv = i
 SELECT CASE ( i )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 140
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 120
 END SELECT
 
!     REAL NUMBER  (1)
 
 100 iout = 1
 IF (llen+13 > 132) GO TO 160
 110 core(nred+1) = ec
 core(nred+2) = ec1
 core(nred+3) = ec2
 nred = nred + 2
 111 llen = llen + 13
 GO TO 40
 
!     ALPHA   (2)
 
 120 iout = 2
 IF (llen+6 > 132) GO TO 160
 130 core(nred+1) = alphc
 core(nred+2) = alphc1
 nred = nred + 1
 GO TO 111
 
!     INTEGER  (3)
 
 140 iout = 3
 IF (llen+13 > 132) GO TO 160
 150 icore(nred+1) = intgc
 GO TO 111
 
!     BUFFER FULL - END RECORD   PRINT LINE
 
 160 core(nred+1) = cparen
 IF (nred == 1) WRITE (otpe,161)
 IF (nred == 1) GO TO 162
161 FORMAT ('0END OF GROUP - NULL GROUP')
IF ( machx == 2 .OR. machx == 5 ) WRITE  (otpe,core) (icore(i),i=ns,np)
IF ( machx /= 2 .AND. machx /= 5 ) CALL wrtfmt (icore(ns), np-ns+1, ccore)
162 line = line + 1
IF (line >= nlpp) CALL page
llen = 0
nred = 1
np   = inum
core(inum+1) = core(ix)
ix   = inum + 1
SELECT CASE ( iout )
  CASE (    1)
    GO TO 110
  CASE (    2)
    GO TO 130
  CASE (    3)
    GO TO 150
  CASE (    4)
    GO TO 20
END SELECT

!     END OF ITEM

170 WRITE  (otpe,180)
180 FORMAT ('0END OF ITEM')
190 RETURN

!     PERFORM FORMATED LISTING OF TABLE

!     EQSS TABLE

2000 CALL sfetch (NAME,item,1,rc)
IF (rc /= 1) RETURN
CALL suread (core(1),4,nout,rc)
IF (rc /= 1) GO TO 3000
nsub = icore(3)
CALL suread (core(1),nz,nout,rc)
IF (rc /= 2) GO TO 3000
ist  = 1  + nout
left = nz - nout
DO  i = 1,nsub
  CALL suread (core(ist),left,nout,rc)
  IF (rc /= 2 .AND. rc /= 3) GO TO 3000
  icomp = 1 + 2*(i-1)
  CALL cmiwrt (1,NAME,core(icomp),ist,nout,core,icore)
END DO
CALL suread (core(ist),left,nout,rc)
IF (rc /= 2 .AND. rc /= 3) GO TO 3000
CALL cmiwrt (8,NAME,0,ist,nout,core,icore)
RETURN

!     BGSS TABLE

2100 CALL sfetch (NAME,item,1,rc)
IF (rc /= 1) RETURN
ngrd = 1
CALL sjump (ngrd)
IF (ngrd < 0) GO TO 3000
ist = 1
CALL suread (core(ist),nz,nout,rc)
IF (rc /= 2 .AND. rc /= 3) GO TO 3000
CALL cmiwrt (2,NAME,NAME,ist,nout,core,icore)
RETURN

!     CSTM TABLE

2200 CALL sfetch (NAME,item,1,rc)
IF (rc /= 1) RETURN
ngrd = 1
CALL sjump (ngrd)
IF (ngrd < 0) GO TO 3000
ist = 1
CALL suread (core(ist),nz,nout,rc)
IF (rc /= 2 .OR. rc /= 3) GO TO 3000
CALL cmiwrt (3,NAME,NAME,ist,nout,core,icore)
RETURN

!     PLTS TABLE

2300 CALL sfetch (NAME,item,1,rc)
IF (rc /= 1) RETURN
CALL suread (core(1),3,nout,rc)
IF (rc /= 1) GO TO 3000
ist = 1
CALL suread (core(ist),nz,nout,rc)
IF (rc /= 2 .AND. rc /= 3) GO TO 3000
CALL cmiwrt (4,NAME,NAME,ist,nout,core,icore)
RETURN

!     LODS TABLE

2400 icode = 5

2410 CALL sfetch (NAME,item,1,rc)
IF (rc /= 1) RETURN
CALL suread (core(1),4,nout,rc)
IF (rc /= 1) GO TO 3000
nsub = icore(4)
CALL suread (core(1),nz,nout,rc)
IF (rc /= 2) GO TO 3000
ist  = 1  + nout
left = nz - nout
DO  i = 1,nsub
  CALL suread (core(ist),left,nout,rc)
  IF (rc /= 2 .AND. rc /= 3) GO TO 3000
  icomp = 1 + 2*(i-1)
  CALL cmiwrt (icode,NAME,core(icomp),ist,nout,core,icore)
  icode = 6
END DO
RETURN

!     LOAP TABLE

2500 icode = 7
GO TO 2410

!     INSUFFICIENT CORE OR ILLEGAL ITEM FORMAT - FORCE PHYSICAL DUMP

3000 WRITE  (otpe,3010) uwm,item,NAME
3010 FORMAT (a25,' 6231, INSUFFICIENT CORE AVAILABLE OR ILLEGAL ITEM ',  &
    'FORMAT REQUIRES AN UNFORMATED', /31X,  &
    'DUMP TO BE PERFORM FOR ITEM ',a4,' OF SUBSTRUCTURE ',2A4)
GO TO 5
END SUBROUTINE itmprt
