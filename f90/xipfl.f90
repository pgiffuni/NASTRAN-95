SUBROUTINE xipfl
     
!     THE PURPOSE OF XIOFL IS TO GENERATE THE INPUT AND OUTPUT FILE
!     SECTIONS FOR AN OSCAR ENTRY.
 
 EXTERNAL        andf,orf
 INTEGER :: dmpcnt,dmppnt,bcdcnt,dmap,osprc,osbot,ospnt,  &
     oscar(1),os(5),andf,orf
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi, &
!                  ** CONTROL CARD NAMES **  &
 ndiag,nsol,ndmap,nestm1,nestm2,nexit, &
!                  ** DMAP CARD NAMES **  &
 nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt,  &
     nchkpt,npurge,nequiv,ncpw,nbpc,nwpc,  &
     maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /passer/ istopf,modnam,icomon
 EQUIVALENCE     (core(1),os(1),loscar),(os(2),osprc),  &
     (os(3),osbot),(os(4),ospnt),(os(5),oscar(1))
 
 
!     SET INPUT FILE FLAG
 
 iofl = 1
 k3   = 0
 istopf = 0
 j = mpl(mplpnt)
 mplpnt = mplpnt + 1
 IF (j /= 0) GO TO 8
 
!     NO INPUT FILES - MAKE ONE NULL ENTRY IN OSCAR
 
 oscar(ospnt+6) = 1
 oscar(ospnt+7) = 0
 oscar(ospnt+8) = 0
 oscar(ospnt+9) = 0
 oscar(ospnt) = oscar(ospnt) + 4
 GO TO 7
 
 
 ENTRY xopfl
!     ===========
 
!     SET O/P FLAG
 
 iofl = 0
 k3   = 0
 istopf = 0
 j = mpl(mplpnt)
 mplpnt = mplpnt + 1
 IF (j /= 0) GO TO 8
 
!     THERE ARE NO O/P FILES - CHANGE OSCAR ENTRY TYPE CODE TO O FORMAT
 
 oscar(ospnt+2) = orf(2,andf(masklo,oscar(ospnt+2)))
 GO TO 7
 
 
!     SCAN INPUT OR OUTPUT SECTION
 
 8 i = ospnt + oscar(ospnt)
 istopf   = i
 oscar(i) = j
 oscar(ospnt) = 1 + oscar(ospnt)
 i = i + 1
 j = i + 3*(j-1)
 oscar(ospnt) = j + 3 - ospnt
 
!     ZERO I/O SECTION
 
 l = j + 2
 DO  k = i,l
   oscar(k) = 0
 END DO
 
!     ENTER FILE NAME IN OSCAR FROM DMAP
 
 DO  k = i,j,3
   CALL xscndm
   SELECT CASE ( irturn )
     CASE (    1)
       GO TO 30
     CASE (    2)
       GO TO 2
     CASE (    3)
       GO TO 20
     CASE (    4)
       GO TO 30
     CASE (    5)
       GO TO 20
   END SELECT
   
!     OK IF NAME RETURNED FROM XSCNDM
   
   2 IF (dmap(dmppnt) == nblank) CYCLE
   
!     ENTER NAME IN OSCAR AND INITIALIZE ORDNAL
   
   oscar(k  ) = dmap(dmppnt  )
   oscar(k+1) = dmap(dmppnt+1)
   oscar(k+2) = 0
 END DO
 7 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 15
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 20
   CASE (    4)
     GO TO 22
   CASE (    5)
     GO TO 20
 END SELECT
 15 IF (dmap(dmppnt+1) == islsh) GO TO 22
 
!     NORMAL EXIT IF DMAP OPERATOR IS /
 
!     ERROR EXIT
!     BLANK ITEM IN O/P SECTION OF TYPE O FORMAT IS OKAY
 
 20 IF (j == 0 .AND. iofl == 0 .AND. dmap(dmppnt) == nblank) GO TO 7
 k1 = 1 + (k-i)/3
 k2 = 1 + (j-i)/3
 IF (k1  <= k2) GO TO 21
 IF (k3  ==  1) GO TO 7
 IF (iofl == 1) CALL xgpidg (62,ospnt,0,0)
 IF (iofl == 0) CALL xgpidg (63,ospnt,0,0)
 k3 = 1
 GO TO 7
 21 irturn = 2
 GO TO 25
 22 irturn = 1
 25 RETURN
 
 
!     DELIMITER OR END OF INSTRUCTION ENCOUNTERED BEFORE ANTICIPATED -
!     CHECK FOR ILLEGAL INPUT FORMAT
 
 30 IF (iofl /= 1 .OR. dmap(dmppnt+1) /= islsh) GO TO 20
 IF (icomon == 0) GO TO 21
 ityp = andf(oscar(ospnt+2),7)
 IF (ityp == 2) GO TO 22
 
!     FIRST INPUT FILE WAS NULL - SHIFT I/P SECTION BY ONE ENTRY AND
!     ZERO FIRST ENTRY
!     ISSUE WARNING MESSAGE
 
 CALL xgpidg (-1,ospnt,0,0)
 IF (i == j) GO TO 22
 i = i + 3
 j = j + 2
 DO  k = i,j
   l = j - k + i
   oscar(l  ) = oscar(l-3)
 END DO
 oscar(i-3) = 0
 oscar(i-2) = 0
 oscar(i-1) = 0
 GO TO 22
END SUBROUTINE xipfl
