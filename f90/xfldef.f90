SUBROUTINE xfldef (name1,name2,nofind)
     
!     THE PURPOSE OF THIS ROUTINE IS TO TURN ON ALL OSCAR ENTRY EXECUTE
!     FLAGS NECESSARY TO DEFINE FILE .
 
!                 DESCRIPTION OF ARGUMENTS
!     NAM1,NAM2 = NAME OF FILE TO BE DEFINED.
!     NOFIND    = INDICATES TO CALLING PROGRAM WHETHER OR NOT FILE WAS
!                 FOUND.
 
 
 INTEGER, INTENT(IN)                      :: name1(1)
 INTEGER, INTENT(IN)                      :: name2(1)
 INTEGER, INTENT(IN OUT)                  :: nofind
 EXTERNAL        andf,orf,complf
 INTEGER :: sol,oscar(1),os(5),ospnt,osbot,  &
     fmed(1),fmedtp,fnm(1),fnmtp,fmdmsk,two,op,ptdtp,  &
     ptdbot,ptdic(1),AND,OR,andf,orf,complf,start, reuse,regen
 COMMON /xmdmsk/ nmskcd,nmskfl,nmskrf,fmdmsk(7)
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp,xgpid1(5), noflgs
 COMMON /system/ bs,op,nogo,dum(78),icpflg
 COMMON /xoldpt/ ptdtp,ptdbot,lptdic,nrlfl,seqno
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
!                  ** CONTROL CARD NAMES **  &
                 ndiag,nsol,ndmap,nestm1,nestm2,nexit, &
!                  ** DMAP CARD NAMES **  &
                  nbegin,nend,njump,ncond,nrept,ntime, &
                  nsave, noutpt,nchkpt,npurge,nequiv,  &
                  ncpw,nbpc,nwpc, maskhi,masklo,isgnon, &
                  nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR, newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xgpi5 / iapp,start,iexit(2),sol,subset,iflag,iestim,  &
     icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /xgpi6 / medtp,fnmtp,cnmtp,medpnt,lmed
 COMMON /two   / two(4)
 EQUIVALENCE     (core(1),os(1),loscar),(osprc,os(2)),  &
     (osbot,os(3)),(iospnt,os(4)), (os(5),oscar(1),fnm(1),fmed(1),ptdic(1)),  &
     (medtp,fmedtp),(two(4),reuse)
 DATA    nxchkp/ 4HXCHK/, ifirst / 0 /
 
 AND(i,j) = andf(i,j)
 OR(i,j)  = orf(i,j)
 
 nam1 = name1(1)
 nam2 = name2(1)
 
!     SCAN OPTDIC FOR FILE NAME
 
 regen  = nofind
 nofind = 1
 IF(ptdbot < ptdtp)  GO TO 200
 DO  ii = ptdtp,ptdbot,3
   i = ptdbot + ptdtp - ii
   IF (ptdic(i) == nam1 .AND. ptdic(i+1) == nam2) GO TO 110
 END DO
 GO TO 200
 
!     FILE IS IN PTDIC - SET REUSE FLAG FOR ALL EQUIVALENCED FILES
 
 110 IF (ptdic(i+2) >= 0) GO TO 130
 DO  j = ptdtp,ptdbot,3
   IF (AND(ptdic(j+2),noflgs) == AND(ptdic(i+2),noflgs))  &
       ptdic(j+2) = OR(ptdic(j+2),reuse)
 END DO
 130 ptdic(i+2) = OR(ptdic(i+2),reuse)
 nofind = 0
 GO TO 1000
 
!     FILE NOT IN PTDIC - CHECK FNM TABLE IF RESTART IS MODIFIED AND
!     APPROACH IS NOT DMAP
 
 200 IF (start == icst .OR. iapp == idmapp) GO TO 1000
 IF (regen < 0)  GO TO 1000
 j = fnmtp + 1
 k = fnmtp + fnm(fnmtp)*3 - 2
 DO  i = j,k,3
   IF (nam1 == fnm(i) .AND. nam2 == fnm(i+1)) GO TO 220
 END DO
 GO TO 1000
 
!     FILE IS IN FNM TABLE - CHECK FOR TABLE ERROR
 
 220 IF (fnm(i+2) <= 0)  GO TO 900
 
!     CLEAR ALL THE MASK WORDS
 
 k = fmed(fmedtp+1)
 DO  l = 1, k
   fmdmsk(l) = 0
 END DO
 
!     SET BIT IN FMDMSK FOR FILE REGENERATION
 
 l = ((fnm(i+2)-1)/31) + 1
 k = fnm(i+2) - 31*(l-1) + 1
 fmdmsk(l) = OR(fmdmsk(l),two(k))
 
!     USE FMDMSK AND FMED TABLE TO TURN ON OSCAR EXECUTE FLAGS
 
 k  = fmed(fmedtp+1)
 j1 = fmedtp + 2
 j2 = j1 + fmed(fmedtp)*fmed(fmedtp+1) - k
 INDEX = 0
 ospnt = 1
 DO  j = j1,j2,k
   DO  k1 = 1,k
     jj = j + k1 - 1
     IF (AND(fmed(jj),fmdmsk(k1)) /= 0)  GO TO 330
   END DO
   CYCLE
   
!     NON-ZERO ENTRY FOUND - COMPUTE DMAP SEQUENCE NUMBER FOR FMED ENTRY
   
   330 n = ((j-j1)/k) + 1
   IF (AND(oscar(iospnt+5),nosgn) < n) GO TO 1000
   
!     SET EXECUTINON FLAG FOR ALL OSCAR ENTRIES WITH SAME DMAP SEQ
!     NUMBER
   
   335 IF (AND(oscar(ospnt+5),nosgn) - n < 0.0) THEN
     GO TO   345
   ELSE IF (AND(oscar(ospnt+5),nosgn) - n == 0.0) THEN
     GO TO   340
   ELSE
     GO TO   350
   END IF
   340 IF (oscar(ospnt+5) < 0 .OR. (oscar(ospnt+3) == nxchkp .AND.  &
       icpflg == 0)) GO TO 345
   IF (ifirst == 1) GO TO 342
   ifirst = 1
   CALL page1
   CALL xgpimw (12,0,0,0)
   342 IF (INDEX == 1) GO TO 344
   INDEX = 1
   CALL xgpimw (3,nam1,nam2,0)
   344 CALL xgpimw (4,0,0,oscar(ospnt))
   nofind = -1
   oscar(ospnt+5) = orf(oscar(ospnt+5),isgnon)
   345 IF (ospnt >= osbot) CYCLE
   ospnt = ospnt + oscar(ospnt)
   GO TO 335
   350 CONTINUE
 END DO
 
!     MAKE SURE SOME MODULES WERE TURNED ON
 
 IF (nofind /= -1)  GO TO 900
 
!     NEGATE FNM TABLE ENTRY FOR THIS FILE
 
 fnm(i+2) = -fnm(i+2)
 
!     TURN OFF REUSE FLAGS IN PTDIC
 
 IF (ptdbot <= ptdtp .OR. iflag /= 0) GO TO 1000
 j = complf(reuse)
 DO  i = ptdtp,ptdbot,3
   ptdic(i+2) = andf(j,ptdic(i+2))
 END DO
 GO TO 1000
 
!     D I A G N O S T I C    M E S S A G E S
 
!     MED OR FILE TABLE INCORRECT FOR REGENERATING FILE
 
 900 CALL xgpidg (41,nam1,nam2,fnm(i+2))
 nofind =-1
 nogo   = 2
 
 1000 RETURN
END SUBROUTINE xfldef
