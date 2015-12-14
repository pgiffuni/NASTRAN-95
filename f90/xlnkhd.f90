SUBROUTINE xlnkhd
     
!     THE PURPOSE OF XLNKHD IS TO GENERATE THE LINK HEADER SECTION FOR
!     AN OSCAR ENTRY
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 DIMENSION       med(1),oscar(1),os(5)
 COMMON /system/ isys(81),cpflg
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     ndiag,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt, nchkpt,npurge,nequiv,  &
     ncpw,nbpc,nwpc, maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xgpi5 / iapp,start,alter(2),sol,subset,iflag,iestim,  &
     icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /xgpi6 / medtp,fnmtp,cnmtp,medpnt,lmed,dummy(5),ifirst
 COMMON /xmdmsk/ nmskcd,nmskfl,nmskrf,medmsk(7)
 COMMON /xoldpt/ xx(4),seqno
 COMMON /autohd/ ihead
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp
 EQUIVALENCE     (core(1),os(1),loscar),(os(2),osprc),  &
     (os(3),osbot),(os(4),ospnt), (oscar(1),med(1),os(5))
 DATA    xchk  / 4HXCHK   /
 
 OR (i,j) = orf(i,j)
 AND(i,j) = andf(i,j)
 mpler = mpl(mplpnt+3)
 IF (ihead == 1) mpler = 4
 
!     CHECK FOR DECLARATIVE INSTRUCTION
 
 IF (ihead == 1) GO TO 20
 IF (mpler /= 5) GO TO 10
 ospnt = oscar(osbot) + osbot
 GO TO 20
 
!     UPDATE OSCAR PARAMETERS
 
 10 osprc = osbot
 osbot = oscar(osbot) + osbot
 ospnt = osbot
 iseqn = oscar(osprc+1) + 1
 
!     LOAD LINK HEADER INFORMATION
 
 oscar(ospnt    ) = 6
 oscar(ospnt + 1) = iseqn
 oscar(ospnt + 2) = mpler + lshift(modidx,16)
 oscar(ospnt + 3) = dmap(dmppnt    )
 oscar(ospnt + 4) = dmap(dmppnt + 1)
 oscar(ospnt + 5) = dmpcnt
 
 mplpnt = mplpnt + 4
 20 oscar(ospnt+5) = OR(isgnon,oscar(ospnt+5))
 
!     ALWAYS RAISE EXECUTE FLAG FOR COLD START RUNS
 
 IF (start == icst)  GO TO 70
 
!     COMPARE SEQ NO. WITH REENTRY SEQ NO.
 
 IF (dmpcnt < rshift(seqno,16)) GO TO 30
 
!     WE ARE BEYOND REENTRY POINT - EXECUTE ALL MODULES HERE ON OUT.
 
 IF (andf(maskhi,seqno) == 0 .AND. mpler /= 5)  &
     seqno = OR(iseqn,AND(masklo,seqno))
 GO TO 70
 
!     WE ARE BEFORE REENTRY POINT - CHECK APPROACH AND TYPE OF RESTART
!     ALWAYS RAISE EXECUTE FLAG FOR INSERT FOR MODIFIED RESTARTS.
 
 30 IF (insert /= 0 .AND. start == imst) GO TO 70
 IF (start == imst) GO TO 40
 
!     LOWER EXECUTE FLAG FOR UNMODIFIED RESTART RUNS.
 
 oscar(ospnt+5) = AND(nosgn,oscar(ospnt+5))
 IF (mpler == 5) GO TO 90
 RETURN
 
!     FOR RIGID FORMAT - CHECK DECISION TABLE FOR MODIFIED RESTART
 
 40 i = med(medtp+1)
 DO  j = 1,i
   k = medpnt + j - 1
   IF (AND(med(k),medmsk(j)) /= 0) GO TO 70
 END DO
 oscar(ospnt+5) = AND(nosgn,oscar(ospnt+5))
 70 IF (oscar(ospnt+3) == xchk .AND. cpflg == 0)  &
     oscar(ospnt+5) = AND(nosgn,oscar(ospnt+5))
 IF (oscar(ospnt+5) >= 0 .AND. mpler /= 5) RETURN
 
!     PRINT COMPILE/EXECUTE FLAG FOR RESTART
 
 90 IF (start == icst   .OR.  ifirst == 0) RETURN
 IF (dmpcnt == iflag .AND. insert == 0) RETURN
 iflag = dmpcnt
 i = 7
 IF (mpler == 5) i = 10
 CALL xgpimw (i,0,0,0)
 RETURN
END SUBROUTINE xlnkhd
