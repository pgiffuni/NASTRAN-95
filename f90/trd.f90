SUBROUTINE trd
     
!     TRANSIENT RESPONSE MODULE DRIVER
 
!     INPUTS   CASEXX,     TRL,NLFT,DIT,KHH,BHH,MHH,PH
!              CASEXX,     TRL,NLFT,DIT,KDD,BDD,MDD,PD
 
!     OUTPUTS  UDVT,PNLD
!              UHVT,PNLH
 
!     PARAMETERS -- MODAL --BCD--INPUT--MODAL=MODAL IMPLIES MODAL
!                   NOUE  --INT--INPUT--NUMBER OF EXTRA POINTS
!                   NONCUP--INT--INPUT--NONCUP=-1 IMPLIES NONCOUPLED
!                   NCOL  --INT--IN/OUT--APPEND FLAG 0  NO APPEND
!                                                    +  COL NUMBER OF
!                                                        LAST TIME STEP
 
!     SCRATCHES   --
 
 
 
 
 
 
 
 INTEGER :: casexx    ,trl      ,nlft     ,dit      ,  &
     kdd       ,bdd      ,mdd      ,pd       ,  &
     udvt      ,pnld     ,modal(2) ,scr1     ,  &
     scr2      ,scr3     ,scr4     ,scr5     , scr6      ,scr7     ,scr8     ,  &
     sr1       ,sr2      ,sr3      ,sr4      , sr5       ,sr6      ,NAME(2)  ,  &
     iz(1)
 
 COMMON /BLANK/  modal     ,noue     ,noncup   ,ncol
 COMMON/system/ibuf, idummy(53), iprec
 COMMON   /zzzzzz /  z(1)
 COMMON   /trdxx /  ik(7)     ,im(7)    ,ib(7)    ,sr1      ,  &
     sr2       ,sr3      ,sr4      ,sr5      ,  &
     sr6       ,iopen    ,isym     ,TO       , nopd      ,ispnl
 
 EQUIVALENCE        ( z(1)    ,iz(1))
 
 DATA               casexx    ,trl      ,nlft     ,dit      /  &
     101       ,102      ,103      ,104      /  &
     ,kdd      ,bdd      ,mdd      ,pd       /  &
     105       ,106      ,107      ,108      /  &
     ,udvt     ,pnld     ,scr1     ,scr2     /  &
     201       ,202      ,301      ,302      /  &
     ,scr3     ,scr4     ,scr5     ,scr6     /  &
     303       ,304      ,305      ,306      /  &
     ,scr7     ,scr8     ,         moda      /  &
     307       ,308      ,         4HMODA    / ,NAME               /  &
     4HTRD     ,4H       /
 
! ----------------------------------------------------------------------
 
!     INITIALIZE
 
 moda1 = -1
 IF ( moda == modal(1)) moda1 = 1
 
!     BUILD INITIAL CONDITIONS
 
 IF (iprec == 1) CALL trd1a  (casexx, trl, scr1, nlftp, ngroup, moda1)
 IF (iprec == 2) CALL trd1a2 (casexx, trl, scr1, nlftp, ngroup, moda1)
 
!     TEST FOR ZERO APPLIED LOAD
 
 ik(1) = scr1
 CALL rdtrl(ik(1))
 IF ( ik(6) /= 0) GO TO 10
 IF ( nlftp /= 0) GO TO 10
 ik(1) = pd
 ik(6) = 0
 CALL rdtrl(ik)
 IF(  ik(6) /= 0) GO TO 10
 IF (ncol > 0) GO TO 10
 CALL mesage(-46,0,0)
 10 CONTINUE
 
!     ESTIMATE CORE
 
 IF( noncup < 0 .AND. modal(1) == moda .AND. nlftp == 0) GO TO 100
 nz = korsz (z)
 igroup = nz- 3*ngroup +1
 ik(1) = kdd
 CALL rdtrl( ik)
 IF( ik(1) < 0) GO TO 20
 nrow = ik(3)
 GO TO 21
 20 ik(1) = 0
 21 ib(1) = bdd
 CALL rdtrl(ib)
 IF( ib(1) < 0) GO TO 30
 nrow = ib(3)
 GO TO 31
 30 ib(1) = 0
 31 im(1) = mdd
 CALL rdtrl(im)
 IF( im(1) < 0) GO TO 35
 nrow = im(3)
 GO TO 36
 35 im(1) = 0
 36 CONTINUE
 icrq = 8*ibuf + 7*iprec*nrow - igroup
 IF(icrq > 0) CALL mesage(-8,icrq,NAME)
 
!     SET UP COMMON
 
 sr1=scr2
 sr2=scr3
 sr3=scr4
 sr4=scr5
 sr5=scr6
 sr6=scr7
 iskip  = 1
 jgroup = igroup
 DO  i = 1, ngroup
   nskip = iz(jgroup+2)
   IF (nskip == 1) GO TO 40
   iskip = 0
   EXIT
   40 jgroup = jgroup + 3
 END DO
 47 DO   i= 1,ngroup
   CALL klock(itime1)
   nstep = iz(igroup)
   delta =  z(igroup+1)
   igroup= igroup +3
   IF (iprec == 1) CALL initl  (3*ngroup, delta)
   IF (iprec == 2) CALL initl2 (3*ngroup, delta)
   CALL klock(itime3)
   IF (iprec == 1)  &
       CALL trd1c  (scr1, pd, ngroup, nlftp, udvt, i, scr8, dit, nlft,  &
       noue, moda1, pnld, iskip)
   IF (iprec == 2)  &
       CALL trd1c2 (scr1, pd, ngroup, nlftp, udvt, i, scr8, dit, nlft,  &
       noue, moda1, pnld, iskip)
   CALL klock (itime2)
   CALL tmtogo(itleft)
   IF( itleft <= 0) GO TO 60
   IF(  i == ngroup) CYCLE
   
!     COMPUTE TIME TO DO NEXT ITERATION
   
   IF (2*(itime3-itime1 + ((itime2-itime3)/nstep)*iz(igroup)) >=  &
       itleft) GO TO 60
 END DO
 55 ik(1) = udvt
 CALL rdtrl (ik(1))
 ncol = ik(2)/3
 RETURN
 
!     UNCOUPLED MODAL
 
 100 CALL trd1e(mdd,bdd,kdd,pd,udvt,ngroup)
 GO TO 55
 
!     INSUFFICIENT TIME LEFT TO FINISH
 
 60 CONTINUE
 ik(1) =udvt
 CALL rdtrl( ik(1))
 ncol = ik(2)/3
 ik(1) = pd
 CALL rdtrl( ik)
 CALL mesage (45, ik(2)-ncol, NAME)
 IF (ncol == 0) CALL mesage(-37,0,NAME)
 RETURN
END SUBROUTINE trd
