SUBROUTINE trht
     
!     TRANSIENT INTEGRATION HEAT TRANSFER MODULE
 
!     INPUTS  CASEXX,USETD,NLFT,DIT,GPTT,KDD,BDD,RDD,PD,TRL (10)
 
!     OUTPUTS  UDVT,PNLD (2)
 
!     SCRATCHES (7)
!     PARAMETERS BETA(R),TABS(R),NORAD(L),RADLIN(L)
 
!     ICR1 IS LLL
!     ICR2 IS ULL
!     ICR5 IS INITIAL CONDITIONS
!     ICR6 IS A MATRIX
!     ICR3,ICR4,ICR7 ARE DECOMP SCRATCH FILES
 
 INTEGER :: casexx,usetd,nlft,dit,gptt,kdd,bdd,rdd,pd,trl,  &
     udvt,pnld,radlin,iscr1,iscr2,iscr3,iscr4,  &
     iz(1),NAME(2),ipnl(7),sysbuf,dit1,pnl1
 COMMON /BLANK / beta,tabs,norad,radlin,sigma
 COMMON /system/ ksystm(65)
 COMMON /zzzzzz/ z(1)
 COMMON /trhtx / ik(7),ib(7),icr1,icr2,icr3,icr4,icr5,isym,icr6, icr7,tim
 COMMON /trdd1 / nlft1,dit1,nlftp1,nout,icount,iloop1,moda1,nz,  &
     icore,iu2,ip4,ipnl,nmodes,nstep,pnl1,ist,more(6)
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(55),iprec),(z(1),iz(1))
 DATA    casexx, usetd,nlft,dit,gptt,kdd,bdd,rdd,pd ,trl/  &
     101   , 102  ,103 ,104,105 ,106,107,108,109,110/,  &
     udvt  , pnld,iscr1 ,iscr2,iscr3,iscr4,iscr5,iscr6,iscr7/  &
     201   , 202 ,301   ,302  ,303  ,304  ,305  ,306  ,307  /
 DATA    NAME  / 4HTRD ,4H    /,   nb   / 8 /
 
!     SET UP FILES
 
 ik(1) = kdd
 CALL rdtrl (ik)
 ib(1) = bdd
 CALL rdtrl (ib)
 icr1 = iscr1
 icr2 = iscr2
 icr3 = iscr3
 icr4 = iscr4
 icr5 = iscr5
 icr6 = iscr6
 icr7 = iscr7
 
!     SET UP NONLINEAR FILES
 
 nlft1 = nlft
 dit1  = dit
 pnl1  = pnld
 IF (ik(1) <= 0) ik(1) = 0
 IF (ib(1) <= 0) ib(1) = 0
 moda1 = -1
 IF (ib(1) /= 0) ik(2) = ib(2)
 
!     OBTAIN PARAMETERS, INITIAL CONDITIONS
 
 CALL trht1a (casexx,usetd,gptt,trl,ngroup)
 
!     ALLOCATE CORE
 
 nz = korsz(z)
 igroup = nz - 3*ngroup + 1
 nv = 4
 IF (nlftp1 /= 0 .OR. norad /= -1) nv = nv + 3
 IF (nz < nv*ik(2)*iprec-nb*sysbuf-3*ngroup) CALL mesage (-8,0,NAME)
 tim = 0.
 DO  i = 1, ngroup
   CALL klock (itim1)
   nstep  = iz(igroup )
   delta  = z(igroup+1)
   igroup = igroup + 3
   
!     FORM  A  MATRIX AND DECOMPOSE
   
   CALL trht1b (3*ngroup,delta)
   CALL klock  (itim3)
   CALL trht1c (ngroup,udvt,pd,rdd,i)
   CALL klock  (itim2)
   CALL tmtogo (itleft)
   IF (i == ngroup) CYCLE
   IF ((itim3-itim1+((itim2-itim3)/nstep)*iz(igroup)) >= itleft) GO TO 30
 END DO
 20 RETURN
 
 30 CALL mesage (45,ngroup-i,NAME)
 GO TO 20
END SUBROUTINE trht
