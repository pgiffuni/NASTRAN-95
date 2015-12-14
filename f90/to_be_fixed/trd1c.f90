SUBROUTINE trd1c(ic,pd,ngroup,nlftp,udv,iloop,scr1,dit,nlft,noue,
!RLBR SPR94003 9/94
!RLBR1                 MODAL,PNL)  &
 modal,pnl,iskip)
 
!     THIS ROUTINE STEPS INTEGRATION PROCEDURE
 
!     THIS ROUTINE IS SUITABLE FOR SINGLE PRECISION OPERATION
 
 
 INTEGER, INTENT(IN)                      :: ic
 INTEGER, INTENT(IN)                      :: pd
 INTEGER, INTENT(IN)                      :: ngroup
 INTEGER, INTENT(IN)                      :: nlftp
 INTEGER, INTENT(IN)                      :: udv
 INTEGER, INTENT(IN)                      :: iloop
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: dit
 INTEGER, INTENT(IN)                      :: nlft
 INTEGER, INTENT(IN)                      :: noue
 LOGICAL :: nopd
 
 INTEGER :: dit1,pnl1,pnl
 INTEGER :: sysbuf,FILE,iz(1),mcb(7),ipnl(7)
!RLBR SPR94003 9/94
!RLBR INTEGER SUBNAM(2)
 INTEGER :: subnam(2), moutpu(7)
 
 COMMON /BLANK /dummy(4), ncol
!RLBR SPR94003 9/94
!RLBR COMMON /SYSTEM/SYSBUF
 COMMON /system/sysbuf, nnout, isystm(79), icpflg
 COMMON /zzzzzz/z(1)
 COMMON /packx /it1,it2,ii,jj,incr
 COMMON /trdxx /ik(7),idum(14),iscr1,iscr2,iscr3,iscr4,iscr5,iscr6,  &
     iopen,isym,TO,nopd,ispnl
 COMMON /unpakx/it3,iii,jjj,incr1
 COMMON /trdd1 /nlft1,dit1,nlftp1,nout,icount,iloop1,modal1,nz,  &
     icore,iu2,ip4,ipnl,nmodes,nstep,pnl1,ist,iu1, deltat,ifrst
 
 EQUIVALENCE    (z(1),iz(1))
 
 DATA   subnam /4HTRD1,1HC/
!RLBNB SPR94003 9/94
 DATA  ioutpu, iscr9 /203, 309/
!RLBNE
 
! ----------------------------------------------------------------------
 
!     INITIALIZE
 
 nrow  = ik(3)
 it1   = 1
 it2   = 1
 ii    = 1
 jj    = nrow
 incr  = 1
 it3   = 1
 iii   = 1
 jjj   = nrow
 incr1 = 1
 nz    = korsz(z)
 igroup= nz -3*ngroup +1
 ibuf1 = igroup -sysbuf
 ibuf2 = ibuf1 -sysbuf
 ibuf3 = ibuf2 -sysbuf
 ibuf4 = ibuf3-sysbuf
 ibuf5 = ibuf4-sysbuf
 ibuf6 = ibuf5-sysbuf
 ibuf7 = ibuf6-sysbuf
 ibuf8 = ibuf7 -sysbuf
!RLBNB SPR94003 9/94
 IF (nlftp == 0) ibuf8 = ibuf7
 ibuf9 = ibuf8 - sysbuf
 ibufa = ibuf9 - sysbuf
 IF (icpflg == 0) ibufa = ibuf8
 IF (icpflg /= 0 .AND. iskip == 1) ibufa = ibuf9
 nz = ibufa - 1
!RLBNE
!RLBD SPR94003 9/94 NZ = IBUF7-1
!RLBD SPR94003 9/94 IF(NLFTP .NE. 0) NZ = IBUF8-1
 iopen = 0
!RLBR SPR94003 9/94 ICRQ = 14*NROW + 1 - NZ
 icrq = 14*(nrow+1) + 1 - nz
 IF(icrq > 0) GO TO 430
!RLBR SPR94003 9/94 IU1=0
 iu1=1
!RLBR SPR94003 9/94 IU2= IU1+NROW
 iu2= iu1+nrow + 1
!RLBR SPR94003 9/94 IU3= IU2+ NROW
 iu3= iu2+ nrow + 1
!RLBR SPR94003 9/94 IP1= IU3+ NROW
 ip1= iu3+ nrow + 1
!RLBR SPR94003 9/94 IP2= IP1+ NROW
 ip2   = ip1+ nrow
 ip3   = ip2+ nrow
 ip4   = ip3+ nrow
 nlft1 = nlft
 dit1  = dit
 nlftp1= nlftp
 iloop1= iloop
 modal1= modal
 ist   = 0
!RLBR SPR94003 9/94 NZ    = NZ - 14*NROW - 1
 nz    = nz - 14*(nrow+1) - 1
 icore = ip4 +nrow
 nmodes= nrow- noue
 pnl1  = pnl
 ASSIGN 60 TO iret1
 nstep = iz(igroup) + 1
 deltat= z(igroup+1)
 nout  = iz(igroup+2)
 IF( iloop /= 1) GO TO 210
 
!     FIRST ENTRY INITIALIZE STUFF
 
 ist  =-1
 FILE = pd
 
!     PUT P0 IN IP2
 
 ipnt = ip2
 nopd = .true.
 ASSIGN 5 TO iretn
 CALL OPEN(*310,pd,iz(ibuf2),0)
 CALL skprec(pd,1)
 nopd = .false.
 GO TO 290
!RLBD SPR94003 9/94     5 FILE = UDV
!RLBD SPR94003 9/94       IAPEND = 0
!RLBR SPR94003 9/94       IF (NCOL .LE. 0) GO TO 8
 5 IF (ncol > 2) GO TO 325
!RLBD SPR94003 9/94      MCB(1) = UDV
!RLBD SPR94003 9/94      CALL RDTRL (MCB)
!RLBD SPR94003 9/94      IF (MCB(2) .NE. 0) GO TO 330
!RLBR SPR94003 9/94    8 CALL GOPEN (UDV,IZ(IBUF3),1)
 CALL gopen (udv, iz(ibuf3), 1)
 CALL makmcb (mcb,udv,nrow,2,2)
!RLBNB SPR94003 9/94
 8 IF (icpflg == 0) GO TO 10
 CALL makmcb (moutpu, ioutpu, nrow+1, iskip, 2)
 CALL gopen (ioutpu, iz(ibuf9), 1)
 IF (iskip == 0) CALL gopen (iscr9, iz(ibufa), 1)
!RLBNE
 10 IF (nlftp == 0) GO TO 20
 
!     CHECK TO SEE IF PNL HAS BEEN PRE-PURGED.
 
 ipnl(1)= pnl1
 CALL rdtrl(ipnl)
 ispnl= 0
 IF(ipnl(1) <= 0) GO TO 20
 ispnl= 1
 CALL gopen(pnl1,iz(ibuf8),1)
 CALL makmcb(ipnl,pnl1,nrow,2,1)
 20 CONTINUE
!RLBR SPR94003 9/94      IF (IAPEND .EQ. 1) GO TO 50
 IF (ncol > 2) GO TO 50
 FILE = ic
 CALL gopen(ic,iz(ibuf1),0)
 ASSIGN 30 TO iretn
 ipnt = iu2
 GO TO 290
 30 ASSIGN 40 TO iretn
 ipnt = iu3
 GO TO 290
 40 CALL CLOSE(ic,1)
 nstep = iz(igroup)+1
 deltat= z(igroup+1)
 nout  = iz(igroup+2)
 
!     FORM  U=1, PO, P-1
 
 CALL form1( z(iu2+1),z(iu3+1),z(iu1+1),z(ip2+1),z(ip1+1),deltat, z(ibuf1))
 
!     START TIME STEP COUNT
 
 50 CONTINUE
 icount = 1
!RLBNB SPR94003 9/94
 mcol = 1
!RLBNE
 60 CONTINUE
 IF (nlftp == 0) GO TO 62
 ifrst=0
 CALL trd1d
 ifrst=1
 62 CONTINUE
 
!     OPEN FBS FILES
 
 FILE = iscr1
 CALL OPEN(*390,iscr1,iz(ibuf4),0)
 FILE = iscr2
 CALL OPEN(*390,iscr2,iz(ibuf5),0)
 FILE = iscr3
!IBMR 5/95
!     CALL OPEN(*390,ISCR3,IZ(IBUF6),0)
 IF ( isym == 1 ) CALL OPEN(*390,iscr3,iz(ibuf6),0)
 FILE = iscr4
 CALL OPEN(*390,iscr4,iz(ibuf7),0)
 
!     ZERO P*
 
 70 CALL tmtogo(itleft)
 IF(itleft <= 0) GO TO 170
 DO  i = 1,nrow
   k = ip4 +i
   z(k) =0.0
 END DO
 IF(nlftp == 0) GO TO 90
 
!     FORM NON-LINEAR LOADS
 
 CALL trd1d
 IF(icount == 1 .OR. icount == nstep .OR. MOD(icount+ist,nout)  &
      == 0) GO TO 85
 GO TO 90
 85 IF (ispnl > 0) CALL pack (z(ip4+1), pnl, ipnl)
 
!     BRING IN NEXT P
 
 90 ipnt = ip3
 FILE = pd
 ASSIGN 100 TO iretn
 IF ( nopd ) GO TO 310
 GO TO 290
 
!     ADD P-S TO FORM P*
 
 100 DO  i=1,nrow
   k = ip4 + i
   l = ip1 + i
   m = ip2 + i
   j = ip3 + i
   z(k) = z(k) +(z(l) + z(m) + z(j))/3.0
 END DO
 IF (iloop /= 1.OR.icount /= 1) GO TO 115
!RLBR SPR94003 9/94      IF (IAPEND .EQ. 1) GO TO 115
 IF (ncol > 2) GO TO 113
 
!     OUTPUT INITIAL DISPLACEMENT
 
 CALL pack (z(iu2 + 1), udv, mcb(1))
 
!     OUTPUT INITIAL VELOCITY
 
 CALL pack (z(iu3 + 1), udv, mcb(1))
!RLBNB SPR94003 9/94
 113 IF (icpflg == 0) GO TO 115
 IF (iskip == 0) CALL WRITE (iscr9, mcol, 1, 0)
!RLBNE
 
!     SOLVE FOR NEXT SOLUTION
 
 115 CALL step  (z(iu3 + 1), z(iu2 + 1), z(iu1 + 1), z(ip4 + 1), iz(ibuf1))
!RLBNB SPR94003 9/94
 IF (icpflg == 0) GO TO 118
 jj = nrow + 1
 z(ip2) = deltat
 IF (iloop /= 1 .AND. icount == 0) z(ip2) = delta1
 CALL pack (z(ip2), ioutpu, moutpu)
 IF (iskip == 1) GO TO 117
 z(iu2) = mcol + 0.1
 CALL pack (z(iu2), ioutpu, moutpu)
 117 jj = nrow
 118 CONTINUE
!RLBNE
 IF (iloop == 1.AND.icount == 1) GO TO 145
 IF (icount == nstep.OR.MOD(icount+ist, nout) == 0) GO TO 130
 IF (icount == 1) GO TO 130
 
!     ROTATE P POINTERS
 
 120 j  = ip1
 ip1= ip2
 ip2= ip3
 ip3= j
 
!     ROTATE U POINTERS
 
 j  = iu1
 iu1= iu2
 iu2= iu3
 iu3= j
 icount = icount +1
!RLBNB SPR94003 9/94
 mcol = mcol + 1
!RLBNE
 IF(icount-nstep < 0) THEN
   GO TO    70
 ELSE IF (icount-nstep == 0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 
!     IT-S OUTPUT TIME -- LUCKY FELLOW
 
 130 CALL pack( z(iu2+1), udv, mcb(1) )
 
!     COMPUTE U DOT
 
 h = 1.0/(2.0*deltat)
 DO  i=1,nrow
   k = ip4 +i
   l = iu3+i
   m = iu1 + i
   z(k) = (z(l)-z(m))*h
 END DO
 CALL pack( z(ip4+1), udv, mcb(1) )
!RLBNB SPR94003 9/94
 IF (icpflg == 0) GO TO 145
 IF (iskip == 0) CALL WRITE (iscr9, mcol, 1, 0)
!RLBNE
 
!     COMPUTE U DOT DOT
 
 145 h = 1.0/(deltat*deltat)
 DO  i=1,nrow
   k = ip4+i
   l = iu3+i
   m = iu1+i
   j = iu2 +i
   z(k) = (z(l)+z(m)- 2.0*z(j))*h
 END DO
 CALL pack( z(ip4+1), udv, mcb(1) )
 GO TO 120
 
!     END OF 1 GROUP
 
 160 IF(iloop /= ngroup) GO TO 200
 GO TO 70
 170 j = 1
 180 CALL CLOSE(udv,j)
 CALL CLOSE(pd, j)
!RLBNB SPR94003 9/94
 IF (icpflg == 0) GO TO 188
 IF (j /= 1 .OR. iskip == 1) GO TO 186
 CALL CLOSE (iscr9, 1)
 
!     COPY THE SINGLE RECORD IN FILE ISCR9 AS THE
!     LAST RECORD IN FILE IOUTPU
 
 CALL gopen (iscr9, iz(ibufa), 0)
 FILE = iscr9
 183 CALL READ (*410, *184, iscr9, z(iu2+1), nrow, 0, iflag)
 CALL WRITE (ioutpu, z(iu2+1), nrow, 0)
 GO TO 183
 184 CALL WRITE (ioutpu, z(iu2+1), iflag, 1)
 CALL CLOSE (iscr9, 1)
 186 CALL CLOSE (ioutpu, j)
 CALL wrttrl (moutpu)
 188 CONTINUE
!RLBNE
 CALL CLOSE(iscr1,1)
 CALL CLOSE(iscr2,1)
!IBMR 5/95
!     CALL CLOSE(ISCR3,1)
 IF ( isym == 1 ) CALL CLOSE(iscr3,1)
 CALL CLOSE(iscr4,1)
 CALL wrttrl(mcb)
 IF( nlftp == 0) GO TO 190
 IF (ispnl == 0) GO TO 190
 CALL CLOSE(pnl,j)
 CALL wrttrl(ipnl)
 190 RETURN
 
!     MORE GROUPS TO COME SAVE STUFF
 
 200 j = 2
 FILE = scr1
 CALL OPEN(*390,scr1,iz(ibuf1),1)
 CALL WRITE(scr1,z(iu3+1),nrow,1)
 CALL WRITE(scr1,z(iu1+1),nrow,1)
 CALL WRITE(scr1,z(iu2+1),nrow,1)
!RLBR SPR94003 9/94
!RLBR CALL WRITE (SCR1,Z(IP1+1),NROW,1)
 CALL WRITE (scr1,z(ip2+1),nrow,1)
 CALL CLOSE(scr1,1)
 GO TO 180
 
!     CHANGE OF TIME STEP--RESTORE POINTERS ETC
 
 210 igroup = igroup +(iloop-1)*3
 delta1 = z(igroup-2)
 nstep  = iz(igroup)
 deltat = z(igroup+1)
 nout   = iz(igroup+2)
 IF (.NOT.nopd) CALL gopen (pd, iz(ibuf2), 2)
 CALL gopen(udv,iz(ibuf3),3)
 mcb(1)= udv
 CALL rdtrl(mcb)
!RLBNB SPR94003 9/94
 IF (icpflg == 0) GO TO 217
 CALL gopen (ioutpu, iz(ibuf9), 3)
 moutpu(1) = ioutpu
 CALL rdtrl (moutpu)
 217 CONTINUE
!RLBNE
 IF(nlftp == 0) GO TO 220
 IF (ispnl > 0) CALL gopen (pnl1, iz(ibuf8), 3)
 220 CONTINUE
 
!     RESTORE STUFF SAVED
 
 FILE = scr1
 CALL OPEN(*390,scr1,iz(ibuf1),0)
 CALL fread(scr1,z(iu1+1),nrow,1)
 CALL fread(scr1,z(iu3+1),nrow,1)
 CALL fread(scr1,z(iu2+1),nrow,1)
 CALL fread(scr1,z(ip2+1),nrow,1)
 CALL CLOSE(scr1,1)
 
!     COMPUTE U DOT
 
!RLBR SPR94003 9/94      H = 1.0D0/DELTA1
 225 h = 1.0D0/delta1
 DO  i=1,nrow
   k =  ip1 +i
   l = iu2 +i
   m = iu3 +i
   z(k) = (z(l)-z(m))*h
 END DO
 
!     COMPUTE U DOT DOT
 
 h = 1.0/(delta1*delta1)
 DO  i=1,nrow
   k = ip4+ i
   l = iu2+ i
   m = iu3+ i
   j = iu1+ i
   z(k) = (z(l)- 2.0*z(m) +z(j))*h
 END DO
!RLBD SPR94003 9/94   250 CONTINUE
 
!     COMPUTE UI PRIME
 
 h = deltat*deltat/2.0
 DO  i=1,nrow
   k =iu1 +i
   l = iu2 +i
   m = ip1+i
   j = ip4 +i
   z(k) = z(l) -deltat*z(m)+ h*z(j)
 END DO
 
!     COMPUTE U DOT PRIME
 
 DO  i=1,nrow
   k = iu3 + i
   l = ip1+i
   m = ip4 + i
   z(k) = z(l) -deltat*z(m)
 END DO
 
!     COMPUTE PI PRIME
 
 DO  i=1,nrow
   k = ip1+i
   z(k) = 0.0
 END DO
 CALL form2(z(ip4+1),z(iu3+1),z(iu1+1),z(ip1+1),z(ibuf1))
 icount = 0
!RLBR SPR94003 9/94      GO TO IRET1, (60,10)
 GO TO iret1, (60,8)
 
!     INTERNAL ROUTINE TO UNPACK VECTORS
 
 290 CALL unpack(*310,FILE,z(ipnt+1))
!RLBR SPR94003 9/94  300 GO TO IRETN, (5,30,40,100,350,360,370)
 300 GO TO iretn, (5,30,40,100,340,350,360,370,385,387)
!RLBR SPR94003 9/94  310 DO 320 INL = 1,NROW
 310 DO  inl = iii, jjj
   k = ipnt +inl
   z(k) = 0.0
 END DO
 GO TO 300
!RLBNB SPR94003 9/94
!     THE FOLLOWING LINES (UNTIL CRPKNE) REPRESENT
!     REPLACEMENTS FOR THE OLD CODE WHICH HAS BEEN
!     DELETED BELOW
 
!     RETRIEVE REQUIRED INFORMATION FROM
!     THE CHECKPOINT RUN
 
 325 mcol = ncol
 CALL gopen (ioutpu, iz(ibuf4), 0)
 moutpu(1) = ioutpu
 CALL rdtrl (moutpu)
 jskip = 1
 IF (moutpu(4) == 1) GO TO 335
 jskip = 2
 CALL skprec (ioutpu, moutpu(2))
 FILE = ioutpu
 nwds = ncol - 1
 327 CALL READ (*410, *330, ioutpu, mcol, -nwds, 0, iflag)
 GO TO 333
 330 nwds = nwds - iflag
 GO TO 327
 333 CALL READ (*410, *333, ioutpu, mcol, 1, 0, iflag)
 CALL REWIND (ioutpu)
 CALL skprec (ioutpu, 1)
 
 335 CALL skprec (ioutpu, jskip*(mcol-1))
 FILE = ioutpu
 jjj = nrow + 1
 
!     GET P SUB I+1
 
 ipnt = ip2 - 1
 ASSIGN 340 TO iretn
 GO TO 290
 340 itype = 1
 delta1 = z(ip2)
 IF (delta1 == deltat) GO TO 345
 itype = 2
 GO TO 350
 345 CALL skprec (ioutpu, -(jskip+1))
 
!     GET P SUB I
 
 ipnt = ip1 - 1
 ASSIGN 350 TO iretn
 GO TO 290
 350 CALL CLOSE (ioutpu, 1)
 
 FILE = udv
 CALL gopen (udv, iz(ibuf3), 0)
 k = 3*(ncol - 1)
 kk = 5
 kkk = 4
 kkp = 0
 jjj = nrow
 CALL skprec (udv, k)
 
!     GET U SUB I+1
 
 ipnt = iu2
 ASSIGN 360 TO iretn
 GO TO 290
 
!     GET U DOT SUB I+1
 
 360 ipnt = ip3
 ASSIGN 370 TO iretn
 GO TO 290
 
 370 IF (mcol == ncol) GO TO 380
 CALL CLOSE (udv, 1)
 FILE = ioutpu
 CALL gopen (ioutpu, iz(ibuf4), 0)
 k = 2*mcol - 3
 kk = 0
 kkk = 3
 kkp = 1
 jjj = nrow + 1
 CALL skprec (ioutpu, k)
 
!     GET U SUB I
 
 380 ipnt = iu1 - kkp
 IF (itype == 2) ipnt = iu3 - kkp
 CALL skprec (FILE, -kk)
 ASSIGN 385 TO iretn
 GO TO 290
 385 IF (itype == 1) GO TO 388
 IF (mcol == ncol) GO TO 386
 itest = z(ipnt+1)
 IF (mcol == itest+1) GO TO 386
 WRITE (nnout, 500)
 CALL mesage (-61, 0, 0)
 386 CALL skprec (FILE, -kkk)
 
!     GET U SUB I-1
 
 ipnt = iu1 - kkp
 ASSIGN 387 TO iretn
 GO TO 290
 387 IF (mcol == ncol) GO TO 388
 itest = z(ipnt+1)
 IF (mcol == itest+2) GO TO 388
 WRITE (nnout, 600)
 CALL mesage (-61, 0, 0)
 388 CALL CLOSE (FILE, 1)
 jjj = nrow
 CALL gopen (udv, iz(ibuf3), 1)
 CALL makmcb(mcb,udv,nrow,2,1)
 
!     OUTPUT INITIAL DISPLACEMENT
 
 CALL pack (z(iu2 + 1), udv, mcb(1))
 
!     OUTPUT INITIAL VELOCITY
 
 CALL pack (z(ip3 + 1), udv, mcb(1))
 IF (itype == 1) GO TO 8
 ASSIGN 8 TO iret1
 GO TO 225
!RLBNE
!RLBDB SPR94003 9/94
!RLBD C
!RLBD C     RETRIEVE LAST VECTOR
!RLBD C
!RLBD   330 CALL GOPEN(UDV,IZ(IBUF3),0)
!RLBD       K = 3*(NCOL - 1)
!RLBD       IAPEND = 1
!RLBD       CALL SKPREC(UDV,K)
!RLBD C
!RLBD C     GET U SUB I+1
!RLBD C
!RLBD       IPNT = IU2
!RLBD       ASSIGN 350 TO IRETN
!RLBD       GO TO 290
!RLBD CP
!RLBD C     GET U SUB I+1 DOT
!RLBD C
!RLBD   350 IPNT = IP1
!RLBD       ASSIGN 360 TO IRETN
!RLBD       GO TO 290
!RLBD C
!RLBD C     GET U SUB I+1 DOT DOT
!RLBD C
!RLBD   360 IPNT = IP4
!RLBD       ASSIGN 370 TO IRETN
!RLBD       GO TO 290
!RLBD   370 CONTINUE
!RLBD       CALL CLOSE(UDV,1)
!RLBD       CALL GOPEN (UDV, IZ(IBUF3), 1)
!RLBD       CALL MAKMCB (MCB, UDV, NROW, 2, 1)
!RLBD C
!RLBD C     OUTPUT INITIAL DISPLACEMENT
!RLBD C
!RLBD       CALL PACK (Z(IU2+1), UDV, MCB(1))
!RLBD C
!RLBD C     OUTPUT INITIAL VELOCITY
!RLBD C
!RLBD       CALL PACK (Z(IP1+1), UDV, MCB(1))
!RLBD C
!RLBD C     FORM P SUB I+1
!RLBD C
!RLBD       DO 380 I =1,NROW
!RLBD       K = IP2+I
!RLBD       Z(K) = 0.0
!RLBD   380 CONTINUE
!RLBD       CALL FORM2(Z(IP4+1),Z(IP1+1),Z(IU2+1),Z(IP2+1),Z(IBUF1))
!RLBD       ASSIGN 10 TO IRET1
!RLBD       GO TO 250
!RLBDE
 
!     ERROR MESAGES
 
 390 ip1 = -1
 400 CALL mesage(ip1,FILE,subnam)
 RETURN
!RLBNB SPR94003 9/94
 410 ip1 = -2
 GO TO 400
!RLBNE
 430 ip1 = -8
 FILE= icrq
 GO TO 400
!RLBNB SPR94003 9/94
 500 FORMAT ('0*** SYSTEM FATAL MESSAGE, LOGIC ERROR 1 IN ',  &
     'SUBROUTINE TRD1C2 WHILE PROCESSING THE RESTART ', 'INFORMATION')
 600 FORMAT ('0*** SYSTEM FATAL MESSAGE, LOGIC ERROR 2 IN ',  &
     'SUBROUTINE TRD1C2 WHILE PROCESSING THE RESTART ', 'INFORMATION')
!RLBNE
END SUBROUTINE trd1c
