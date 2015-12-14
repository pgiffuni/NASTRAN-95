SUBROUTINE cdetm (method,eed,mdd,bdd,kdd,lama,phid,oceigs,nfound,  &
        scr1,scr2,scr3,scr4,scr5,scr6,scr7,scr8)
     
!     SOLVES COMPLEX EIGENVALUE PROBLEM BY DETERMINANT METHOD
 
 
 INTEGER, INTENT(IN)                      :: method
 INTEGER, INTENT(IN)                      :: eed
 INTEGER, INTENT(IN)                      :: mdd
 INTEGER, INTENT(IN)                      :: bdd
 INTEGER, INTENT(IN)                      :: kdd
 INTEGER, INTENT(IN)                      :: lama
 INTEGER, INTENT(IN)                      :: phid
 INTEGER, INTENT(IN)                      :: oceigs
 INTEGER, INTENT(OUT)                     :: nfound
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER, INTENT(IN)                      :: scr5
 INTEGER, INTENT(IN)                      :: scr6
 INTEGER, INTENT(IN)                      :: scr7
 INTEGER, INTENT(IN OUT)                  :: scr8
 INTEGER :: sysbuf,iz(1),fa,fl,fu,sr1,sr2,sr3,powr,  &
     ipole(2),eigc(2),ihead(50),poin,dretn,ips(3),  &
     ipdet(3),FILE,NAME(2), otpe,ilusp(2), iscr2(7),isp(2)
 REAL :: lu
 DOUBLE PRECISION :: zd(1),d1,d2,d3,d4,d5,d6,minda,rl,pr,dt1,pi,dr,di,  &
     psr(3),psi(3),dsr(3),dsi(3),pkr(3),pki(3),  &
     detr(3),deti(3),dt2,dt3,d7,d8,ptr,pti,d9,d10,  &
     hk1r,hk1i,hkr,hki,lamdkr,lamdki,deltkr,deltki,  &
     gkr,gki,rootr,rooti,lamk1r,lamk1i,hkp1r,hkp1i,  &
     h1bar,h2bar,h3bar,test,amcb(2),bmcb(2),cmcb(2),  &
     dpi,d2pi,raddeg,degrad,d4pisq, ddistx,ddist2,zz,zdkm1,zdk
 COMMON /machin/  mach
 COMMON /system/  ksystm(65)
 COMMON /zzzzzz/  z(1)
 COMMON /condad/  dpi,d2pi,raddeg,degrad,d4pisq
 COMMON /output/  head(1)
 COMMON /msgx  /  nmsgx,maxgx
 COMMON /saddx /  nomat,lcadd,mcba(12),mcbb(12),mcbc(12),mcbd(12),  &
     mcbe(12),mc(7)
 COMMON /cdcmpx/  fa(7),fl(7),fu(7),sr1,sr2,sr3,dr,di,powr,nx, minda,ib
 EQUIVALENCE      (ksystm(1),sysbuf), (ksystm(2),otpe),  &
     (amcb(1),mcba(9)), (bmcb(1),mcbb(9)), (cmcb(1),mcbc(9)), (z(1),iz(1),zd(1))
 DATA    ipole ,  eigc,  ihead,        poin,   NAME            /  &
     257,4 ,  207,2, 0,1009,1,47*0, 4HPOIN, 4HCDET, 4HM     /
 DATA    nit,im1, SIGN, numint, iz2,iz3, iz4,iz5, iz6,iz7, iz8 /  &
     20, 1,   -1.0, 4,      2,  3,   4,  5,   6,  7,   8   /
 
!     DEFINITION OF VARIABLES
 
!     METHOD   SELECT SET OF POLES AND EIGC CARDS
!     EED      EIGENVALUE EXTRACTION DATA BLOCK
!     MDD      MASS MATRIX      - DIRECT OR MODAL
!     BDD      DAMPING MATRIX   - DIRECT OR MODAL
!     KDD      STIFFNESS MATRIX - DIRECT OR MODAL
!     LAMA     EIGENVALUE FILE
!     PHID     EIGENVECTOR FILE
!     OEIGS    EIGENVALUE SUMMARY FILE
!     NFOUND   TOTAL NUMBER OF EIGENVALUES FOUND IN ALL REGIONS
!     SCR      SCRATCHES
!     IPOLE    LOCATE WORDS FOR POLES
!     EIGC     LOCATE WORDS FOR EIGC CARDS
!     EPSI     CONVERGENCE CRITERION
!     IBUF     POINTER TO BUFFER
!     NPOLE    NUMBER OF POLES TO BE USED
!     IPOLES   POINTER TO START OF POLES - 4 WORDS PER POLE ID,X,Y,MUL
!     NREGN    NUMBER OF REGIONS
!     IREG1    POINTER TO WORDS DESCRIBING REGIONS
!     INORM    NORMALIZATION METHOD - 1 = MAX,  0 = POINT
!     ISIL     POINTER FOR NORM IF  NORM = 0
!     NE       ESTIMATED NUMBER OF ROOTS IN REGION
!     ND       DESIRED NUMBER OF ROOTS
!     LREGN    LENGTH OF BLOCK DESCRIBING REGION
!     NPASS    NUMBER OF PASSES THROUGH STARTING POINTS
!     NCHANG   NUMBER OF CHANGES OF CONVERGENCE CRITERIA
!     NMOVES   NUMBER OF STARTING POINT MOVES
!     NDCOMP   NUMBER OF DECOMPOSITIONS
!     NFAIL    NUMBER OF FAILURES TO INTERATE TO A ROOT
!     NOUTSD   NUMBER OF PREDICTIONS OUTSIDE REGION
!     ITERM    REASON FOR TERMINATION  0 - FOUND REQUESTED NUMBER
!                                      1 - NO MORE IN REGIONS
!     IRGP     POINTER FOR CURRENT REGION
!     ICNT     NUMBER OF INTERATIONS
!     NIT      MAXIMUM NUMBER OF INTERATIONS/ROOT
!     NUMINT   MAXIMUM NUMBER OF CONVERGENCE  CHANGES
!     NROW     ORDER OF PROBLEM
!     ICMPX    SWITCH IF ROOTS FOUND ARE COMPLEX CONJUGATE -0  NOT-1
!     ISPNT    POINTER TO CURRENT 3 STARTING POINTS
!     PS       SORTED 3 STARTING POINTS
!     DS       SORTED 3 DET OF STARTING POINTS
!     IPS      POWERS OF STARTING POINTS
!     P        TRIAL EIGENVALUE
!     D        SCALED SWEPT DETERMINANT AT P
!     IFPASS   FLAG TO SIGNAL ROOT FROUND ON PASS 1, 0 IF NOT
!     IPOINT   NUMBER OF STARTING POINTS USED IN CURRENT REGION
!     ILUSP    INDEX TO LAST USED STARTING POINT (1ST OF 3) IN EACH
!              SUBRGN
!     NSBRGN   NUMBER OF SUBREGIONS IN PROBLEM THIS REGION
!     NSBDON   FLAG MARKING COMPLETED SUBREGION
!     ISP      POINTS NEAREST AND NEXT NEAREST THE ORIGIN
 
!     STRUCTURE OF REGION
 
!     A1,B1,A2,B2,XL,NE,ND,NF,POINTER TO NEXT REGION (ZERO IF LAST),RL,
!     X  (12  WORDS)
 
!     STARTING POINTS  8NE + 8 WORDS  - D.P. COMPLEX
!     DETERMINANTS     8NE + 8 WORDS    D.P. COMPLEX
!     SCALE FACTORS    4NE + 4 WORDS    2 INTEGERS PER STARTING POINT
 
!     ROOTS   4ND WORDS                 D.P. COMPLEX
 
 
!     DEFINE EPSI (CONVERGENCE CRITERION)
 
 epsi = 1.0E-16
 IF (mach == 5 .OR. mach == 21) epsi = 1.0E-12
 
 lc   = korsz(z)
 ibuf = lc - sysbuf - 1
 lc   = (ibuf/2)*2 - 1
 nosing = 1
 CALL sswtch (7,iprt)
 ising = 0
 fa(1) = kdd
 CALL rdtrl (fa(1))
 IF (fa(1) <= 0) GO TO 1
 idd = kdd
 GO TO 9
 1 fa(1) = mdd
 CALL rdtrl (fa(1))
 IF (fa(1) <= 0) GO TO 2
 idd = mdd
 GO TO 9
 2 fa(1) = bdd
 CALL rdtrl (fa(1))
 IF (fa(1) <= 0) GO TO 990
 idd   = bdd
 9 fa(1) =-scr2
 fa(5) = 4
 fl(1) = idd
 CALL rdtrl (fl(1))
 fl(4) = 4
 fl(5) = 4
 fu(1) = idd
 CALL rdtrl (fu(1))
 fu(4) = 5
 fu(5) = 4
 sr1   = scr3
 sr2   = scr4
 sr3   = scr5
 fl(1) = scr6
 fu(1) = scr7
 DO  i = 1,7
   mcba(i) = 0
   mcbb(i) = 0
   mcbc(i) = 0
   mc(i)   = 0
 END DO
 mcba(1) = kdd
 mcbb(1) = bdd
 mcbc(1) = mdd
 CALL rdtrl (mcba(1))
 CALL rdtrl (mcbb(1))
 CALL rdtrl (mcbc(1))
 
!     MUST HAVE  B OR M MATRICES
 
 IF (mcbb(1) < 0) mcbb(1) = 0
 IF (mcbc(1) < 0) mcbc(1) = 0
 IF (mcbb(1)+mcbc(1) == 0) GO TO 990
 nrow  = MAX0(mcba(3),mcbb(3),mcbc(3))
 icmpx = 0
 IF (mcba(5) > 2 .OR. mcbb(5) > 2 .OR. mcbc(5) > 2) icmpx = 1
 amcb(1) = 1.0D0
 amcb(2) = 0.d0
 mc(2)   = mcba(2)
 mc(3)   = mcba(3)
 mc(4)   = mcba(4)
 mc(5)   = 4
 mcba(8) = 4
 mcbc(8) = 4
 mcbb(8) = 4
 nomat   = 3
 mc(1)   = scr2
 ndesrd  = 0
 
!     PICK UP AND STORE ANY POLES
 
 FILE = eed
 CALL preloc (*950,iz(ibuf),eed)
 npole = 0
 CALL locate (*40,iz(ibuf),ipole(1),iflag)
 
!     FOUND POLE CARDS
 
 20 lc = lc - 4
 30 CALL READ (*970,*40,eed,iz(lc),4,0,iflag)
 IF (iz(lc) /= method) GO TO 30
 npole = npole + 1
 GO TO 20
 40 ipoles = lc + 4
 
!     STORE REGIONS
 
 nregn = 0
 CALL locate (*990,iz(ibuf),eigc(1),iflag)
 50 CALL READ (*970,*990,eed,iz(1),10,0,iflag)
 IF (method == iz(1) .OR. method == -1) GO TO 70
 
!     SKIP REMAINDER OF EIGC CARD
 
 60 CALL READ (*970,*990,eed,iz(1),7,0,iflag)
 IF (iz(iz6) /= -1) GO TO 60
 GO TO 50
 
!     EIGC CARD FOUND - ALLOCATE CORE + BUILD UP REGIONS
 
 70 inorm = 0
 IF (iz(iz4) /= poin) inorm = 1
 isil = iz(iz6)
 IF (z(iz8) /= 0.0) epsi = z(iz8)
 
!     PROCESS EACH REGION DEFINITION
 
 80 CALL READ (*970,*980,eed,z(1),7,0,iflag)
 IF (iz(iz7) < 0) GO TO 130
 nregn = nregn + 1
 alph1 = z(  1)
 w1    = z(iz2)
 alph2 = z(iz3)
 w2    = z(iz4)
 xl    = z(iz5)
 NE    = iz(iz6)
 nd    = iz(iz7)
 IF (nd == 0) nd = 3*NE
 lregn = 20*NE + 4*nd + 32
 ndesrd = ndesrd + nd
 IF (nregn == 1) GO TO 90
 iz(lc+8) = lc - lregn
 90 lc = lc - lregn
 IF (lc    <= 0) GO TO 1000
 IF (nregn /= 1) GO TO 100
 ireg1 = lc
 
!     ZERO REGION
 
 100 k = lc - 1
 DO  i = 1,lregn
   k = k + 1
   IF (i == 9) CYCLE
   iz(k) = 0
 END DO
 
!     STORE CONSTANTS
 z(lc  ) = alph1
 z(lc+1) = w1
 z(lc+2) = alph2
 z(lc+3) = w2
 z(lc+4) = xl
 iz(lc+5)= NE
 iz(lc+6)= nd
 
!     DISTRIBUTE  STARTING POINTS
 
 d1 = alph2 - alph1
 d2 = w2 - w1
 rl = DSQRT(d1*d1+d2*d2)
 z(lc+9) = rl
 d1 = d1/rl
 d2 = d2/rl
 j  = (lc+1)/2 + 6
 d3 = rl/FLOAT(4*NE +4)
 zd(j  ) = d1*d3 + alph1
 zd(j+1) = d2*d3 + w1
 k  = 2*NE + 1
 d3 = rl/FLOAT(k+1)
 d1 = d1*d3
 d2 = d2*d3
 DO  i = 1,k
   j  = j + 2
   zd(j  ) = zd(j-2) + d1
   zd(j+1) = zd(j-1) + d2
 END DO
 GO TO 80
 130 lcadd = lc - 1
 nx = lcadd
 iz(lc+8) = 0
 CALL CLOSE (eed,1)
 IF (lc-4*nrow <= 0) GO TO 1000
 
!     INITIALIZE CUMULATIVE POINTERS
 ifail  = 0
 nfound = 0
 npass  =-1
 nchang = 0
 nmoves = 0
 ndcomp = 0
 nfail  = 0
 noutsd = 0
 iterm  = 1
 ifpass = 1
 
!     RETURN HERE TO SEARCH ALL REGIONS AGAIN
 
 140 npass = npass + 1
 IF (ifpass == 0) GO TO 690
 ifpass = 0
 irgp   = ireg1
 
!     FIND REGION WHICH LACKS ROOTS
 
 DO  i = 1,nregn
   IF (iz(irgp+6) > iz(irgp+7)) GO TO 170
   irgp = iz(irgp+8)
 END DO
 
!     ALL REGIONS HAVE ENOUGH ROOTS - EXIT
 
 GO TO 710
 
!     PICK UP REGION POINTERS AND PARAMETERS
 
 170 alph1 =  z(irgp  )
 w1    =  z(irgp+1)
 alph2 =  z(irgp+2)
 w2    =  z(irgp+3)
 xl    =  z(irgp+4)
 NE    = iz(irgp+5)
 nd    = iz(irgp+6)
 nf    = iz(irgp+7)
 rl    = z(irgp+9)
 xvr   = (alph2-alph1)/rl
 yvr   = (w2-w1)/rl
 ipoint= 0
 ispnt = (irgp+1)/2 + 4
 
!     FIND POINTS CLOSEST TO AND NEXT CLOSEST TO ORIGIN THUS DIVIDING
!     REGION INTO TWO SUBREGIONS
 
 ispt1 = (irgp+13)/2
 lspt  =  ispt1 + 2*(2*NE+2) - 2
 ddistx= 0.
 nxorg = 0
 nrorg = ispt1
 ddist2= zd(ispt1)*zd(ispt1) + zd(ispt1+1)*zd(ispt1+1)
 ispt2 = ispt1 + 2
 DO  i = ispt2,lspt,2
   zz = zd(i)*zd(i) + zd(i+1)*zd(i+1)
   IF (zz > ddist2) EXIT
   nxorg  = nrorg
   nrorg  = i
   ddistx = ddist2
   ddist2 = zz
 END DO
 175 IF (zz > ddistx) GO TO 178
 177 nxorg  = i
 ddistx = zz
 GO TO 179
 178 IF (nxorg == 0) THEN
   GO TO   179
 ELSE
   GO TO   177
 END IF
 
!     CALCULATE THE NUMBER OR SUBREGIONS, NSBRGN. THERE MUST BE AT LEAST
!     3 POINTS EACH SIDE OF BISECTOR IN ORDER TO HAVE 2 SUBREGIONS.
!          ISPT2+2 .LE. (NRORG+NXORG)/2 .LE. LSPT-4
 
 179 CONTINUE
 IF (2*(ispt2+2) <= nrorg+nxorg .AND. nrorg+nxorg <= (lspt-4)*2) GO TO 185
 
!     ONLY ONE SUBREGION
!     FIND FIRST UNEVALUATED POINT
 
 nsbrgn = 1
 k = irgp + 16*NE + 32
 l = 2*NE
 DO  j = 1,l
   IF (iz(k) == 0) GO TO 196
   k = k + 2
   ipoint = ipoint + 1
   ispnt  = ispnt  + 2
 END DO
 
!     ALL TRIED  GO TO BEGINNING
 
 ipoint = 0
 ispnt  = (irgp+1)/2 + 4
 GO TO 196
 
!     TWO SUBREGIONS EXIST. DETERMINE STARTING POINTS FOR EACH
 
 185 nsbrgn = 2
 nsbdon = 0
 isp(1) = nrorg
 IF (nrorg < nxorg) isp(1) = nxorg
 isp(2) = nrorg + nxorg - isp(1)
 kreg   = 2
 ilusp(1) = isp(1)
 ilusp(2) = isp(2) - 2
 GO TO 192
 
!     RETURN HERE TO GET NEW STARTING POINT (OR NEW REGION IF NECESSARY)
!     DETINES  ISPNT
 
 190 IF (nsbrgn == 1) GO TO 196
 IF (nsbdon-1 < 0) THEN
   GO TO   192
 ELSE IF (nsbdon-1 == 0) THEN
   GO TO   195
 ELSE
   GO TO  1945
 END IF
 
!     CHANGE SUBREGIONS
 
 192 kreg = 3 - kreg
 SELECT CASE ( kreg )
   CASE (    1)
     GO TO 194
   CASE (    2)
     GO TO 195
 END SELECT
 
!     PROCESS FIRST SUBREGION
 
 194 ispnt = ilusp(1)
 ls    = isp  (2)
 aloc1 = zd(ls  )
 wloc1 = zd(ls+1)
 ilusp(1) = ilusp(1) + 2
 IF (ispnt+4 == lspt) GO TO 1942
 pr = .45*zd(ispnt+4) + .55*zd(ispnt+6) - aloc1
 pi = .45*zd(ispnt+5) + .55*zd(ispnt+7) - wloc1
 1940 ipoint = (ispnt -ispt1)/2 + 1
 GO TO 220
 
!     PROCESS LAST SET OF STARTING IN FIRST SUBREGION
 
 1942 nsbdon = nsbdon + 1
 pr = .45*zd(ispnt+4) + .55*alph1 - aloc1
 pi = .45*zd(ispnt+5) + .55*w1    - wloc1
 GO TO 1940
 
!     SUBREGION 2 IS COMPLETE.  IS SUBREGION 1 FINISHED AS WELL
 
 1945 IF (nsbdon == 3) GO TO 680
 GO TO 194
 
!     PROCESS SUBREGION 2
 
 195 ispnt = ilusp(2) - 2
 ilusp(2) = ilusp(2) - 2
 ls    = isp(1)
 aloc1 = zd(ls)
 wloc1 = zd(ls+1)
 IF (ispnt == ispt1) GO TO 1952
 pr = -(.45*zd(ispnt-2) + .55*zd(ispnt  )) + aloc1
 pi = -(.45*zd(ispnt-1) + .55*zd(ispnt+1)) + wloc1
 GO TO 1940
 
!     LAST SET OF STARTING POINTS IN SUBREGION2 TO PROCESS
 
 1952 nsbdon = nsbdon + 2
 pr = -(.45*zd(ispt1  ) + .55*zd(ispnt  )) + aloc1
 pi = -(.45*zd(ispt1+1) + .55*zd(ispnt+1)) + wloc1
 GO TO 1940
 
!     ONLY ONE SUBREGION PROCESS FROM END TO END
 
 196 ipoint = ipoint + 1
 ispnt  = ispnt+2
 IF (ipoint > 2*NE) GO TO 680
 
!     FIND OUT IF ANY DETERMINT EVALUATIONS ARE NECESSARY
 
!     COMPUTE LOCAL SEARCH REGION DESCRITIONS
 
 aloc1 = alph1
 wloc1 = w1
 IF (ipoint == 2*NE) GO TO 210
 pr = .45*zd(ispnt+4) + .55*zd(ispnt+6) - aloc1
 pi = .45*zd(ispnt+5) + .55*zd(ispnt+7) - wloc1
 GO TO 220
 210 pr = alph2 - aloc1
 pi = w2 - wloc1
 220 rll= DSQRT(pr*pr+pi*pi)
 k  = irgp + 16*NE + 24 + 2*ipoint
 i  = 1
 ising = 0
 230 k  = k + 2
 IF (iz(k) /= 0) GO TO 250
 
!     EVALUATE DETERMINANT
 
 j  = ispnt + 2*i - 2
 pr = zd(j  )
 pi = zd(j+1)
 ASSIGN 240 TO dretn
 GO TO 810
 240 iz(k  ) = 1
 iz(k+1) = powr
 m = 4*NE + 2 + ispnt + 2*i
 zd(m  ) = dr
 zd(m+1) = di
 250 i = i + 1
 IF (i <= 3) GO TO 230
 IF (ising == 3 .AND. npass == 0) GO TO 701
 
!     SORT STARTING POINTS BY MAGNITUDE OF DET
 
 260 CALL klock (itime1)
 k = ispnt + 4*NE + 4
 l = irgp + 16*NE + 26 + 2*ipoint
 CALL cdetm2 (zd(ispnt),zd(k),iz(l),psr(1),psi(1),dsr(1),dsi(1), ips(1))
 
!     LOAD STARTING POINTS INTO TRAIL EIGENVALUES
 
 DO  i = 1,3
   pkr(i)  = psr(i)
   pki(i)  = psi(i)
   detr(i) = dsr(i)
   deti(i) = dsi(i)
   ipdet(i)= ips(i)
 END DO
 dt2 = 1.0D38
 
!     START INTERATION LOOP
 
 icnt = 1
 280 hk1r = pkr(2) - pkr(1)
 hk1i = pki(2) - pki(1)
 hkr  = pkr(3) - pkr(2)
 hki  = pki(3) - pki(2)
 IF (hkr == 0.0D0 .AND. hki == 0.0D0) GO TO 550
 d1     = hk1r*hk1r + hk1i*hk1i
 lamdkr = (hkr*hk1r + hki*hk1i)/d1
 lamdki = (hki*hk1r - hkr*hk1i)/d1
 deltkr = 1.0D0 + lamdkr
 deltki = lamdki
 
!     COMPUTE GK
 
 d1 = lamdkr*lamdkr - lamdki*lamdki
 d2 = 2.0*lamdkr*lamdki
 d3 = d1*detr(1) - d2*deti(1)
 d4 = d2*detr(1) + d1*deti(1)
 d1 = deltkr*deltkr - deltki*deltki
 d2 = 2.0*deltkr*deltki
 d5 =-d1*detr(2) + d2*deti(2)
 d6 =-d2*detr(2) - d1*deti(2)
 CALL csumm (d3,d4,ipdet(1),d5,d6,ipdet(2),d1,d2,id1)
 d3 = lamdkr + deltkr
 d4 = lamdki + deltki
 d5 = d3*detr(3) - d4*deti(3)
 d6 = d4*detr(3) + d3*deti(3)
 CALL csumm (d1,d2,id1,d5,d6,ipdet(3),gkr,gki,igk)
 
!     COMPUTE TERM UNDER RADICAL IN EQ. 11
 
 d1 = detr(1)*lamdkr - deti(1)*lamdki
 d2 = deti(1)*lamdkr + detr(1)*lamdki
 d3 =-detr(2)*deltkr + deti(2)*deltki
 d4 =-deti(2)*deltkr - detr(2)*deltki
 CALL csumm (d1,d2,ipdet(1),d3,d4,ipdet(2),d5,d6,id1)
 CALL csumm (d5,d6,id1,detr(3),deti(3),ipdet(3),d1,d2,id2)
 d3 = deltkr*lamdkr - deltki *lamdki
 d4 = deltki*lamdkr + deltkr *lamdki
 d5 = d1*d3 - d2*d4
 d6 = d2*d3 + d1*d4
 d1 =-4.0*(detr(3)*d5 - deti(3)*d6)
 d2 =-4.0*(deti(3)*d5 + detr(3)*d6)
 
!     COMPUTE  GK*GK
 
 d3 = gkr*gkr - gki*gki
 d4 = 2.0*gkr*gki
 CALL csumm  (d3,d4,2*igk,d1,d2,ipdet(3)+id2,d5,d6,id1)
 CALL csqrtn (d5,d6,id1,rootr,rooti,iroot)
 CALL csumm  (gkr,gki,igk,rootr,rooti,iroot,d9,d10,id3)
 CALL csumm  (gkr,gki,igk,-rootr,-rooti,iroot,d7,d8,id4)
 IF (icnt == 1) GO TO 290
 d1  = d9
 d2  = d10
 id1 = id3
 d5  = d9*d9 + d10*d10
 d6  = d7*d7 + d8*d8
 IF (d5 >= d6) GO TO 310
 d1  = d7
 d2  = d8
 id1 = id4
 GO TO 310
 
!     COMPUTE  NUMERATOR  EQ. 11
 
 290 d1 = d9
 d2 = d10
 id1= id3
 m  = 2
 GO TO 310
 300 d1 = d7
 d2 = d8
 id1= id4
 m  = 1
 310 d3 =-2.0*(detr(3)*deltkr - deti(3)*deltki)
 d4 =-2.0*(deti(3)*deltkr + detr(3)*deltki)
 d5 = d1*d1 + d2*d2
 d6 = 10.0**(ipdet(3) - id1)
 lamk1r = d6*(d3*d1 + d4*d2)/d5
 lamk1i = d6*(d4*d1 - d3*d2)/d5
 hkp1r  = lamk1r*hkr - lamk1i*hki
 hkp1i  = lamk1i*hkr + lamk1r*hki
 pr = pkr(3) + hkp1r
 pi = pki(3) + hkp1i
 IF (icnt /= 1) GO TO 370
 dt3 = 0.0D0
 DO  i = 1,3
   dt3 = dt3+DSQRT((pkr(i)-pr)**2 + (pki(i)-pi)**2)
 END DO
 IF (dt3 > dt2) GO TO 340
 ptr = pr
 pti = pi
 dt2 = dt3
 340 IF (m == 2) GO TO 300
 pr  = ptr
 pi  = pti
 
!     DO RANGE CHECKS
 
 
!     COMPUTE U VECTOR
 
 370 xu = pr - alph1
 yu = pi - w1
 lu = SQRT(xu*xu + yu*yu)
 IF (lu == 0.0) GO TO 380
 xu = xu/lu
 yu = yu/lu
 x  = lu*(xu*xvr + yu*yvr)
 y  = lu*(yu*xvr - xu*yvr)
 IF (ABS(y) > xl/2.0 .OR. x < 0.0 .OR. x > rl) GO TO 400
 
!     SEE IF POINT IS IN LOCAL REGION
 
 380 xu = pr - aloc1
 yu = pi - wloc1
 lu = SQRT(xu*xu + yu*yu)
 IF (lu == 0.0) GO TO 390
 xu = xu/lu
 yu = yu/lu
 y  = lu*(yu*xvr-xu*yvr)
 x  = lu*(xu*xvr+yu*yvr)
 IF (ABS(y) > xl/2.0 .OR. x < 0.0 .OR. x > rll) GO TO 190
 
!     TRY FOR CONVERGENCE
 
 390 ASSIGN 450 TO dretn
 GO TO 810
 
!     PREDICTED OUTSIDE BIG REGION
 
 400 noutsd = noutsd + 1
 GO TO 190
 
!     BEGIN CONVERGENCE TESTS
 
 450 IF (icnt <= 2) GO TO 520
 h1bar = DSQRT(hk1r*hk1r + hk1i*hk1i)
 h2bar = DSQRT(hkr*hkr + hki*hki)
 h3bar = DSQRT(hkp1r*hkp1r + hkp1i*hkp1i)
 460 test  = epsi*rl
 IF (h1bar > test*1.0E7) GO TO 480
 IF (h2bar > test*1.0E4) GO TO 480
 IF (h3bar >      h2bar) GO TO 470
 IF (h3bar >       test) GO TO 480
 GO TO 550
 470 IF (h2bar <=  1.0E-7*rl) GO TO 550
 480 icnt = icnt + 1
 IF (icnt-nit < 0) THEN
   GO TO   530
 ELSE IF (icnt-nit == 0) THEN
   GO TO   500
 END IF
 490 ifail = 1
 nfail = nfail + 1
 GO TO 190
 500 IF (nchang < numint .AND. ifail == 1) GO TO 510
 GO TO 490
 510 epsi   = epsi*10.0
 nchang = nchang + 1
 GO TO 460
 
!     CONTINUE INTERATIONS
 
 520 icnt = icnt + 1
 530 DO  i = 1,2
   pkr(i)  = pkr(i+1)
   pki(i)  = pki(i+1)
   ipdet(i)= ipdet(i+1)
   detr(i) = detr(i+1)
   deti(i) = deti(i+1)
 END DO
 pkr(3)  = pr
 pki(3)  = pi
 detr(3) = dr
 deti(3) = di
 ipdet(3)= powr
 GO TO 280
 
!     ACCEPT CURRENT EIGENVALUE
 
 550 FILE   = lama
 nfound = nfound + 1
 ifpass = 1
 IF (nfound > 1) im1 = 3
 CALL OPEN (*950,lama,iz(ibuf),im1)
 zd(1  ) = pr
 zd(iz2) = pi
 CALL WRITE (lama,zd(1),4,1)
 CALL CLOSE (lama,2)
 
!     BUILD LOAD FOR FBS
 
 IF (minda == 0.0D0) minda = 1.0D-8
 SIGN =-SIGN
 d1   = nrow
 d2   = nfound
 j    = 2*nrow
 DO  i = 1,j,2
   k    = (i+1)/2
   zd(i  ) = SIGN*minda/(1.0D0+(1.0D0-FLOAT(k)/d1)*d2)
   zd(i+1) = 0.0D0
 END DO
 iscr2(1) = sr2
 iscr2(7) = fu(7)
 CALL cdtfbs (zd(1),zd(j+1),iz(ibuf),iscr2,nrow)
 
!     NORMALIZE
 
 d1 = 0.0D0
 DO  i = 1,j,2
   d2 = zd(i)*zd(i) + zd(i+1)*zd(i+1)
   IF (d2 < d1) CYCLE
   d3 = zd(i  )
   d4 = zd(i+1)
   d1 = d2
 END DO
 IF (inorm == 0) GO TO 600
 580 DO  i = 1,j,2
   d5 = (zd(i)*d3  + zd(i+1)*d4)/d1
   zd(i+1) = (d3*zd(i+1) - d4*zd(i))/d1
   zd(i  ) = d5
 END DO
 GO TO 610
 600 jj = 2*isil
 d2 = zd(jj)*zd(jj) + zd(jj-1)*zd(jj-1)
 IF (d2 == 0.0D0 .OR. d1/d2 > 1.0D6) GO TO 580
 d3 = zd(jj-1)
 d4 = zd(jj  )
 d1 = d2
 GO TO 580
 
!     WRITE OUT NORMALIZED VECTOR
 
 610 FILE = phid
 CALL OPEN  (*950,phid,iz(ibuf),im1)
 CALL WRITE (phid,zd(1),4*nrow,1)
 CALL CLOSE (phid,2)
 
!     STORE ACCEPTED VALUE
 
 iz(irgp+7) = iz(irgp+7) + 1
 nf = nf + 1
 j  = (irgp+1)/2 + 2*nf + 10*NE + 14
 zd(j  ) = pr
 zd(j+1) = pi
 ifail   = 0
 
!     CHECK FOR STARTING POINT MOVES
 
 j   = ireg1
 i   = 1
 620 dt1 = 200.0*epsi*epsi*z(j+9)
 m   = 2*iz(j+5) + 2
 k   = (j+1)/2 + 5
 l   = 1
 630 k   = k + 2
 kkk = j + 16*iz(j+5) + 26 + 2*l
 IF (DSQRT((zd(k)-pi)**2+(zd(k-1)-pr)**2) >= dt1) GO TO 650
 
!     SHIFT STARTING POINT
 
 d2  = 1000.0*epsi*epsi*z(j+9)
 zd(k-1) = DSIGN((z(j+2)-z(j  ))/z(j+9)*d2+zd(k-1),zd(k-1))
 zd(k  ) = DSIGN((z(j+3)-z(j+1))/z(j+9)*d2+zd(k  ),zd(k  ))
 nmoves  = nmoves + 1
 
!     IF  DETERMINANT EVALUATED - REEVALUATE FOR SHIFT
 
 IF (iz(kkk) == 0) GO TO 670
 dt2 = pr
 dt3 = pi
 pr  = zd(k-1)
 pi  = zd(k  )
 ASSIGN 640 TO dretn
 GO TO 810
 640 pr = dt2
 pi = dt3
 kk = k + 4*iz(j+5) + 4
 zd(kk  ) = di
 zd(kk-1) = dr
 iz(kkk+1)= powr
 GO TO 660
 
!     SWEEP ACCEPTED VALUE FROM STORED  DETM-S
 
 650 kk = k + 4*iz(j+5) + 4
 d2 = zd(k-1) - pr
 d3 = zd(k  ) - pi
 d4 = d2*d2 + d3*d3
 d5 = (zd(kk-1)*d2 + zd(kk)*d3)/d4
 zd(kk  ) = (zd(kk)*d2 - zd(kk-1)*d3)/d4
 zd(kk-1) = d5
 
!     SWEEP CONJUGATES S
 
 IF (icmpx == 1 .OR. DABS(pi) < 1000.0*z(j+9)*epsi) GO TO 660
 d3 = zd(k) + pi
 d4 = d2*d2 + d3*d3
 d5 = (zd(kk-1)*d2 + zd(kk)*d3)/d4
 zd(kk) = (zd(kk)*d2 - zd(kk-1)*d3)/d4
 zd(kk-1) = d5
 660 zdkm1 = zd (kk-1)
 zdk   = zd (kk  )
 izk   = iz (kkk+1)
 CALL cdetm3 (zdkm1,zdk,izk)
 zd(kk-1) = zdkm1
 zd(kk  ) = zdk
 iz(kkk+1)= izk
 670 l = l + 1
 IF (l <= m) GO TO 630
 j = iz(j+8)
 i = i + 1
 IF (i <= nregn) GO TO 620
 CALL klock  (itime2)
 CALL tmtogo (itleft)
 IF (2*(itime2-itime1) > itleft .AND. nfound /= ndesrd) GO TO 700
 IF (nf < nd) GO TO 260
 
!     FIND NEXT REGION LACKING ROOTS
 
 680 IF (iz(irgp+8) == 0) GO TO 140
 irgp = iz(irgp+8)
 IF (iz(irgp+6) > iz(irgp+7)) GO TO 170
 GO TO 680
 690 iterm = 2
 GO TO 710
 
!     INSUFFICIENT TIME
 
 700 IF (nmsgx >= maxgx) nmsgx = maxgx - 1
 CALL mesage (45,ndesrd-nfound,NAME)
 iterm = 3
 GO TO 710
 
!     SINGULAR MATRIX
 
 701 iterm = 4
 GO TO 710
 
!     END OF ROUTINE  PUT OUT SUMMARY
 
 710 CALL gopen (oceigs,iz(ibuf),1)
 CALL WRITE (oceigs,ihead(1),10,0)
 iz(  1) = nfound
 iz(iz2) = npass
 iz(iz3) = nchang
 iz(iz4) = nmoves
 iz(iz5) = ndcomp
 iz(iz6) = nfail
 iz(iz7) = noutsd
 iz(iz8) = iterm
 CALL WRITE (oceigs,iz(1),40,0)
 CALL WRITE (oceigs,head(1),96,1)
 ihead(3) = 3
 ihead(10)= 6
 CALL WRITE (oceigs,ihead,50,0)
 CALL WRITE (oceigs,head,96,1)
 j  = ireg1
 DO  i = 1,nregn
   NE = iz(j+5)
   k  = (j+1)/2+6
   kk = 4*NE + 4
   kd = j + 27 + 16*NE
   NE = 2*NE + 2
   DO  l = 1,NE
     iz(1)  = l
     z(iz2) = zd(k  )
     z(iz3) = zd(k+1)
     m  = k  + kk
     kd = kd + 2
     iz(iz6) = iz(kd)
     
!     CONVERT TO MAGNITUDE AND PHASE  SCALE ON MAGNITIDE
!     PHASE IN DEGRESS BETWEEN 0 AND 360
     
     d1 = DSQRT(zd(m)*zd(m) + zd(m+1)*zd(m+1))
     IF (d1 ==  0.0D0) GO TO 760
     720 IF (d1 > 10.0D0) GO TO 740
     730 IF (d1 <  1.0D0) GO TO 750
     GO TO 770
     740 d1 = d1*0.1D0
     iz(iz6) = iz(iz6) + 1
     GO TO 720
     750 d1 = d1*10.0D0
     iz(iz6) = iz(iz6) - 1
     GO TO 730
     
!     NOT  EVALUATED
     
     760 z(iz4) = 0.0
     z(iz5) = 0.0
     GO TO 780
     770 z(iz4) = d1
     
!     COMPUTE PHASE
     
     z(iz5) = DATAN2(zd(m+1),zd(m))*raddeg
     
!     DETERMINE QUADRANT
     
     IF (z(iz5) < 0.) z(iz5) = z(iz5) + 360.0
     780 CONTINUE
     CALL WRITE (oceigs,iz(1),6,0)
     k = k + 2
   END DO
   j = iz(j+8)
 END DO
 CALL CLOSE (oceigs,1)
 fa(1) = oceigs
 CALL wrttrl (fa(1))
 RETURN
 
!     INTERNAL SUBROUTINE TO EVALUATE DR,DI AT PR,PI
 
 810 ndcomp = ndcomp + 1
 
!     SET UP FOR ADD
 
 bmcb(1) = pr
 bmcb(2) = pi
 cmcb(1) = pr*pr - pi*pi
 cmcb(2) = 2.*pr*pi
 CALL sadd (z(1),z(1))
 fa(1) = -IABS(fa(1))
 IF (nosing == 0) GO TO 821
 isave = sr2
 sr2   = scr8
 scr8  = isave
 821 CALL tmtogo (kk)
 IF (kk <= 0) GO TO 700
 ib = 0
 CALL cdcomp (*930,z(1),z(1),z(1))
 nosing = 1
 IF (iprt /= 0) WRITE (otpe,831) pr,pi,dr,di,powr
 831 FORMAT (10X,4D16.7,i8)
 
!     SCALE DETERMINANT BY POLES AND EIGENVALUES PREVIOUSLY FOUND
 
 id1 = ireg1
 DO  id = 1,nregn
   id2 = iz(id1+5)
   kk  = iz(id1+7)
   IF (kk == 0) GO TO 870
   kd = 14 + 10*id2 + (id1+1)/2
   DO  ll = 1,kk
     kd = kd + 2
     d1 = pr - zd(kd  )
     d2 = pi - zd(kd+1)
     d3 = d1*d1  + d2*d2
     d4 = (dr*d1 + di*d2)/d3
     d5 = (di*d1 - dr*d2)/d3
     dr = d4
     di = d5
     IF (icmpx == 1) GO TO 850
     
!     SWEEP COMPLEX CONJUGATE ROOTS
     
     IF (DABS(zd(kd+1)) < 1000.0*z(id1+9)*epsi) GO TO 850
     d2 = pi + zd(kd+1)
     d3 = d1*d1  + d2*d2
     d4 = (dr*d1 + di*d2)/d3
     d5 = (di*d1 - dr*d2)/d3
     dr = d4
     di = d5
     850 CALL cdetm3 (dr,di,powr)
   END DO
   870 id1 = iz(id1+8)
 END DO
 
!     SWEEP POLES
 
 IF (npole == 0) GO TO 940
 id1 = ipoles
 DO  id = 1,npole
   d1 = pr - z(id1+1)
   d2 = pi - z(id1+2)
   d3 = 1.0D0
   d4 = 0.0D0
   kd = iz(id1+3)
   DO  id2 = 1,kd
     d5 = d1*d3 - d2*d4
     d6 = d2*d3 + d1*d4
     d3 = d5
     d4 = d6
   END DO
   d1 = d3*d3  + d4*d4
   d2 = (dr*d3 + di*d4)/d1
   d5 = (di*d3 - dr*d4)/d1
   dr = d2
   di = d5
   id1= id1 + 4
   
!     SCALE AGAIN
   
   CALL cdetm3 (dr,di,powr)
 END DO
 GO TO 940
 
!     SINGLULAR MATRIX
 
 930 dr    = 0.0D0
 di    = 0.0D0
 powr  = 0
 ising = ising + 1
 minda = 1.0E-11
 IF (nosing == 0) GO TO 940
 nosing= 0
 isave = sr2
 sr2   = scr8
 scr8  = isave
 
!     RETURN
 
 940 IF (iprt /= 0) WRITE (otpe,831) pr,pi,dr,di,powr
 GO TO dretn, (240,450,640)
 
!     ERROR  MESAGES
 
 950 ip1 = -1
 960 CALL mesage (ip1,FILE,NAME)
 970 ip1 = -2
 GO TO 960
 980 ip1 = -3
 GO TO 960
 990 ip1 = -7
 GO TO 960
 1000 ip1 = -8
 GO TO 960
END SUBROUTINE cdetm
