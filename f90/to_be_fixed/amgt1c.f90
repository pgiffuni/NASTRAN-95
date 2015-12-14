SUBROUTINE amgt1c (q,nstns2,c1sbar,c2sbar)
     
!     SUPERSONIC CASCADE CODE FOR SWEPT TURBOPROPS.
 
 
 COMPLEX, INTENT(OUT)                     :: q(nstns2,nstns2)
 INTEGER, INTENT(IN)                      :: nstns2
 REAL, INTENT(IN)                         :: c1sbar
 REAL, INTENT(IN)                         :: c2sbar
 INTEGER :: sln
 REAL :: m2sbar
 COMPLEX :: sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,  &
     f1s,c1,c2p,c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,  &
     sum2,f2,blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,  &
     c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,  &
     bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,  &
     pres1,pres2,pres3,pres4,qres4,fqa,fqb,fq7,presu, presl, gusamp
 DIMENSION       gye(29,29),gee(29,80),presu(29),presl(29),xup(29),  &
     xtemp(29),geetmp(29,40),xlow(29),aye(10,29),  &
     INDEX(29,3), pres1(21),pres2(21),  &
     pres3(21),pres4(21),qres4(21),sbkde1(201),  &
     sbkde2(201),sumsv1(201),sumsv2(201),svkl1(201),  &
     svkl2(201),xlsv1(21),xlsv2(21),xlsv3(21),xlsv4(21)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /tamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refswp,sln,nstnsx,stg,  &
     chord,dcbdzb,bspace,mach,den,vel,sweep,amachd, redfd,blspc,amachr,tsonic
 COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigm,rfreq
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 COMMON /blk3  / sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,a,b,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1,c2p,  &
     c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,  &
     blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,c1p,c1n,  &
     bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,bc2,bc3,  &
     bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,pres1,pres2,  &
     pres3,pres4,qres4,fqa,fqb,fq7
 COMMON /blk4  / i,r,y,a1,b1,c4,c5,gl,i6,i7,jl,nl,ri,rt,r5,sn,sp,  &
     xl,y1,amu,gam,idx,inx,nl2,rl1,rl2,rq1,rq2,xl1,  &
     alp1,alp2,gamn,gamp,iner,iout,redf,stag,step,  &
     amach,betnn,betnp,bkap1,xlsv1,xlsv2,xlsv3,xlsv4,  &
     alpamp,amoaxs,gusamp,disamp,pitaxs,pitcor
 
 
!     THEORY DEPENDENT RESTRICTION OF NO MORE THAN 10 COMPUTING
!     STATIONS PER STREAMLINE IS REFLECTED IN CODING.
 
 IF (nstns > 10) GO TO 9993
 
 redf  = redfd
 amach = amachd
 ai    = CMPLX(0.0,1.0)
 pi    = 3.1415927
 pitcor= blspc
 stag  = 90.0 - stg
 sigma = -sigm*pi/180.0
 beta  = SQRT(amach**2-1.0)
 scrk  = redf*amach/(beta**2)
 del   = scrk*amach
 amu   = redf/(beta**2)
 sp    = pitcor*COS(stag*pi/180.0)*2.0
 sn    = pitcor*SIN(stag*pi/180.0)*2.0
 sps   = sp
 sns   = sn*beta
 dstr  = SQRT(sps**2-sns**2)
 sps1  = ABS(sps-sns)
 IF (sps1 < .00001)  GO TO 9991
 
!     PARAMETERS RELATED TO SWEEP CHANGES
 
 csbar  = .25*(den*vel**2*chord**2)/(refden*refvel**2)
 csbar1 = 2.0/chord
 m2sbar = -dcbdzb/chord
 c2ssch = csbar1*c2sbar
 csblsb = csbar*csbar1
 csbm2s = csbar*m2sbar
 tanlam = TAN(sweep*pi/180.)
 dlsdzb = dcbdzb/2.0
 td     = tanlam*dlsdzb
 
!     ZERO OUT GEE
 
 nstns4 = 4*nstns
 nstns8 = 8 * nstns
 DO  i = 1,29
   DO  j = 1,nstns8
     gee(i,j) = 0.0
   END DO
 END DO
 pitaxs   = 0.0
 amoaxs   = 0.
 CALL asycon
 CALL akp2
 rl1 = 9
 s1  = sps - sns
 aa  = s1/rl1
 xlsv1(1) = 0.0
 DO  jl = 1,9
   xlsv1(jl+1) = jl*aa
 END DO
 aa  = sps - sns
 rl2 = 19
 s1  = 2.0 + sns - sps
 temp= s1/rl2
 xl  = aa
 DO  jl = 1,20
   xlsv2(jl) = xl
   xlsv3(jl) = xl + sns - sps
   xl = xl + temp
 END DO
 xl = sns + 2.0 - sps
 temp = (sps-sns)/rl1
 DO  jl = 1,10
   xlsv4(jl) = xl
   xl = xl + temp
 END DO
 
!     ACCUMULATE PRESSURE VECTORS INTO G-MATRIX
 
 DO  nm = 1,nstns
   ntimes = 1
   IF (nm > 2) ntimes = 2
   DO  nmm = 1,ntimes
     
     jndx = 0
     5000 IF (jndx == 0) GO TO 5010
     
     IF (nm > 2) GO TO 5020
     
     gl = 0.0
     IF (nm == 1) a = tanlam/csbar1
     IF (nm == 1) b = 0.0
     IF (nm == 2) a = 0.0
     IF (nm == 2) b = tanlam/csbar1
     
     GO TO 2047
     
     5020 IF (nmm == 1) gusamp =-ai*tanlam/csbar1/2.0
     IF (nmm == 1) gl = (nm-2)*pi/2.0
     IF (nmm == 2) gusamp = ai*tanlam/csbar1/2.0
     IF (nmm == 2) gl =-(nm-2)*pi/2.0
     
     a = gusamp
     b = 0.0
     
     GO TO 2047
     
!     DEFINE -----------------------------
!              ALPAMP - PITCHING AMP
!              DISAMP - PLUNGING AMP
!              GUSAMP - GUST AMP
!              GL -GUST WAVE NUMBER
     5010 alpamp = 0.0
     IF (nm == 2) alpamp = 1.0
     disamp = 0.0
     IF (nm == 1) disamp = 1.0
     gusamp = 0.0
     gl     = 0.0
     IF (nm > 2 .AND. nmm == 1) gusamp =-(redf+ai*td)/2.+(nm-2)*pi/4.
     IF (nm > 2 .AND. nmm == 1) gl = (nm-2)*pi/2.0
     IF (nm > 2 .AND. nmm == 2) gusamp = (redf+ai*td)/2.+(nm-2)*pi/4.
     IF (nm > 2 .AND. nmm == 2) gl =-(nm-2)*pi/2.0
     
     a = (1.0+ai*redf*pitaxs)*alpamp - (ai*redf-td)*disamp
     b =-(ai*redf-td)*alpamp
     IF (gl == 0.0) GO TO 2047
     a = gusamp
     b = 0.0
     2047 CONTINUE
     
     CALL suba
     
!     FIND  DELTA P(LOWER-UPPER)
     
     DO  nx = 1,10
       presu(nx) = pres1(nx)
       xup(nx)   = xlsv1(nx)
       IF (nx == 10) GO TO 55
       nxx = nx + 20
       presl(nxx) = pres4(nx+1)
       xlow( nxx) = xlsv4(nx+1)
       GO TO 610
       55 presu(nx) = (pres1(10) + pres2(1))/2.0
       xup(10)   = (xlsv1(10) + xlsv2(1))/2.0
       610  CONTINUE
     END DO
     DO  nx = 1,20
       nxx = nx + 10
       IF (nx == 20) GO TO 65
       presu(nxx) = pres2(nx+1)
       xup  (nxx) = xlsv2(nx+1)
       presl(nx)  = pres3(nx  )
       xlow( nx)  = xlsv3(nx  )
       GO TO 710
       65 presl(20) = (pres3(20) + pres4(1))/2.0
       xlow(20)  = (xlsv3(20) + xlsv4(1))/2.0
       710 CONTINUE
     END DO
     
     jx   = jndx*4*nstns
     nmz  = nm + jx
     nm2z = nm + nstns   + jx
     nm3z = nm + 2*nstns + jx
     nm4z = nm + 3*nstns + jx
     
     DO  nmmm = 1,29
       gee(nmmm,nmz ) = gee(nmmm,nmz ) + REAL(presl(nmmm))
       gee(nmmm,nm2z) = gee(nmmm,nm2z) + AIMAG(presl(nmmm))
       gee(nmmm,nm3z) = gee(nmmm,nm3z) + REAL(presu(nmmm))
       gee(nmmm,nm4z) = gee(nmmm,nm4z) + AIMAG(presu(nmmm))
       
     END DO
     
     IF (jndx /= 0) CYCLE
     jndx = 1
     GO TO 5000
     
   END DO
 END DO
 
!     NOW DEFINE  I-MATRIX (NSTNS X 29)
 
 aye(1,1) = c1sbar*2.0 + c2ssch*2.0
 aye(1,2) = c1sbar*8.0/3.0 + c2ssch*2.0
 aye(2,1) = c1sbar*8.0/3.0 + c2ssch*2.0
 aye(2,2) = c1sbar*4.0 + c2ssch*8.0/3.0
 
 conz1 = 1.0
 
 DO  i = 3,nstns
   conz4 = (1.+conz1 )*2./(pi*(j-2))
   conz5 = conz1*4./ (pi*(j-2))
   conz6 = conz1*8./(pi*(j-2)) - (1.+conz1)*16./(pi*(j-2))**3
   
   aye(i,1) = c1sbar*conz5 + c2ssch*conz4
   aye(i,2) = c1sbar*conz6 + c2ssch*conz5
   conz1    = -conz1
 END DO
 
 conz1 = 1.0
 
 DO  j = 3,29
   conz4 = (1.+conz1)*2./(pi*(j-2))
   conz5 = conz1*4./(pi*(j-2))
   conz6 = conz1*8./(pi*(j-2)) - (1.+conz1)*16./(pi*(j-2))**3
   
   aye(1,j) = c1sbar*conz5 + c2ssch*conz4
   aye(2,j) = c1sbar*conz6 + c2ssch*conz5
   conz1    = -conz1
 END DO
 
 DO  i = 3, nstns
   
   DO  j = 3,29
     conz1 = 0.0
     IF (j == i) GO TO 286
     IF ((i+j)/2*2 == (i+j)) GO TO 285
     conz1 = -16.*(i-2)*(j-2)/(pi*pi*(i-j)*(i-j)*(i+j-4)**2)
     285 conz2 = 0.0
     GO TO 284
     286 conz1 = 1.0
     conz2 = 1.0
     284 aye(i,j) = c1sbar*conz1 + c2ssch*conz2
   END DO
 END DO
 284   CONTINUE
 
 
!     Q DUE TO PRESL ONLY
 
!     NOW DEFINE LARGE G MATRIX
 
 DO  i = 1,29
   gye(1,i) = 0.0
   gye(i,1) = 1.0
 END DO
 
!     PUT XLOW IN XTEMP
 
 DO  i = 1,29
   xtemp(i) = xlow(i)
 END DO
 DO  j = 3,29
   const = (j-2)*pi/2.0
   DO  i = 2,29
     gye(i,j) = SIN(const*xtemp(i))
   END DO
 END DO
 DO  i = 2,29
   gye(i,2) = xtemp(i)
 END DO
 
!     PUT PRESL PARTS OF GEE IN GEETMP (UNPRIMED AND PRIMED TERMS)
 
 DO  i = 1,29
   DO  j = 1,nstns2
     geetmp(i,j) = gee(i,j)
     geetmp(i,j+nstns2) = gee(i,j+nstns4)
   END DO
 END DO
 
!     SOLVE FOR G-INVERSE G IN GEE MATRIV
!     ISING = 1  NON-SINGULAR (GYE)
!     ISING = 2  SIGULAR      (GYE)
!     INDEX IS WORK STORAGE FOR ROUTINE INVERS
 
 ising = -1
 CALL invers (29,gye,29,geetmp,nstns4,determ,ising,INDEX)
 IF (ising == 2) GO TO 9992
 
!     NOW  MULTIPLY  I*G-INVERSE*G(DELTA P'S)
 
 DO    j = 1,nstns
   DO    k = 1,nstns
     
     sumr1 = 0.0
     sumi1 = 0.0
     sumr2 = 0.0
     sumi2 = 0.0
     
     DO  i = 1,29
       sumr1 = sumr1 + aye(j,i)*geetmp(i,k)
       sumi1 = sumi1 + aye(j,i)*geetmp(i,k+nstns)
       sumr2 = sumr2 + aye(j,i)*geetmp(i,k+nstns4)
       sumi2 = sumi2 + aye(j,i)*geetmp(i,k+nstns+nstns4)
     END DO
     
     conz1 = csblsb*sumr1 + csbm2s*sumr2
     conz2 = csblsb*sumi1 + csbm2s*sumi2
     conz3 = csbar*sumr2
     conz4 = csbar*sumi2
     
     q(j,k      ) = 2.0*CMPLX(conz1,-conz2)
     q(j,k+nstns) = 2.0*CMPLX(conz3,-conz4)
     q(j+nstns,k) = (0.0,0.0)
     q(j+nstns,k+nstns) = (0.0,0.0)
   END DO
 END DO
 
!     FINALLY, Q DUE TO (PRESL-PRESU) IS COMPUTED BY SUBTRACTING Q DUE
!     TO PRESU FROM Q DUE TO PRESL ABOVE
 
!     LARGE G MATRIX
 
 DO  i = 1,29
   gye(1,i) = 0.0
   gye(i,1) = 1.0
 END DO
 
!     PUT XUP IN XTEMP
 
 DO  i = 1,29
   xtemp(i) = xup(i)
 END DO
 DO  j = 3,29
   const = (j-2)*pi/2.0
   DO  i = 2,29
     gye(i,j) = SIN(const*xtemp(i))
   END DO
 END DO
 DO  i = 2,29
   gye(i,2) = xtemp(i)
 END DO
 
!     PUT PRESU PARTS OF GEE IN GEETMP (UNPRIMED AND PRIMED TERMS)
 
 DO  i = 1,29
   DO  j = 1,nstns2
     
     nsns2 = nstns2 + j
     geetmp(i,j) = gee(i,nsns2)
     geetmp(i,nsns2) = gee(i,nsns2+nstns4)
   END DO
 END DO
 
!     SOLVE FOR G-INVERSE G IN GEETMP MATRIX
!     ISING = 1  NON-SINGULAR (GYE)
!     ISING = 2  SINGULAR GYE
!     INDEX IS WORK STORAGE FOR ROUTINE INVERS
 
 ising = -1
 CALL invers (29,gye,29,geetmp,nstns4,determ,ising,INDEX)
 
 IF (ising == 2) GO TO 9992
 
!    MULTIPLY I*G-INVERS*G
 
 DO  j = 1,nstns
   DO  k = 1,nstns
     
     sumr1 = 0.0
     sumi1 = 0.0
     sumr2 = 0.0
     sumi2 = 0.0
     
     DO  i = 1, 29
       sumr1 = sumr1 + aye(j,i)*geetmp(i,k)
       sumi1 = sumi1 + aye(j,i)*geetmp(i,k+nstns)
       sumr2 = sumr2 + aye(j,i)*geetmp(i,k+nstns4)
       sumi2 = sumi2 + aye(j,i)*geetmp(i,k+nstns+nstns4)
     END DO
     
     conz1 = csblsb*sumr1 + csbm2s*sumr2
     conz2 = csblsb*sumi1 + csbm2s*sumi2
     conz3 = csbar*sumr2
     conz4 = csbar*sumi2
     
     q(j,k      ) = q(j,k) - 2.0*CMPLX(conz1,-conz2)
     q(j,k+nstns) = q(j,k+nstns) - 2.0*CMPLX(conz3,-conz4)
   END DO
 END DO
 RETURN
 
 9991 WRITE (ibbout,3000) ufm
 GO TO 9999
 9992 WRITE (ibbout,3001) ufm
 GO TO 9999
 9993 WRITE (ibbout,3002) ufm,sln,nstns
 9999 CALL mesage (-61,0,0)
 RETURN
 
 3000 FORMAT (a23,' - AMG MODULE -SUBROUTINE AMGT1C', /39X,  &
     'AXIAL MACH NUMB. IS EQUAL TO OR GREATER THAN ONE.')
 3001 FORMAT (a23,' - AMG MODULE - LARGE G-MATRIX IS SINGULAR IN ',  &
     'ROUTINE AMGT1C.')
 3002 FORMAT (a23,' - AMG MODULE - NUMBER OF COMPUTING STATIONS ON ',  &
     'STREAMLINE',i8,4H is ,i3,1H. ,/39X,'SUPERSONIC CASCADE ',  &
     'ROUTINE AMGT1C ALLOWS ONLY A MAXIMUM OF 10.')
END SUBROUTINE amgt1c
