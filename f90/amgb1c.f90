SUBROUTINE amgb1c (q)
     
    !     UNSTEADY FLOW ANALYSIS OF A SUPERSONIC CASCADE
 
 
    COMPLEX, INTENT(OUT)                     :: q(nstns,nstns)
    INTEGER :: sln
    COMPLEX :: sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
        am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
        am5t,ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,  &
        f1s,c1,c2p,c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,  &
        sum2,f2,blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,  &
        c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,  &
        bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,  &
        pres1,pres2,pres3,pres4,qres4,fqa,fqb,fq7,presu, presl, gusamp
    DIMENSION       gye(29,29),gee(29,40),presu(29),presl(29),xup(29),  &
        xtemp(29),geetmp(29,20),xlow(29),aye(10,29),  &
        INDEX(29,3), pres1(21),pres2(21),  &
        pres3(21),pres4(21),qres4(21),sbkde1(201),  &
        sbkde2(201),sumsv1(201),sumsv2(201),svkl1(201),  &
        svkl2(201),xlsv1(21),xlsv2(21),xlsv3(21),xlsv4(21)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /system/ sysbuf,ibbout
    COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigm ,rfreq
    COMMON /bamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
        refmac,refden,refvel,refflo,sln,nstnsx,stg,  &
        chord,radius,bspace,mach,den,vel,flowa,amachd, redfd,blspc,amachr,tsonic
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
 
    IF (nstns > 10) GO TO 420
 
    redf  = redfd
    amach = amachd
    ai    = CMPLX(0.0,1.0)
    pi    = 3.1415927
    pitcor= blspc
    stag  = 90.0 - stg
    sigma = -sigm*pi/180.0
    beta  = SQRT(amach**2 - 1.0)
    scrk  = redf*amach/(beta**2)
    del   = scrk*amach
    amu   = redf/(beta**2)
    sp    = pitcor*COS(stag*pi/180.0)*2.0
    sn    = pitcor*SIN(stag*pi/180.0)*2.0
    sps   = sp
    sns   = sn*beta
    dstr  = SQRT(sps**2 - sns**2)
    sps1  = ABS(sps - sns)
    IF (sps1 < .00001) GO TO 400
 
    !     ZERO OUT GEE
 
    nstns2 = 2*nstns
    nstns4 = 4*nstns
    DO  i = 1,29
        DO  j = 1,nstns4
            gee(i,j) = 0.0
        END DO
    END DO
    pitaxs = 0.0
    amoaxs = 0.
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
        xl  = xl  + temp
    END DO
    xl  = sns + 2.0 - sps
    temp= (sps-sns)/rl1
    DO  jl = 1,10
        xlsv4(jl) = xl
        xl = xl + temp
    END DO
 
    !     ACCUMULATE PRESSURE VECTORS INTO G-MATRIX
 
    DO  nm = 1,nstns
        ntimes = 1
        IF (nm > 2) ntimes = 2
        DO  nmm = 1,ntimes
     
            !     DEFINE -----------------------------
            !            ALPAMP - PITCHING AMP
            !            DISAMP - PLUNGING AMP
            !            GUSAMP - GUST AMP
            !            GL -GUST WAVE NUMBER
     
            alpamp = 0.0
            IF (nm == 2) alpamp = 1.0
            disamp = 0.0
            IF (nm == 1) disamp = 1.0
            gusamp = 0.0
            gl     = 0.0
            IF (nm > 2 .AND. nmm == 1) gusamp =-redf/2.0 + (nm-2)*pi/4.0
            IF (nm > 2 .AND. nmm == 1) gl = (nm-2)*pi/2.0
            IF (nm > 2 .AND. nmm == 2) gusamp = redf/2.0 + (nm-2)*pi/4.0
            IF (nm > 2 .AND. nmm == 2) gl =-(nm-2)*pi/2.0
     
            a = (1.0+ai*redf*pitaxs)*alpamp - ai*redf*disamp
            b =-ai*redf*alpamp
            IF (gl == 0.0) GO TO 50
            a = gusamp
            b = 0.0
50      CONTINUE
        CALL suba
     
        !     FIND  DELTA P(LOWER-UPPER)
     
        DO  nx = 1,10
            presu(nx) = pres1(nx)
            xup(nx)   = xlsv1(nx)
            IF (nx == 10) GO TO 60
            nxx = nx + 20
            presl(nxx) = pres4(nx+1)
            xlow( nxx) = xlsv4(nx+1)
            GO TO 70
60          presu(nx) = (pres1(10) + pres2(1))/2.0
            xup(10)   = (xlsv1(10) + xlsv2(1))/2.0
70      CONTINUE
        END DO
        DO  nx = 1,20
            nxx = nx + 10
            IF (nx == 20) GO TO 90
            presu(nxx) = pres2(nx+1)
            xup  (nxx) = xlsv2(nx+1)
            presl(nx)  = pres3(nx)
            xlow( nx)  = xlsv3(nx)
            GO TO 100
90          presl(20) = (pres3(20) + pres4(1))/2.0
            xlow(20)  = (xlsv3(20) + xlsv4(1))/2.0
100     CONTINUE
        END DO
        nm2 = nm + nstns
        nm3 = nm + 2*nstns
        nm4 = nm + 3*nstns
        DO  nmmm = 1,29
            gee(nmmm,nm)  = gee(nmmm,nm ) + REAL(presl(nmmm))
            gee(nmmm,nm2) = gee(nmmm,nm2) + AIMAG(presl(nmmm))
            gee(nmmm,nm3) = gee(nmmm,nm3) + REAL(presu(nmmm))
            gee(nmmm,nm4) = gee(nmmm,nm4) + AIMAG(presu(nmmm))
        END DO
    END DO
END DO
 
!     NOW DEFINE  I-MATRIX (NSTNS X 29)
 
aye(1,1) = 2.0
con = 1.0
aye(1,2) = 2.0
n1n = 27
DO  j = 1,n1n
    aye(1,j+2) = con*4.0/j/pi
    con = 1.0 - con
END DO
aye(2,1) = 2.0
aye(2,2) = 2.66666667
con = 1.0
DO  j = 1,n1n
    aye(2,j+2) = con*4/j/pi
    con = -con
END DO
DO  i = 3,nstns
    DO  j = 2,28
        con = 0.0
        IF ((i-1) == j) con = 1.0
        aye(i,j+1) = con
    END DO
END DO
DO  j = 3,nstns
    aye(j,1) = aye(1,j)
    aye(j,2) = aye(2,j)
END DO
 
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
 
!     PUT PRESL PART OF GEE IN GEETMP
 
DO  i = 1,29
    DO  j = 1,nstns2
        geetmp(i,j) = gee(i,j)
    END DO
END DO
 
!     SOLVE FOR G-INVERSE G IN GEE MATRIV
!     ISING = 1 NON-SINGULAR (GYE)
!     ISING = 2  SIGULAR     (GYE)
!     INDEX IS WORK STORAGE FOR ROUTINE INVERS
 
ising = -1
CALL invers (29,gye,29,geetmp,nstns2,determ,ising,INDEX)
IF (ising == 2) GO TO 410
 
!     NOW  MULTIPLY  I*G-INVERSE*G(DELTA P'S)
 
DO  j = 1,nstns
    DO  k = 1,nstns
        nf   = k + nstns
        sumi = 0.0
        sumr = 0.0
        DO  i = 1,29
            sumr = aye(j,i)*geetmp(i,k ) + sumr
            sumi = aye(j,i)*geetmp(i,nf) + sumi
        END DO
     
        !  NOTE - NOTE THAT DUE TO CEXP( - I*OMEGA*T) TYPE OF TIME DEPENDENCE
        !         IN UCAS DEVELOPMENT, Q IS DEFINED AS THE COMPLEX CONJUGATE
        !         OF 'USUAL' Q
     
        q(j,k) = 2.0*CMPLX(sumr,-sumi)
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
DO  i = 2, 29
    gye(i,2) = xtemp(i)
END DO
 
!     PUT PRESU PART OF GEE IN GEETMP
 
DO  i = 1,29
    DO  j = 1,nstns2
     
        nsns2 = nstns2 + j
        geetmp(i,j) = gee(i,nsns2)
    END DO
END DO
 
!     SOLVE FOR G-INVERSE G IN GEETMP MATRIX
!     ISING = 1  NON-SINGULAR (GYE)
!     ISING = 2  SINGULAR GYE
!     INDEX IS WORK STORAGE FOR ROUTINE INVERS
 
ising = -1
CALL invers (29,gye,29,geetmp,nstns2,determ,ising,INDEX)
 
IF (ising == 2) GO TO 410
 
!     MULTIPLY I*G-INVERS*G
 
DO  j = 1,nstns
    DO  k = 1,nstns
        nf = k + nstns
        sumi = 0.0
        sumr = 0.0
        DO  i = 1,29
       
            sumr = aye(j,i)*geetmp(i,k ) + sumr
            sumi = aye(j,i)*geetmp(i,nf) + sumi
       
        END DO
     
        q(j,k) = q(j,k) - 2.0*CMPLX(sumr,-sumi)
    END DO
END DO
 
RETURN
 
400 WRITE (ibbout,500) ufm
GO TO 430
410 WRITE (ibbout,510) ufm
GO TO 430
420 WRITE (ibbout,520) ufm,sln,nstns
430 CALL mesage (-61,0,0)
RETURN
 
500 FORMAT (a23,' - AMG MODULE -SUBROUTINE AMGB1C', /39X,  &
    'AXIAL MACH NUMB. IS EQUAL TO OR GREATER THAN ONE.')
510 FORMAT (a23,' - AMG MODULE - LARGE G-MATRIX IS SINGULAR IN ',  &
    'ROUTINE AMGBIC.')
520 FORMAT (a23,' - AMG MODULE - NUMBER OF COMPUTING STATIONS ON ',  &
    'STREAMLINE',i8,4H is ,i3,1H. , /39X,'SUPERSONIC CASCADE',  &
    ' ROUTINE AMGB1C ALLOWS ONLY A MAXIMUM OF 10.')

END SUBROUTINE amgb1c
