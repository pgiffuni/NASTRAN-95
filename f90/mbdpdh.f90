SUBROUTINE mbdpdh (ajjl,f,df,f1,df1,f2,df2,xwte,ywte,parea,capphi,  &
        dphite,dss,q,q1,q2,ndn,nd1,nw1,nwn,kte,kte1,  &
        kte2,nte,nncb,nnsbd,in17,ibuf,a)
     
 
 REAL, INTENT(IN OUT)                     :: ajjl
 REAL, INTENT(OUT)                        :: f(1)
 REAL, INTENT(OUT)                        :: df(1)
 REAL, INTENT(OUT)                        :: f1(1)
 REAL, INTENT(OUT)                        :: df1(1)
 REAL, INTENT(OUT)                        :: f2(1)
 REAL, INTENT(OUT)                        :: df2(1)
 REAL, INTENT(IN)                         :: xwte(1)
 REAL, INTENT(IN OUT)                     :: ywte(1)
 REAL, INTENT(IN)                         :: parea(50,50,3)
 COMPLEX, INTENT(IN)                      :: capphi(1)
 COMPLEX, INTENT(OUT)                     :: dphite(3,nnsbd)
 COMPLEX, INTENT(OUT)                     :: dss(nncb,nnsbd)
 COMPLEX, INTENT(OUT)                     :: q(1)
 COMPLEX, INTENT(OUT)                     :: q1(1)
 COMPLEX, INTENT(OUT)                     :: q2(1)
 INTEGER, INTENT(IN)                      :: ndn(1)
 INTEGER, INTENT(IN)                      :: nd1(1)
 INTEGER, INTENT(IN OUT)                  :: nw1(1)
 INTEGER, INTENT(IN)                      :: nwn(1)
 INTEGER, INTENT(IN OUT)                  :: kte(1)
 INTEGER, INTENT(IN OUT)                  :: kte1(1)
 INTEGER, INTENT(IN OUT)                  :: kte2(1)
 INTEGER, INTENT(OUT)                     :: nte(1)
 INTEGER, INTENT(IN)                      :: nncb
 INTEGER, INTENT(IN)                      :: nnsbd
 INTEGER, INTENT(IN OUT)                  :: in17
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 COMPLEX, INTENT(IN OUT)                  :: a(1)
 LOGICAL :: cntrl2,cntrl1,crank1,crank2,asym,surf,lphi,tebox
 COMPLEX :: dphi,tdh, ws,wf1,wf2,temphi,wphi,sumphi,traile  
 
 COMMON /mboxa/ x(12),y(12),tang(10),ang(10),cotang(10)
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2,nbox,npts0,npts1,  &
     npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,boxl,boxw,  &
     boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 DATA    nhcont,nhdss /4HCONT,4HDSS /
 
 nskp   = 0
 n1     = npts0 + npts1
 CALL gopen (in17,ibuf,0)
 DO  mood = 1,njj
   DO  i = 1,nsbd
     nte(i) = 0
     dphite(1,i) = (0.0, 0.0)
     dphite(2,i) = (0.0, 0.0)
     dphite(3,i) = (0.0, 0.0)
   END DO
   DO  i = 1,ncb
     DO  j = 1,nsbd
       dss(i,j) = (0.0, 0.0)
     END DO
   END DO
   DO  j = 1,kct
     q (j)  = (0.0, 0.0)
     f(j)   = 0.0
     df(j)  = 0.0
   END DO
   IF (.NOT.cntrl1) GO TO 116
   DO  j = 1,kc1t
     q1(j)  = (0.0, 0.0)
     f1(j)  = 0.0
     df1(j) = 0.0
   END DO
   116 IF (.NOT.cntrl2) GO TO 118
   DO  j = 1,kc2t
     q2(j)  = (0.0, 0.0)
     f2(j)  = 0.0
     df2(j) = 0.0
   END DO
   118 CALL fread (in17,z,-nskp,0)
   jj = mood
   IF (jj > npts0) GO TO 140
   CALL fread (in17,f,kct,0)
   CALL fread (in17,df,kct,0)
   nskp = nskp + 2*kct
   GO TO 160
   140 IF (jj > n1) GO TO 150
   CALL fread (in17,f1 ,kc1t,0)
   CALL fread (in17,df1,kc1t,0)
   nskp = nskp + kc1t*2
   GO TO 160
   150 CALL fread (in17,f2 ,kc2t,0)
   CALL fread (in17,df2,kc2t,0)
   nskp = nskp + kc2t*2
   160 CALL bckrec (in17)
   
!     START LOOP FOR ROWS ON PLANFORM
   
   kc  = 0
   kc1 = 0
   kc2 = 0
   DO  i = 1,ncb
     ixr = i - 1
     xb  = boxl*(FLOAT(ixr) + 0.5)
     xbb = xb + boxl/2.0
     
!     BOXES ON PLANE OF MAIN
     
     DO  j = 1,nsbd
       IF (.NOT.(i >= nd1(j) .AND. i <= ndn(j))) CYCLE
       dphi  = (0.0, 0.0)
       wphi  = dphi
       tdh   = (0.0 ,0.0)
       lphi  = .false.
       surf  = .false.
       tebox = .false.
       IF (i >= (nw1(j)+nwn(j))/2) tebox = .true.
       iyr   = j - 1
       yb    = boxw*FLOAT(iyr)
       k     = 1
       IF (yb > y(2)) k = 2
       paw   = parea(i,j,1)
       paf1  = parea(i,j,2)
       paf2  = parea(i,j,3)
       pawf  = paw + paf1 + paf2
       IF (.NOT.tebox .AND. beta > tang(k)) pawf = 1.0
       pad   = 1.0 - pawf
       ws    = (0.0, 0.0)
       wf1   = (0.0, 0.0)
       wf2   = (0.0, 0.0)
       IF (j == 1 .AND. asym) GO TO  800
       IF (j > nsb) GO TO 500
       IF (pad >= 0.995) GO TO  400
       IF (paw < 0.005) GO TO  200
       
       kc = kc + 1
       ws = 2.0*paw*CMPLX(df(kc), ek*f(kc))
       
       200 IF (paf1 < 0.005) GO TO  250
       
       kc1 = kc1 + 1
       wf1 = 2.0*paf1*CMPLX(df1(kc1), ek*f1(kc1))
       
       250 IF (paf2 < 0.005) GO TO 300
       
       kc2 = kc2 + 1
       wf2 = 2.0*paf2*CMPLX(df2(kc2), ek*f2(kc2))
       
       300 tdh    = (ws+wf1+wf2)/(pawf*cr)
       lphi   = .true.
       temphi = sumphi(ixr,iyr,nd1,ndn,capphi,dss,nncb,nnsbd,asym)
       dphi   = tdh*capphi(1) + temphi
       IF (pawf >= .005) surf = .true.
       IF (.NOT.surf .OR. .NOT.tebox) GO TO 350
       nte(j) =  i
       dphite(3,j) = dphite(2,j)
       dphite(2,j) = dphite(1,j)
       dphite(1,j) = dphi
       
       350 IF (pawf > 0.995) GO TO  800
       400 IF (.NOT.tebox) GO TO 500
       xt = xwte(j)
       IF (xt >= xbb) GO TO 420
       dphite(1,j) = traile(xt,j,nte,dphite,nnsbd,boxl)
       IF (xt <= xb) GO TO 450
       420 IF (xt >= xbb+boxl) GO TO 500
       wphi = dphi
       GO TO 800
       450 ex   = ek*(xb-xt)/boxl
       wphi = dphite(1,j)*CMPLX(COS(ex), -SIN(ex))
       GO TO 800
       500 dphi = pawf*dphi
       wphi = (0.0, 0.0)
       IF (.NOT.lphi) temphi = sumphi(ixr,iyr,nd1,ndn,capphi,dss,nncb,  &
           nnsbd,asym)
       tdh  = pad*(wphi-temphi)/capphi(1) + pawf*tdh
       IF (.NOT.surf) dphi = wphi
       800 IF (surf) CALL mbgaw (boxl,dphi,ws,paw,paf1,paf2,q,q1,q2,j,kc,kc1,  &
           kc2)
       
       dss(i,j) = tdh
     END DO
   END DO
   CALL mbgate (ntote,dphite,nnsbd,ywte,q,q1,q2,kte,kte1,kte2)
   CALL mbgae (ajjl,in17,a,f,df,f1,df1,f2,df2,q,q1,q2,mood)
   CALL bug (nhcont,3000,njj,30)
   CALL bug (nhdss ,3000,dss,4)
 END DO
 CALL CLOSE (in17,1)
 RETURN
END SUBROUTINE mbdpdh
