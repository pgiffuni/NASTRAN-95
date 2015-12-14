SUBROUTINE detdet(deta,ipowr1,p,sml1,oldd,iold)
     
 
 DOUBLE PRECISION, INTENT(OUT)            :: deta(1)
 INTEGER, INTENT(OUT)                     :: ipowr1(1)
 DOUBLE PRECISION, INTENT(IN)             :: p(1)
 REAL, INTENT(OUT)                        :: sml1
 DOUBLE PRECISION, INTENT(IN)             :: oldd
 INTEGER, INTENT(IN)                      :: iold
 
 
 DOUBLE PRECISION :: prod,det,core,sml2, det1,sml21
 
 DOUBLE PRECISION :: mindd
 
 INTEGER :: fa,fl,fc,sr1,sr2,sr3,fa1,fl1,fc1,sr11,sr21
 INTEGER :: option,sdet
 INTEGER :: otpe,b,c,r,sr31
 
 COMMON /machin/mach
 COMMON /regean/im(21),ia,in(4),lc,in1(2),mz,in2(2),rminr,in3(2),  &
     nevm,il1,il2,nfound,lama,ibuck,nsym
 COMMON /detmx /isew(33),ipaav,isew1(5),ndcmp,isew2(20),npole ,ising
 COMMON /dcompx/fa(7),fl(7),fc(7),sr1,sr2,sr3,det,ipowr,nz ,sml2
 COMMON /zzzzzz /core(1)
 COMMON /reigkr/option
 COMMON /sfact /fa1(7)   ,fc1(7)   ,fl1(7)   ,sr11     ,sr21  &
     ,nz1      ,det1(2)  ,ipowra   ,sr31     ,mindd  &
     ,ichol    ,b        ,c        ,r        ,sml21
 COMMON /system/ksystm(65)
 
 EQUIVALENCE ( ksystm( 2) , otpe   )
 
 DATA sdet /4HSDET/
 
! ----------------------------------------------------------------------
 
 CALL sswtch (7, iprt)
 isave = in(4)
 in(4) = il2
 il2 = isave
 nzz = (korsz(core)/2)*2 -lc
 ndcmp = ndcmp+1
 IF(option == sdet) GO TO 5
 fa(1)=ia
 CALL rdtrl(fa)
 
!     SET UP FOR UNSYMMETRIC
 
 nz = nzz
 
 
!     PUT IN TO PREVENT REWRITE
 
!     FA1(1) = -FA1(1)
 
!WKBD 10/94 SPR94011 FA(1) = -FA(1)
 fl(1)=in(1)
 fc(1)=in(2)
 DO  i=2,5
   fl(i)=fa(i)
   fc(i)=fa(i)
 END DO
 sr1 = in(3)
 sr2 = in(4)
 sr3 = il1
 CALL decomp(*60,core,core,core)
!WKBD 10/94 SPR94011     FC(1) = SR2
 CALL wrttrl(fc)
 GO TO 14
 
!     SET UP FOR SYMMETRIC DECOMPOSITION
 
 5 fa1(1) = ia
 CALL rdtrl(fa1)
 fl1(1) = in(1)
 fc1(1) = in(4)
 ichol = 0
 IF(ndcmp == 1) b=0
 nz1 = nzz
 DO  i = 2,5
   fl1(i) = fa1(i)
   fc1(i) = fa1(i)
 END DO
 sr11= in(3)
 sr21 = in(2)
 sr31 = il1
 IF(mach == 4 .OR. mach == 12) fl1(5) = 1
 CALL sdcomp(*60,core,core,core)
 fc1(5) = fl1(5)
 CALL wrttrl(fc1)
 ipowr=ipowra
 det = det1(1)
 sml1= sml21
 14 prod = 1.0D0
 IF(iprt == 0) GO TO 15
 WRITE(otpe,99) p(1),det,ipowr
 99 FORMAT(2D16.7,i8)
 15 CONTINUE
 iprod = 0
 IF( mz  == 0) GO TO 12
 ii =  IABS(mz)
 DO   i=1,ii
   prod = prod* p(1)
   CALL detm6(prod,iprod)
 END DO
 
!     TAKE OUT  POLE AT  RMINR
 
 12 IF (npole == 0) GO TO 20
 DO    i = 1,npole
   prod = prod*(p(1)- rminr)
   CALL detm6( prod,iprod)
 END DO
 20 IF(nfound == 0) GO TO 40
 DO  i=1,nfound
   ii = ipaav +i
   IF(p(1) == core(ii)) GO TO 70
   prod = prod*(p(1)-core(ii))
   CALL detm6(prod,iprod)
 END DO
 40 deta(1) = det/prod
 sml1= sml2
 ipowr1(1)= ipowr-iprod
 CALL detm6(deta(1),ipowr1(1))
 50 IF(iprt == 0) GO TO 51
 WRITE(otpe,99) p(1),deta(1),ipowr1(1)
 51 RETURN
 60 deta(1) = 0.0D0
 ipowr1(1)=1
 sml1 = 1.0E-8
 ising = ising+1
 isave = in(4)
 in(4) = il2
 il2 = isave
 GO TO 50
 
!     SET DK = DK-1
 
 70 deta(1) = oldd
 sml1 = sml2
 ipowr1(1) = iold
 GO TO 50
END SUBROUTINE detdet
