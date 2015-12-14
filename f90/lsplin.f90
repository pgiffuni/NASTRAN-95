SUBROUTINE lsplin(ni,xyi,nd,xyd,ky,kd,kt,dz,dtx,dty,dtor,  &
        g,ncore,isng)
     
 INTEGER, INTENT(IN)                      :: ni
 REAL, INTENT(IN)                         :: xyi(1)
 INTEGER, INTENT(IN)                      :: nd
 REAL, INTENT(IN)                         :: xyd(1)
 INTEGER, INTENT(IN)                      :: ky
 INTEGER, INTENT(IN)                      :: kd
 INTEGER, INTENT(IN OUT)                  :: kt
 REAL, INTENT(IN)                         :: dz
 REAL, INTENT(IN)                         :: dtx
 REAL, INTENT(IN)                         :: dty
 REAL, INTENT(IN)                         :: dtor
 REAL, INTENT(OUT)                        :: g(1)
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(OUT)                     :: isng
 LOGICAL :: ikt
 LOGICAL :: bnone,bone
 LOGICAL :: spec
 LOGICAL :: oxr,oyr
 LOGICAL :: both
 LOGICAL :: stwo,sone,nethr
 INTEGER :: size
 
 DIMENSION  NAME(2)
 DATA NAME/4HLSPL,4HIN  /
 
 spec = .false.
 bone = .false.
 bnone = .false.
 sone = .false.
 stwo = .false.
 oxr = .false.
 oyr = .false.
 both = .false.
 ikt = .false.
 nethr = .true.
 ey = FLOAT(ky)
 
!     KY EFFECTS RIGID BODY ROWS AND COLUMNS OF A AND ROWS OF B
 
!     DTX AND DTY EFFECT ROWS AND COLUMNS OF A AND ROWS OF B
 
!     KD EFFECTS COLUMNS OF B
 
!     SPEC GUARDS AGAINST SINGULAR MATRIX FOR ROTATIONS WITHOUT Y
!     CONSTRAINED
 
 IF(kd == 0) bnone = .true.
 IF(kd == 1) bone = .true.
 IF(ky == -1) sone = .true.
 IF(ky == 1) stwo = .true.
 IF(kt == 1) ikt = .true.
 IF(dtx < 0.0) oxr = .true.
 IF(dty < 0.0) oyr = .true.
 IF(oxr.AND.oyr) both = .true.
 IF(.NOT.oxr.AND..NOT.oyr) nethr = .false.
 dtor2 =  dtor/2.0
 nsc = 3
 IF(ky == 1) nsc = 2
 IF(ky == -1) nsc = 1
 size = 3
 IF(oxr) size = size-1
 IF(oyr) size = size-1
 IF(oyr.AND. ky > -1) GO TO 5
 GO TO 7
 5 temp = xyi(1)
 spec = .true.
 nii = 2*ni
 DO  i = 1,nii,2
   IF(xyi(i) /= temp) spec = .false.
 END DO
 7 CONTINUE
 nca = size*ni + nsc
 IF(spec) nca = nca -1
 nca2= 2*nca
 ncb = (kd+1) * nd
 ncc = size*ni
 
!     CORE NEEDED
!                A         G        INVERS
 needed =nca*nca+  ncb*ncc  +    3*nca
!                                   B
 IF(ikt) needed = needed +  ncb*nca
!                                        C
 IF(.NOT.ikt) needed = needed +  ncc*nca
 IF(needed > ncore) CALL mesage(-8,0,NAME)
 is= ncore -3*nca-1
 ig= 1
 
!     IF KT = 1  COMPUTE  B THEN A  THEN  C IN THE SPACE OF A
 
!     IF KT = 0 COMPUTE  C THEN A  THEN  B IN THE SPACE OF A
 
 IF(ikt) GO TO 100
 
!     FILL IN C MATRIX
 
 ic = ncb * ncc
 mp = ic +1
 10 DO  i = 1,ncc
   DO  j = 1,nca
     ic = ic +1
     g(ic) = 0.0
     IF(i == j) g(ic) = 1.0
   END DO
 END DO
 IF(ikt) GO TO 300
 nc = ncc
 ia =  ic
 GO TO 200
 
!     B MATRIX
 
 100 ib =  ncb * ncc
 mp = ib + 1
 110 nj = 2*nd
 nii= 2*ni
 DO  j=1,nj,2
   DO  i=1,nii,2
     ym = xyd(j+1) - xyi(i+1)
     aym= ABS(ym)
     aymd = aym*dtor2
     yp =(xyd(j+1) + xyi(i+1))
     ayp= ABS(yp) * ey
     aypd = ayp*dtor2
     ib= ib +1
     g(ib) =    aym**3/12.0  -  xyd(j)*xyi(i)*aymd  &
         + ayp**3/12.0 - xyd(j)*xyi(i)*aypd
     IF(bnone) GO TO 111
     g(ib+nca) = aym*ym/4.0  + ayp*yp/4.0
     IF(bone) GO TO 111
     g(ib+nca2)=  xyi(i)*aymd  + xyi(i)*aypd
     111 IF(both) CYCLE
     IF(oxr) GO TO 113
     ib=ib+1
     g(ib) =    -aym*ym/4.0  + ayp*yp/4.0
     IF(bnone) GO TO 112
     g(ib+nca) = -aym/2.0 + ayp/2.0
     IF(bone) GO TO 112
     g(ib+nca2) = 0.0
     112 IF(oyr) CYCLE
     113 ib = ib +1
     g(ib)=   xyd(j)*aymd   + xyd(j)*aypd
     IF(bnone) CYCLE
     g(ib+nca)  =  0.0
     IF(bone) CYCLE
     g(ib+nca2) =  -aymd   - aypd
   END DO
   ib = ib +1
   IF(sone) GO TO 123
   g(ib) = 1.0
   IF(bnone) GO TO 121
   g(ib+nca) = 0.0
   IF(bone) GO TO 121
   g(ib+nca2) = 0.0
   121 IF(stwo) GO TO 122
   ib = ib +1
   g(ib) = xyd(j+1)
   IF(bnone) GO TO 122
   g(ib+nca) = 1.0
   IF(bone) GO TO 122
   g(ib+nca2) = 0.0
   122 IF(spec) GO TO 128
   ib = ib +1
   g(ib) =-xyd(j)
   IF(bnone) GO TO 128
   g(ib+nca) = 0.0
   IF(bone) GO TO 128
   g(ib+nca2) = 1.0
   GO TO 128
   123 g(ib) = xyd(j+1)
   IF(bnone) GO TO 128
   g(ib+nca) = 1.0
   IF(bone) GO TO 128
   g(ib+nca2) = 0.0
   128 ib = ib + kd*nca
 END DO
 IF(.NOT.ikt) GO TO 400
 ia = ib
 nc = ncb
 
!     A MATRIX
 
 200 nii= 2*ni
 k = ia
 
!     ZERO A
 
 ii = k+1
 ik = ii + nca*nca
 DO  i = ii,ik
   g(i) = 0.0
 END DO
 ii = 0
 DO  i = 1,nii,2
   DO  j = i,nii,2
     k = k+1
     yp   =(xyi(i+1) + xyi(j+1))
     ayp  = ABS(yp) * ey
     aypd = ayp*dtor2
     ym   = xyi(i+1) - xyi(j+1)
     aym  = ABS(ym)
     aymd = aym*dtor2
     g(k)       = aym**3/12.0 - xyi(i)*xyi(j)*aymd  &
         + ayp**3/12.0 - xyi(i)*xyi(j)*aypd
     IF(i == j) g(k) = g(k) + dz
     IF(both) CYCLE
     IF(oxr) GO TO 212
     g(k+nca)  = aym*ym/4.0  + ayp*yp/4.0
     IF(oyr) GO TO 214
     g(k+nca2)  = xyi(i)*aymd  + xyi(i)*aypd
     k = k+1
     g(k)  =    -aym*ym/4.0  + ayp*yp/4.0
     g(k+nca)  = -aym/2.0 + ayp/2.0
     IF(i == j) g(k+nca) = g(k+nca) + dtx
     g(k+nca2)  = 0.0
     k = k+1
     g(k)       = xyi(j)*aymd  + xyi(j)*aypd
     g(k+nca)   = 0.0
     g(k+nca2) = -aymd  - aypd
     IF(i == j) g(k+nca2) = g(k+nca2) + dty
     CYCLE
     212 g(k+nca) = xyi(i)*aymd  + xyi(i)*aypd
     k = k+1
     g(k) = xyi(j)*aymd  + xyi(j)*aypd
     g(k+nca) = -aymd   - aypd
     IF(i == j) g(k+nca) = g(k+nca) + dty
     CYCLE
     214 k = k+1
     g(k)  =    -aym*ym/4.0  + ayp*yp/4.0
     g(k+nca)  = -aym/2.0 + ayp/2.0
     IF(i == j) g(k+nca) = g(k+nca) + dtx
   END DO
   k = k+1
   IF(sone) GO TO 234
   g(k) = 1.0
   IF(both) GO TO 231
   g(k+nca) = 0.0
   IF(nethr) GO TO 231
   g(k+nca2) = 0.0
   231 IF(stwo) GO TO 232
   k = k+1
   g(k) = xyi(i+1)
   IF(both) GO TO 232
   IF(oxr) g(k+nca) = 0.0
   IF(oyr) g(k+nca)  = 1.0
   IF(nethr) GO TO 232
   g(k+nca) = 1.0
   g(k+nca2) = 0.0
   232 IF(spec) GO TO 238
   k = k+1
   g(k) = -xyi(i)
   IF(both) GO TO 238
   IF(oxr) g( k+nca) = 1.0
   IF(oyr) g( k+nca) = 0.0
   IF(nethr) GO TO 238
   g( k+nca) = 0.0
   g( k+nca2) = 1.0
   GO TO 238
   234 g(k) = xyi(i+1)
   IF(both) GO TO 238
   IF(oxr) g( k+nca) = 0.0
   IF(oyr) g( k+nca) = 1.0
   IF(nethr) GO TO 238
   g( k+nca) = 1.0
   g( k+nca2) = 0.0
   238 ii = ii+1
   k = k + size*ii + (size-1)*nca
 END DO
 
!     LOWER TRIANGLE IF A STORED TRANSPOSE INTO UPPER TRIANGLE
 
 k = ia
 DO  i = 1,nca
   DO  j = i,nca
     k = k+1
     kk = k + (nca-1)*(j-i)
     g(kk) = g(k)
   END DO
   k = k+i
 END DO
 
!     CALL  INVERSE   A-1 C  OR  A-1 B
 
!     REPLACE CALLS TO INVAER WITH CALLS TO INVERS.
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 isng = -1
 CALL invers(nca,g(ia+1),nca,g(mp),nc,det,isng,g(is))
 IF(isng == 2) GO TO 1000
 
!     ADJUST INDEXES TO A AND A-1 RESULT
 
 ib = ia
 icb= ib+1
 IF(.NOT.ikt) GO TO 110
 ic = ia
 icc = ic +1
 GO TO 10
 300 CALL gmmats(g(mp),ncb,nca,0,g(icc),ncc,nca,1,g(ig))
 GO TO 1000
 400 CALL gmmats(g(mp),ncc,nca,0,g(icb),ncb,nca,1,g(ig))
 1000 RETURN
END SUBROUTINE lsplin
