SUBROUTINE q4bmgd (dshp,gpth,bgpdt,gpnorm,phi,bmatrx)
     
!     THIS ROUTINE ASSEMBLES PORTIONS OF B-MATRIX FOR QUAD4
 
 
 DOUBLE PRECISION, INTENT(IN)             :: dshp(1)
 DOUBLE PRECISION, INTENT(IN)             :: gpth(1)
 REAL, INTENT(IN OUT)                     :: bgpdt(4,1)
 REAL, INTENT(IN)                         :: gpnorm(4,1)
 DOUBLE PRECISION, INTENT(IN)             :: phi(9)
 DOUBLE PRECISION, INTENT(OUT)            :: bmatrx(1)
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth,badj
 INTEGER :: rowflg
 
 DOUBLE PRECISION :: psitrn(9),bbar(120),atrans(6),  &
      deriv,thick,hzta,term,detj, uev,unv,anglei,edgel,edgshr,bb1,bb2,bb3,  &
     bsbar1(6),bsbar(48),tee(9)
 COMMON /q4dt  /  detj,hzta,psitrn,nnode,badj,n1
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /q4comd/  anglei(4),edgshr(3,4),edgel(4),unv(3,4),  &
     uev(3,4),rowflg,iorder(4)
!*****
!     INITIALIZE
!*****
 ndof =nnode*6
 ndof3=nnode*3
 nd2=ndof*2
 nd3=ndof*3
 nd4=ndof*4
 nd5=ndof*5
 nd6=ndof*6
 
!     SET THE SIZE OF B-MATRIX BASED ON THE ROW FLAG.
!     ROWFLG = 1      OUT OF PLANE SHEAR (LAST 2 ROWS)   ND2
!     ROWFLG = 2      IN-PLANE SHEAR (THIRD  ROW)        NDOF
!     ROWFLG = 3      THE FIRST SIX (THREE) ROWS         ND6 (ND3)
 
 nn = nd6
 IF (norpth) nn = nd3
 IF (rowflg == 1) nn = nd2
 IF (rowflg == 2) nn = nd2
!*****
!     SET UP TERMS TO BE FILLED IN B-MATRIX
!*****
 DO  k=1,nnode
   kpoint=6*(k-1)
   thick =gpth(k)
   
!     COMPUTE THE TERMS WHICH GO IN THE FIRST 6(3) ROWS.
   
   IF (rowflg == 1) GO TO 20
   atrans(1)=-psitrn(2)*gpnorm(4,k)+psitrn(3)*gpnorm(3,k)
   atrans(2)= psitrn(1)*gpnorm(4,k)-psitrn(3)*gpnorm(2,k)
   atrans(3)=-psitrn(1)*gpnorm(3,k)+psitrn(2)*gpnorm(2,k)
   atrans(4)=-psitrn(5)*gpnorm(4,k)+psitrn(6)*gpnorm(3,k)
   atrans(5)= psitrn(4)*gpnorm(4,k)-psitrn(6)*gpnorm(2,k)
   atrans(6)=-psitrn(4)*gpnorm(3,k)+psitrn(5)*gpnorm(2,k)
   
   DO  i=1,2
     ipoint=nd3*(i-1)
     itot =ipoint+kpoint
     deriv=dshp(n1*(i-1)+k)
     bbar(     1+itot)=deriv*psitrn(1)
     bbar(     2+itot)=deriv*psitrn(2)
     bbar(     3+itot)=deriv*psitrn(3)
     bbar(ndof+1+itot)=deriv*psitrn(4)
     bbar(ndof+2+itot)=deriv*psitrn(5)
     bbar(ndof+3+itot)=deriv*psitrn(6)
     term=hzta*thick*deriv
     bbar(     4+itot)=term*atrans(1)
     bbar(     5+itot)=term*atrans(2)
     bbar(     6+itot)=term*atrans(3)
     bbar(ndof+4+itot)=term*atrans(4)
     bbar(ndof+5+itot)=term*atrans(5)
     bbar(ndof+6+itot)=term*atrans(6)
   END DO
   CYCLE
   
!     COMPUTE THE TERMS WHICH GO IN THE LAST 2 ROWS.
   
   20 IF (.NOT.bendng) RETURN
   tee(1)= 0.0D0
   tee(2)=-gpnorm(4,k)
   tee(3)= gpnorm(3,k)
   tee(4)=-tee(2)
   tee(5)= 0.0D0
   tee(6)=-gpnorm(2,k)
   tee(7)=-tee(3)
   tee(8)=-tee(6)
   tee(9)= 0.0D0
   
   kp1=kpoint*2
   kp2=kp1+7
   j=iorder(k)
   i=j-1
   IF (i == 0) i=4
   
   ib=0
   30 ib=ib+1
   bb1=-unv(ib,j)*edgshr(1,j)/edgel(j) +unv(ib,i)*edgshr(1,i)/edgel(i)
   bb2=-unv(ib,j)*edgshr(2,j)/edgel(j) +unv(ib,i)*edgshr(2,i)/edgel(i)
   bb3=-unv(ib,j)*edgshr(3,j)/edgel(j) +unv(ib,i)*edgshr(3,i)/edgel(i)
   bsbar(kp1+ib  )=psitrn(1)*bb1+psitrn(2)*bb2+psitrn(3)*bb3
   bsbar(kp1+ib+3)=psitrn(4)*bb1+psitrn(5)*bb2+psitrn(6)*bb3
   IF (ib < 3) GO TO 30
   
   ib=0
   40 ib=ib+1
   bb1=-(uev(ib,j)*edgshr(1,j)+uev(ib,i)*edgshr(1,i))*0.5D0
   bb2=-(uev(ib,j)*edgshr(2,j)+uev(ib,i)*edgshr(2,i))*0.5D0
   bb3=-(uev(ib,j)*edgshr(3,j)+uev(ib,i)*edgshr(3,i))*0.5D0
   bsbar1(ib  )=psitrn(1)*bb1+psitrn(2)*bb2+psitrn(3)*bb3
   bsbar1(ib+3)=psitrn(4)*bb1+psitrn(5)*bb2+psitrn(6)*bb3
   IF (ib < 3) GO TO 40
   CALL gmmatd (bsbar1,2,3,0,tee,3,3,0,bsbar(kp2))
   
!*****
!     FILL IN B-MATRIX FOR THE NORMAL PATH
!*****
   
 END DO
 IF (.NOT.norpth) GO TO 200
 SELECT CASE ( rowflg )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 100
 END SELECT
 
!     ROWFLG = 3       FIRST THREE ROWS
 
 100 DO  kbar=1,ndof
   bmatrx(kbar     )=phi(1)*bbar(kbar     )+phi(2)*bbar(kbar+nd3 )
   bmatrx(kbar+ndof)=phi(4)*bbar(kbar+ndof)+phi(5)*bbar(kbar+nd4 )
   bmatrx(kbar+nd2 )=phi(4)*bbar(kbar     )+phi(1)*bbar(kbar+ndof)  &
       +phi(5)*bbar(kbar+nd3 )+phi(2)*bbar(kbar+nd4 )
 END DO
 GO TO 300
 
!     ROWFLG = 2       IN-PLANE SHEAR (3RD ROW)
 
 120 DO  kbar=1,ndof
   bmatrx(kbar)=phi(4)*bbar(kbar    )+phi(1)*bbar(kbar+ndof)  &
       +phi(5)*bbar(kbar+nd3)+phi(2)*bbar(kbar+nd4 )
 END DO
 GO TO 300
 
!     ROWFLG = 1       OUT-OF-PLANE SHEAR (LAST 2 ROWS)
 
 140 DO  kbar=1,ndof
   ibar=((kbar-1)/3)*3+kbar
   bmatrx(kbar+ndof)=bsbar(ibar  )
   bmatrx(kbar     )=bsbar(ibar+3)
 END DO
 GO TO 300
 
!*****
!     FILL IN B-MATRIX FOR THE MIDI PATH
!*****
 
 200 DO  iji=1,nn
   bmatrx(iji)=0.0D0
 END DO
 SELECT CASE ( rowflg )
   CASE (    1)
     GO TO 280
   CASE (    2)
     GO TO 260
   CASE (    3)
     GO TO 220
 END SELECT
 
!     ROWFLG = 3       FIRST SIX ROWS
 
 220 IF (.NOT.membrn) GO TO 240
 DO  ka=1,nnode
   kk=(ka-1)*6
   DO  m=1,3
     bmatrx(m+kk     )=phi(1)*bbar(m+kk     )+phi(2)*bbar(m+kk+nd3)
     bmatrx(m+kk+ndof)=phi(4)*bbar(m+kk+ndof)+phi(5)*bbar(m+kk+nd4)
   END DO
 END DO
 
 240 IF (.NOT.bendng) GO TO 300
 DO  ka=1,nnode
   kk=(ka-1)*6
   DO  n=4,6
     bmatrx(n+kk+nd3)=phi(1)*bbar(n+kk     )+phi(2)*bbar(n+kk+nd3)
     bmatrx(n+kk+nd4)=phi(4)*bbar(n+kk+ndof)+phi(5)*bbar(n+kk+nd4)
   END DO
 END DO
 GO TO 300
 
!     ROWFLG = 2       IN-PLANE SHEAR (3RD AND 6TH ROWS)
 
 260 DO  ka=1,nnode
   kk=(ka-1)*6
   DO  m=1,3
     n=3+m
     bmatrx(m+kk     )=phi(4)*bbar(m+kk    )+phi(1)*bbar(m+kk+ndof)  &
         +phi(5)*bbar(m+kk+nd3)+phi(2)*bbar(m+kk+nd4)
     bmatrx(n+kk+ndof)=phi(4)*bbar(n+kk    )+phi(1)*bbar(n+kk+ndof)  &
         +phi(5)*bbar(n+kk+nd3)+phi(2)*bbar(n+kk+nd4)
   END DO
 END DO
 GO TO 300
 
!     ROWFLG = 1       OUT-OF-PLANE SHEAR (LAST 2 ROWS)
 
 280 DO  ka=1,nnode
   kk=(ka-1)*6
   DO  m=1,3
     n=3+m
     kkk=kk*2
     bmatrx(m+kk+ndof)=bsbar(m+kkk  )
     bmatrx(n+kk+ndof)=bsbar(m+kkk+6)
     bmatrx(m+kk     )=bsbar(n+kkk  )
     bmatrx(n+kk     )=bsbar(n+kkk+6)
   END DO
 END DO
 
 300 CONTINUE
 RETURN
END SUBROUTINE q4bmgd
