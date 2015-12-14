SUBROUTINE ortho (u,v,x1,x2,x3,x4,x5,nz,ibuf1,ibuf2,ibuf3,ibuf4)
     
!     ORTHO WILL ORTHOGONALIZE THE CURRENT ITERANT WITH RESPECT TO
!     THE PREVIOUSLY EXTRACTED EIGENVECTORS
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: u(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v(1)
 DOUBLE PRECISION, INTENT(IN)             :: x1(1)
 DOUBLE PRECISION, INTENT(OUT)            :: x2(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: x3(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: x4(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: x5(1)
 INTEGER, INTENT(IN OUT)                  :: nz
 INTEGER, INTENT(IN OUT)                  :: ibuf1(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf2(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf3(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf4(1)
 INTEGER :: filem     ,filek    ,fileb    ,filelm   ,  &
     filevc    ,sr0fil   ,sr5fil   ,REAL     , sub(2)    ,  &
     sqr       ,cdp
 DOUBLE PRECISION :: pj(2)    , const1(2) ,const2(2),alpha(2),beta(2)
 COMMON /cinvpx/  filek(7)  ,filem(7) ,fileb(7),filelm(7) ,  &
     filevc(7),dmpfil    ,scrfil(10)
 COMMON /cinvxx/  dum(17)   ,REAL     ,xxxx     ,northo
 COMMON /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp       ,csp      ,cdp      ,sqr
 EQUIVALENCE      (sr0fil,scrfil(10)) ,(sr5fil,scrfil(5)) , (ncol,filek(2))
 DATA    sub   /  4HORTH ,4HO   /
 
 ncol2 = ncol  + ncol
 ncol4 = ncol2 + ncol2
 IF (fileb(1) == 0) GO TO 5
 CALL cmtimu (v,x1,0,ibuf4)
 CALL cmtimu (u,x2,fileb,ibuf4)
 GO TO 7
 5 DO  i = 1,ncol2
   x2(i) = 0.d0
 END DO
 7 CONTINUE
 CALL cmtimu (u,x3,0,ibuf4)
 CALL sswtch (12,l7)
 const1(1) = 1.0D0
 const1(2) = 0.
 const2(1) =-1.0D0
 const2(2) = 0.
 CALL csub (x1,x2,x2,const1,const2)
 
!     REWIND EIGENVALUE AND EIGENVECTOR FILES
 
 ifile = filelm(1)
 CALL OPEN (*1000,ifile,ibuf1,rdrew)
 ifile = filevc(1)
 CALL OPEN (*1000,ifile,ibuf2,rdrew)
 ifile = sr0fil
 CALL OPEN (*1000,ifile,ibuf3,rdrew)
 DO  k = 1,northo
   
!     READ AN EIGENVALUE
   
   ifile = filelm(1)
   CALL READ (*1010,*1020,ifile,pj(1),4,1,flag)
   const1(1) = -1.d0
   const1(2) = 0.
   CALL csub (x3,x2,x5,pj,const1)
   
!     READ THE RIGHT EIGENVECTOR
   
   ifile = filevc(1)
   CALL READ (*1010,*1020,ifile,x1(1),ncol4,1,flag)
   
!     READ THE LEFT EIGENVECTOR
   
   ifile = sr0fil
   CALL READ (*1010,*1020,ifile,x4(1),ncol4,1,flag)
   
   IF (fileb(1) /= 0) GO TO 40
   
!    COMPUTE ALPHA USING REAL FORMULA
   
   CALL cx trn y (x4,x3,const1)
   GO TO 55
   40 CALL cx trn y (x4(1),x5(1),const1(1))
   55 alpha(1) = const1(1)
   alpha(2) = const1(2)
   beta(1)  = alpha(1)*pj(1) - alpha(2)*pj(2)
   beta(2)  = alpha(1)*pj(2) + alpha(2)*pj(1)
   IF (l7 == 0) GO TO 1901
   WRITE  (6,500) const1,const2,alpha
   500 FORMAT (4H num ,2D12.5,6H denom ,2D12.5,6H alpha ,2D12.5 )
   1901 CONTINUE
   DO  i = 1,ncol2,2
     u(i  ) = u(i  ) - alpha(1)*x1(i) + alpha(2)*x1(i+1)
     u(i+1) = u(i+1) - alpha(2)*x1(i) - alpha(1)*x1(i+1)
     IF (fileb(1) == 0) CYCLE
     v(i  ) = v(i  ) - beta(1)*x1(i  ) + beta(2)*x1(i+1)
     v(i+1) = v(i+1) - beta(1)*x1(i+1) - beta(2)*x1(i  )
   END DO
 END DO
 CALL CLOSE (filelm,norew)
 CALL CLOSE (filevc,norew)
 CALL CLOSE (sr0fil,norew)
 RETURN
 
 1000 no = -1
 GO TO 1500
 1010 no = -2
 GO TO 1500
 1020 no = -3
 1500 CALL mesage (no,ifile,sub)
 RETURN
END SUBROUTINE ortho
