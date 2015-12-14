SUBROUTINE detm5
     
!     WRITES EIGENVALUE SUMMARY FOR DETERMINANT METHOD
 
 DOUBLE PRECISION :: psave(1),det(1), ps(1)
 INTEGER :: ipdet(8),sysbuf
 DIMENSION core(5)
 
 COMMON /condas/ consts(5)
 COMMON /detmx/ p(32),n2ev,ipsav,ips,idet,ipdeta,prec,nstart,  &
     ndcmp,ic,nsmove,iterm,is,nd,iadd,sml1,ipdetx(4),ipdet1(4),  &
     ifail,k,fact1,iffnd,nfail
 COMMON /regean/ im(26),lcore,rmax,rmin,mz, nev, epsi, rminr, NE,  &
     nit, nevm, scr6, ipout, nfound, lama
 COMMON /zzzzzz/ psave
 COMMON /system/sysbuf
 
 EQUIVALENCE  (psave(1),ps(1),det(1),ipdet(1),core(1))
 EQUIVALENCE ( consts(2) , tphi   )
 
! ----------------------------------------------------------------------
 
 nz = korsz(psave) -lcore-sysbuf
 CALL gopen(ipout,ipdet(nz+1),1)
 ipdet(1) = 1
 ipdet(2) = nfound
 IF(mz > 0) ipdet(2) = ipdet(2) +mz
 ipdet(3) = nstart
 ipdet(4) = ic
 ipdet(5) = nsmove
 ipdet(6) = ndcmp
 ipdet(7) = nfail
 ipdet(8) = iterm
 DO  i=9,12
   ipdet(i) = 0
 END DO
 CALL WRITE(ipout,ipdet(1),12,0)
 IF(ndcmp == 0) GO TO 61
 n2ev2 = iadd+nd
 DO  i=1,n2ev2
   nnd = i+idet
   nnp = i+ips
   nni = i+ipdeta
   
!     PUT UUT STRRTING POINT SUMMARY
   
   ipdet(1) = i
   core(2) = psave(nnp)
   core(3) = SQRT(ABS(core(2)))
   core(4) = core(3)/tphi
   core(5) = psave(nnd)
   ipdet(6) = ipdet(nni)
   
!     SCALE DETERMINANTE FOR PRETTY PRINT
   
   IF(core(5) == 0.0) GO TO 50
   20 IF(ABS(core(5)) >= 10.0) GO TO 40
   30 IF(ABS(core(5)) >= 1.0) GO TO 50
   core(5) = core(5)*10.0
   ipdet(6) = ipdet(6)-1
   GO TO 30
   40 core(5) = core(5)*0.1
   ipdet(6) = ipdet(6)+1
   GO TO 20
   50 CALL WRITE(ipout,core(1),6,0)
 END DO
 61 CONTINUE
 CALL WRITE(ipout,core(1),0,1)
 CALL CLOSE(ipout,1)
 RETURN
END SUBROUTINE detm5
