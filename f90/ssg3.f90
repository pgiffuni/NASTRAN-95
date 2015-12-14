SUBROUTINE ssg3
     
!     DMAP FOR STATIC SOLUTION GENERATOR 3
 
!     SSG3   LLL,KLL,PL,LOO,KOOB,PO /ULV,UOV,RULV,RUOV/ V,N,OMIT/
!            V,Y,IRES/V,N,SKIP/V,N,EPSI $
 
 INTEGER :: lll,kll,pl,loo,po,ulv,uov,sr1,sr2,omit,rulv,ruov
 COMMON /BLANK/ omit,ires,nskip,epsi
 DATA    lll  , kll,  pl, loo, koob, po, ulv, uov, rulv, ruov /  &
     101  , 102, 103, 104, 105, 106, 201, 202, 203 , 204  /
 DATA    sr1  , sr2  / 301  , 302  /
 
 CALL ssg3a (kll,lll,pl,ulv,sr1,sr2,0,rulv)
 IF (omit >= 0) CALL ssg3a (koob,loo,po,uov,sr1,sr2,0,ruov)
 RETURN
END SUBROUTINE ssg3
