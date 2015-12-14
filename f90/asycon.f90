SUBROUTINE asycon
     
    !     SUBROUTINE FOR COMPUTING CONSTANT TERM IN KAPPA MINUS
 
    COMPLEX :: bsycon,ai,c1,c1test,alp,aln,arat1,arat2
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /blk2  / bsycon
    COMMON /system/ sysbuf,ibbout
    COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 
    c1   = 1.0
    pi2  = 2.0*pi
    a1   = pi2/(sps-sns)
    gam0 = sps*del - sigma
    a2   =-a1
    b1   = gam0/(sps-sns)
    s1   = sps/(dstr**2)
    s2   = sns/dstr
    c1test = 0.0
    DO  i = 1,200
        r    = i
        gamp = pi2*r + gam0
        gamn =-pi2*r + gam0
        c2p  = gamp/dstr - scrk
        c2q  = gamp/dstr + scrk
        c2n  = gamn/dstr - scrk
        c3q  = gamn/dstr + scrk
        nn   = 0
        csec = c2p*c2q
        IF (csec < 0.0) nn = 1
        t1   = gamp*s1
        t2   = s2*SQRT(ABS(csec))
        IF (c2p < 0.0 .AND. c2q < 0.0) t2 =-t2
        IF (nn == 0) alp = t1 + t2
        IF (nn == 1) alp = CMPLX(t1,t2)
        nn   = 0
        csec = c2n*c3q
        IF (csec < 0.0) nn = 1
        t1   = gamn*s1
        t2   = s2*SQRT(ABS(csec))
        IF (c2n < 0.0 .AND. c3q < 0.0) t2 =-t2
        IF (nn == 0) aln = t1 + t2
        IF (nn == 1) aln = CMPLX(t1,t2)
        arat1 = (a1*r+b1)/alp
        arat2 = (a2*r+b1)/aln
        c1    = c1*arat1*arat2
        IF (cabs((c1-c1test)/c1) < 0.0001)  GO TO 60
        c1test = c1
    END DO
    GO TO 9999
60 CONTINUE
   bsycon = c1
   RETURN
 
9999 WRITE  (ibbout,1000) ufm
1000 FORMAT (a23,' - AMG MODULE - SUBROUTINE ASYCON')
   CALL mesage (-61,0,0)

   RETURN
END SUBROUTINE asycon
