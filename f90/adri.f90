SUBROUTINE adri (fl,nfreq,ncore,qhhl,scr2,scr1,scr3,scr4,nrow,  &
    ncol,nogo)
     
 

    REAL, INTENT(OUT)           :: fl(1)
    INTEGER, INTENT(IN)         :: nfreq ,ncore ,qhhl  ,ncol
    INTEGER, INTENT(IN OUT)     :: scr1  ,scr2  ,scr3  ,scr4
    INTEGER, INTENT(OUT)        :: nrow  ,nogo
    INTEGER                     :: trl(7),out

    DIMENSION  mcb(7),name(2)

    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm
    COMMON /blank / bov   ,rm
    COMMON /condas/ pi    ,twopi
    COMMON /system/ isys  ,out   ,dum(52),iprec
    COMMON /unpakx/ iout  ,inn   ,nnn   ,incr1
    COMMON /packx / iti   ,ito   ,ii    ,nn    ,incr
    COMMON /type  / p(2)  ,iwc(4)

    DATA    nhfrdi, NAME /4HFRDI,4HADRI,4H    /
 
    ibuf1 = ncore - isys
    ibuf2 = ibuf1 - isys
    nrow  = 0
    incr  = 1
    incr1 = 1
    ii    = 1
    inn   = 1
    mcb(1)= qhhl
    CALL rdtrl (mcb)
    IF (mcb(1) < 0) GO TO 250
    nrow  = mcb(3)
    CALL OPEN (*250,qhhl,fl(ibuf2),0)
    CALL gopen (scr1,fl(ibuf1),1)
    CALL READ (*220,*220,qhhl,fl(1),-2,0,flag)
    CALL READ (*220,*220,qhhl,ncol,1,0,flag)
    CALL READ (*220,*220,qhhl,n,1,0,flag)
    n    = n + n
    ni   = (mcb(2)/ncol)*2
    ni   = MIN0(ni,n)
    nnn  = nrow
    nn   = ncol*nrow
    iti  = 3
    ito  = iti
    iout = iti
    nwc  = iwc(iti)
    CALL makmcb (trl,scr1,nn,mcb(4),ito)
 
    !     MAKE   DEPENDENT FREQ LIST
 
    ipd  = 1
    nl   = 2*nfreq
    n    = nfreq + 1
    ipi  = ipd + nl
    DO  i = 1,nfreq
        fl(nl  ) = fl(n-i)*twopi*bov
        fl(nl-1) = 0.0
        nl   = nl -2
    END DO
 
    !     MAKE INDEPENDENT FREQ LIST
 
    CALL READ (*220,*220,qhhl,fl(ipi),ni,1,flag)
 
    !     FIND M"S CLOSEST TO RM
 
    icp = ipi + ni
    rmi = 1.e20
    rms = 0.0
    DO  i = 1,ni,2
        rmx = ABS(fl(ipi+i-1) - rm)
        rmi = AMIN1(rmi,rmx)
        IF (rmx > rmi) CYCLE
        rms = fl(ipi+i-1)
    END DO
    rmi = rms
 
    !     DO ALL K"S ASSOCIATED WITH RMI
 
    k = 0
    DO  i = 1,ni,2
        IF (fl(ipi+i-1) == rmi) GO TO 120
   
        !     SKIP MATRIX
   
        CALL skprec (qhhl,ncol)
        CYCLE
   
        !     MAKE MATRIX INTO COLUMN
   
120     fl(ipi+k+1) = fl(ipi+i)
        k  = k + 2
        ji = icp
        n  = nrow*nwc
        DO  j = 1,ncol
            CALL unpack (*131,qhhl,fl(ji))
            GO TO 135
131         CALL zeroc (fl(ji),n)
135         ji = ji + n
        END DO
   
        !     DIVIDE IMAG PART OF QHHL BY FREQUENCY
   
        jj = icp + 1
        kk = ji  - 1
        DO  j = jj,kk,2
            fl(j) = fl(j)/fl(ipi+i)
        END DO
        CALL pack (fl(icp),scr1,trl)
    END DO
    CALL CLOSE (qhhl,1)
    CALL CLOSE (scr1,1)
    CALL wrttrl (trl)
    CALL bug (nhfrdi,150,k ,1)
    CALL bug (nhfrdi,150,nfreq,1)
    CALL bug (nhfrdi,150,fl(1),icp)
 
    !     SETUP TO CALL MINTRP
 
    ni   = k/2
    nogo = 0
    nc   = ncore - icp
    CALL dmpfil (-scr1,fl(icp),nc)
    im   = 0
    ik   = 1
    CALL mintrp (ni,fl(ipi),nfreq,fl(ipd),-1,im,ik,0.0,scr1,scr2,  &
        scr3,scr4,fl(icp),nc,nogo,iprec)
    IF (nogo == 1) GO TO 200
    CALL dmpfil (-scr2,fl(icp),nc)
    RETURN
 
200 WRITE  (out,210) ufm
210 FORMAT (a23,' 2271, INTERPOLATION MATRIX IS SINGULAR')
    !IBMR 6/93  GO TO 240                                                 !*
    GO TO 240
220 CALL mesage (3,qhhl,NAME)
240 nogo = 1
250 CALL CLOSE (qhhl,1)

    RETURN
END SUBROUTINE adri
