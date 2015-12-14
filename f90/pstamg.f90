SUBROUTINE pstamg (INPUT,ajjl,skj)
     
!     DRIVER FOR PISTON THEORY
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: ajjl
 INTEGER, INTENT(IN OUT)                  :: skj
 INTEGER :: sysbuf, NAME(2),iz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,iout
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /pstonc/ nnj,nmach,nthry,nthick,nalpha,nxis,ntaus,nstrip, seclam
 EQUIVALENCE     (z(1),iz(1))
 DATA    nhacpt, nhcomm,NAME /4HACPT,4HCOMM,4HPSTA,4HMG  /
 
 icore = korsz(iz) - 4*sysbuf
 
!     BRING IN DATA AND ALLOCATE CORE
 
 CALL fread (INPUT,nnj,9,0)
 idel = 1
 ib   = idel + nstrip
 ica  = ib  + nstrip
 ipalp= ica + nstrip
 
!     READ FIXED ARRAYS
 
 nw = 3*nstrip
 CALL fread (INPUT,z,nw,0)
 
!     READ ALPHA ARRAY AND STUFF  AT END (INTEGRALS OR TAUS)
 
 iend = 0
 DO  i = 1,nmach
   CALL fread (INPUT,rm,1,0)
   IF (rm /= fmach) GO TO 10
   iend = 1
   CALL fread (INPUT,z(ipalp),nalpha,0)
   CYCLE
   10 CALL fread (INPUT,z,-nalpha,0)
 END DO
 IF (iend == 0) GO TO 50
 ipt = ipalp + nalpha
 CALL READ (*30,*30,INPUT,z(ipt),icore,1,n)
 30 nt = ipt + n
 CALL bug (nhacpt,30,z,nt)
 CALL bug (nhcomm,30,nnj,9)
 
!     OUTPUT SKJ
 
 rm  = 1.0
 iti = 1
 ito = 3
 ii  = isk
 nsk = nsk + 1
 nn  = nsk
 DO  i = 1,nnj
   CALL pack (rm,skj,tskj)
   ii  = ii + 1
   IF (i == nnj) CYCLE
   nn  = nn + 1
 END DO
 isk = ii
 nsk = nn
 iti = 3
 ito = 3
 CALL psta (z(idel),z(ib),z(ica),z(ipalp),z(ipt),ajjl)
 nrow = nrow + nnj
 GO TO 70
 
!     ERROR MESSAGE
 
 50 WRITE  (iout,60) ufm,fmach
 60 FORMAT (a23,' 2428, MACH NUMBER ',f10.5,' WAS NOT FOUND IN ',  &
     'PISTON THEORY ALPHA ARRAY.')
 CALL mesage (-61,0,NAME)
 
 70 RETURN
END SUBROUTINE pstamg
