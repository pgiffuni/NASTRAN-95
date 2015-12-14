SUBROUTINE alg
     
    !     THIS IS THE DRIVER SUBROUTINE FOR THE ALG MODULE
 
    INTEGER :: apress,atemp,strml,pgeom,NAME(2),sysbuf, title1(18),wd(2),algdb
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm,uim
    COMMON /BLANK / apress,atemp,strml,pgeom,iprtk,ifail,SIGN,zorign,  &
        fxcoor,fycoor,fzcoor
    COMMON /system/ sysbuf,nout
    COMMON /algino/ iscr3,algdb
    COMMON /udstr2/ nbldes,stag(21),chordd(21)
    COMMON /ud3prt/ iprtc,istrml,ipgeom
    COMMON /zzzzzz/ iz(1)
    COMMON /contrl/ nanal,naero,narbit,log1,log2,log3,log4,log5,log6
    DATA    NAME  / 4HALG ,4H    /
    DATA    wd    / 2HNO  ,2HAN  /
    DATA    iscr1 , iscr2 / 301,302 /
 
    iscr3  = 303
    iscr4  = 304
    istrml = strml
    ipgeom = pgeom
    IF (ipgeom == 3) ipgeom = 1
    iprtc = iprtk
    nz    = korsz(iz)
    ibuf1 = nz - sysbuf + 1
    ibuf2 = ibuf1 - sysbuf
    ibuf3 = ibuf2 - sysbuf
    IF (3*sysbuf > nz) CALL mesage (-8,0,NAME)
    CALL algpr (ierr)
    IF (ierr < 0) GO TO 400
    algdb = iscr1
    IF (ierr == 1) algdb = iscr2
    log1  = algdb
    log2  = nout
    log3  = 7
    log4  = algdb
    log5  = iscr4
    log6  = 9
    CALL gopen (log1,iz(ibuf1),0)
    CALL fread (log1,title1,18,1)
    CALL fread (log1,nanal,1,0)
    CALL fread (log1,naero,1,1)
    narbit = 0
    IF (iprtc == 1) WRITE (log2,20) title1,nanal,wd(naero+1)
    IF (iprtc == 0) WRITE (log2,40) uim
20  FORMAT (1H1,/40X,'ALG module - compressor design - control section' ,&
        /40X,48(1H*), //10X,'TITLE = ',18A4, /10X,'NUMBER of analytic &
           &mealine bladerows ='    ,i3, /10X,'THERE will be ',a2,' ENTRY to &
           &the aerodynamic section '    )
40  FORMAT (a29,' - MODULE ALG ENTERED.')
 
    IF (nanal == 0) GO TO 200
    ifile = log5
    CALL OPEN (*500,log5,iz(ibuf2),1)
    CALL algan
    CALL CLOSE (log5,1)
200 IF (naero == 0) GO TO 300
    ifile = log5
    CALL OPEN (*500,log5,iz(ibuf2),0)
    ifile = iscr3
    CALL OPEN (*500,iscr3,iz(ibuf3),1)
    CALL algar
    CALL CLOSE (iscr3,1)
    CALL CLOSE (log5,1)
300 CALL CLOSE (log1,1)
    CALL algpo (iscr3)
400 GO TO 600
500 CALL mesage(-1,ifile,NAME)
 
600 RETURN
END SUBROUTINE alg
