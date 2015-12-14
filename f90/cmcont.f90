SUBROUTINE cmcont
     
    !     THIS ROUTINE DEFINES THE CONNECTION ENTRIES IN TERMS OF IP
    !     NUMBERS.
 
    EXTERNAL        lshift,rshift,andf
    LOGICAL :: odd
    INTEGER :: scsfil,scmcon,buf3,buf4,ofile,scr1,scr2,buf1,buf2,  &
        score,istrt(100),ilen(100),ii(9),io(9),andf,  &
        rshift,dof(6),ip(6),aaa(2),scconn,combo,outt
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
    COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
    COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub
    COMMON /cmbfnd/ inam(2),ierr
    COMMON /BLANK / step,idry
    COMMON /zzzzzz/ z(1)
    DATA    aaa   / 4HCMCO,4HNT   /
 
    icor  = score
    iclen = lcore
    mfile = scsfil
    CALL OPEN (*200,scsfil,z(buf3),0)
    ofile = scr2
    ifile = scr1
    nwd   = 2 + npsub
    odd   = .false.
 
    DO  i = 1,npsub
        odd   = .NOT.odd
        ncsub = combo(i,5)
   
        !     READ IN EQSS FOR ITH PSEUDO-STRUCTURE
   
        mfile = ifile
        CALL OPEN (*200,ifile,z(buf1),0)
        mfile = ofile
        CALL OPEN (*200,ofile,z(buf2),1)
   
        !     MOVE TO FIRST COMPONENT EQSS
   
        DO  j = 1,ncsub
            mfile = scsfil
            CALL READ (*210,*10,scsfil,z(score),lcore,1,nnn)
            GO TO 220
10          istrt(j) = score
            ilen(j)  = nnn
            score = score + nnn
            lcore = lcore - nnn
        END DO
        CALL skpfil (scsfil,1)
   
        !     CONNECTION ENTRIES IN TERMS OF GRID POINT ID ARE ON SCR1
        !     IN THE FORM...
        !        C/CC/G1/G2/G3/G4/G5/G6/G7
   
        !     READ CONNECTION ENTRY..
   
        mfile = ifile
30      CALL READ (*110,*40,ifile,ii,10,1,nnn)
40  CONTINUE
    icomp = ii(2+i)/1000000
    igrid = ii(2+i) - 1000000*icomp
    IF (igrid == 0) GO TO 100
   
    !     THE ABOVE RETRIEVED THE ORIGINAL GRID PT. NO., NOW FIND OUT
    !     IF IT HAS SEVERAL IP NO.
   
    IF (ilen(icomp) == 0) GO TO 50
    CALL gridip (igrid,istrt(icomp),ilen(icomp),ip,dof,nip,z,nnn )
    IF (ierr /= 1) GO TO 70
50  idry = -2
    WRITE  (outt,60) ufm,igrid,combo(i,1),combo(i,2),icomp
60  FORMAT (a23,' 6535, A MANUAL CONNECTION SPECIFIES GRID ID ',i8,  &
        ' OF PSEUDOSTRUCTURE ',2A4, /30X,  &
        'COMPONENT STRUCTURE,I4,22H WHICH DOES NOT EXIST.')
    GO TO 30
    70 DO  j = 1,nip
        ii2 = rshift(dof(j),26)
        ii2 = lshift(ii2,26)
        dof(j) = dof(j) - ii2
        io(1)  = andf(ii(1),dof(j))
        IF (io(1) == 0) CYCLE
        io(2) = ii(2)
        DO  jj = 1,nwd
            io(2+jj) = ii(2+jj)
        END DO
        io(2+i) = ip(j)
        CALL WRITE (ofile,io,nwd,1)
    END DO
    GO TO 30
100 CALL WRITE (ofile,ii,nwd,1)
    GO TO 30
110 CALL CLOSE (ifile,1)
    IF (i == npsub ) CALL CLOSE (ofile,2)
    IF (i < npsub ) CALL CLOSE (ofile,1)
    isave = ifile
    ifile = ofile
    ofile = isave
END DO
scconn = scr1
IF (odd) scconn = scr2
IF (scconn == scr1) scr1 = 305
IF (scconn == scr2) scr2 = 305
score = icor
lcore = iclen
CALL CLOSE (scsfil,1)
RETURN
 
200 imsg = -1
GO TO 230
210 imsg = -2
GO TO 230
220 imsg = -8
230 CALL mesage (imsg,mfile,aaa)

RETURN
END SUBROUTINE cmcont
