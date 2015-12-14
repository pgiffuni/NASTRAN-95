SUBROUTINE cmtoc
     
    !     THIS SUBROUTINE GENERATES A TABLE OF CONTENTS FOR A COMBINE
    !     OPERATION. FOR EACH PSEUDO-STRUCTURE IT LISTS THE NAME, NUMBER
    !     OF COMPONENTS, AND EACH COMPONENT BASIC SUBSTRUCTURE NAME.
    !     THIS DATA IS THEN WRITTEN ON SCRATCH FILE SCTOC.
 
    EXTERNAL        rshift,andf
    LOGICAL :: PRINT,tocopn
    INTEGER :: sctoc,buf5,combo,NAME(2),z,score,aaa(2),outt,  &
        ihed(96),xxx,andf,rshift
    COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
    COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
    COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
        mcon,restct(7,7),isort,origin(7,3),iprint,tocopn
    COMMON /zzzzzz/ z(1)
    COMMON /output/ ititl(96),ihdr(96)
    COMMON /system/ xxx
    DATA    ihed  / 7*4H     ,  &
        4HP s , 4HE u , 4HD o , 4HS t , 4HR u , 4HC t , 4HU r ,  &
        4HE   , 4HT a , 4HB l , 4HE   , 4HO f , 4H  c , 4HO n ,  &
        4HT e , 4HN t , 4HS   , 15*4H         ,  &
        4H pse, 4HUDO-, 4H    , 4H   n, 4HO. o, 4HF   ,26*2H  ,  &
        4HSTRU, 4HCTUR, 4HE   , 4H com, 4HPONE, 4HNTS , 4H   -,  &
        4H----, 4H----, 4H- co, 4HMPON, 4HENT , 4HNAME, 4HS --,  &
        4H----, 4H----, 4H-   , 8*4H     /
    DATA    aaa   / 4HCMTO, 4HC   /
    DATA    nheqss/ 4HEQSS/
 
    PRINT = .false.
    IF (andf(rshift(iprint,1),1) == 1) PRINT = .true.
    tocopn = .true.
    itot = 0
    DO  i = 1,96
        ihdr(i) = ihed(i)
    END DO
    IF (PRINT) CALL page
    CALL OPEN (*60,sctoc,z(buf5),1)
    DO  i = 1,npsub
        NAME(1) = combo(i,1)
        NAME(2) = combo(i,2)
        CALL sfetch (NAME,nheqss,1,itest)
        CALL suread (z(score),-1,nwds,itest)
        z(score  ) = NAME(1)
        z(score+1) = NAME(2)
        CALL WRITE (sctoc,z(score),3,0)
        itot = itot + 3
        ia   = score
        ib   = score+2
        IF (PRINT) WRITE(outt,30) (z(kdh),kdh=ia,ib)
30      FORMAT (34X,2A4,6X,i4)
        combo(i,5) = z(score+2)
        nwds = nwds - 4
        ia   = score+4
        ib   = ia+nwds-1
        nt   = (ib - ia + 1)/8
        IF (nt == 0) nt = 1
        IF (PRINT) CALL page2 (nt)
        IF (PRINT) WRITE (outt,40) (z(kdh),kdh=ia,ib)
        itot = itot + nwds
40      FORMAT (1H+,57X,2X,2A4,2X,2A4,2X,2A4,2X,2A4,/  &
            (58X,2X,2A4,2X,2A4,2X,2A4,2X,2A4))
        CALL WRITE (sctoc,z(score+4),nwds,1)
    END DO
    CALL CLOSE (sctoc,1)
    CALL OPEN (*60,sctoc,z(buf5),0)
 
    !     DETERMINE WHETHER TO CLOSE FILE
 
    IF (itot <= xxx) RETURN
    tocopn = .false.
    CALL CLOSE (sctoc,1)
    RETURN
 
60  CALL mesage (-1,sctoc,aaa)

    RETURN
END SUBROUTINE cmtoc
