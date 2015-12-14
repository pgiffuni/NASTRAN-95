SUBROUTINE amgsba(ajjl,a0,ar,nsbe,a,yb,zb)
     
    !     BUILD AJJL FOR DOUBLET LATTICE WITH BODIES
 
 
    INTEGER, INTENT(IN OUT)                  :: ajjl
    REAL, INTENT(IN)                         :: a0(1)
    REAL, INTENT(IN OUT)                     :: ar(1)
    INTEGER, INTENT(IN)                      :: nsbe(1)
    REAL, INTENT(OUT)                        :: a(1)
    REAL, INTENT(IN OUT)                     :: yb(1)
    REAL, INTENT(IN OUT)                     :: zb(1)
    INTEGER :: sysbuf,ecore, NAME(2),scr1,scr2,scr5
 
 
    COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
        inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
        iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
        icg,ixij,ix,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria  &
        ,inasb,ifla1,ifla2,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5
    COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
    COMMON /system/ sysbuf
    COMMON /zzzzzz / z(1)
    COMMON /packx/ iti,ito,ii,nn,incr
    COMMON /condas/ pi,twopi
    DATA NAME /4HAMGS,4HBA  /
 
    ii = nrow +1
    nn = nrow + nj1
    IF(next+2*nj1 > ecore) CALL mesage(-8,0,NAME)
    incr = 1
    iti = 3
    ito = 3
    IF(nt0 == 0) GO TO 100
    nbuf = 1
    IF(ntzs /= 0) nbuf = nbuf +1
    IF(ntys /= 0) nbuf = nbuf +1
    ibuf1 = ecore - nbuf*sysbuf
    ibuf2 = ibuf1
    IF(ntzs /= 0) ibuf2 = ibuf1+sysbuf
    ibuf3 = ibuf2
    IF(ntys /= 0) ibuf3 = ibuf2+sysbuf
    ns5 = nt0*2
    ns1 = ntzs*2
    ns2 = ntys*2
    ntot = ns5+ns1+ns2
    IF(next+ntot > ibuf3) CALL mesage(-8,0,NAME)
 
    !     BUILD PANEL AND BODY PART OF AJJL
 
    CALL gopen(scr5,z(ibuf1),0)
    IF(ntzs /= 0) CALL gopen(scr1,z(ibuf2),0)
    IF(ntys /= 0) CALL gopen(scr2,z(ibuf3),0)
    DO  i=1,nt0
        CALL fread(scr5,a,ns5,0)
        IF(ntzs /= 0) CALL fread(scr1,a(ns5+1),ns1,0)
        IF(ntys /= 0) CALL fread(scr2,a(ns5+ns1+1),ns2,0)
        CALL pack(a,ajjl,mcb)
    END DO
    CALL CLOSE(scr5,1)
    CALL CLOSE(scr1,1)
    CALL CLOSE(scr2,1)
100 CALL zeroc(a,2*nj1)
 
    !     ADD DIAGIONAL TERMS OF AJJL FOR SLENDER BODIES
 
    IF(ntzs == 0.AND.ntys == 0) GO TO 1000
    i = nt0*2+1
    den=twopi*2.0
    IF(ntzs == 0) GO TO 200
    nfsbeb = 1
    nlsbeb = 0
    DO  ib = 1,nbz
        nlsbeb = nlsbeb + nsbe(ib)
        DO  it = nfsbeb,nlsbeb
            a(i) = 1.0 / (den*a0(it)**2)
            IF(ABS(yb(ib)) < .00001) a(i) = (1.0+FLOAT(nd))*a(i)
            IF(ABS(zb(ib)) < .00001) a(i) = (1.0+FLOAT(NE))*a(i)
            CALL pack(a,ajjl,mcb)
            a(i) = 0.0
            i = i+2
        END DO
        nfsbeb = nfsbeb + nsbe(ib)
    END DO
200 IF(ntys == 0) GO TO 1000
    nfyb = nb+1-nby
    nfsbeb = 1
    nlsbeb = 0
    nl = nfyb-1
    IF(nl == 0) GO TO 220
    DO  j=1,nl
        nlsbeb = nlsbeb+nsbe(j)
        nfsbeb = nfsbeb+nsbe(j)
    END DO
    220 DO  ib = nfyb,nb
        nlsbeb = nlsbeb+nsbe(ib)
        DO  it = nfsbeb,nlsbeb
            a(i) = 1.0 / (den*a0(it)**2)
            IF(ABS(yb(ib)) < .00001) a(i) = (1.0-FLOAT(nd))*a(i)
            IF(ABS(zb(ib)) < .00001) a(i) = (1.0-FLOAT(NE))*a(i)
            CALL pack(a,ajjl,mcb)
            a(i) = 0.0
            i = i+2
        END DO
        nfsbeb = nfsbeb+nsbe(ib)
    END DO
1000 RETURN
END SUBROUTINE amgsba
