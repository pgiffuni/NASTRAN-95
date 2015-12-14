SUBROUTINE amgbfs (skj,ee,delx,nc,nba,xis2,xis1,a0,a0p,nsbe)
     
    !     BUILD SKJ CALL BFSMAT THEN SHUFFEL AND DEAL
 
 
    INTEGER, INTENT(IN OUT)                  :: skj
    REAL, INTENT(IN)                         :: ee(1)
    REAL, INTENT(IN)                         :: delx(1)
    INTEGER, INTENT(IN)                      :: nc(1)
    INTEGER, INTENT(IN)                      :: nba(1)
    REAL, INTENT(IN)                         :: xis2(1)
    REAL, INTENT(IN)                         :: xis1(1)
    REAL, INTENT(IN)                         :: a0(1)
    REAL, INTENT(IN)                         :: a0p(1)
    INTEGER, INTENT(IN)                      :: nsbe(1)
    INTEGER :: NAME(2),tskj,sysbuf,scr1,scr2,ecore
 
    COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
    COMMON /dlbdy / nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
        inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,  &
        izb,iyb,iavr,iarb,infl,ixle,ixte,int121,int122,  &
        izs,iys,ics,iee,isg,icg,ixij,ix,idelx,ixic,ixlam,  &
        ia0,ixis1,ixis2,ia0p,iria,inasb,ifla1,ifla2,ith1a,  &
        ith2a,ecore,next,scr1,scr2,scr3,scr4,scr5
    COMMON /system/ sysbuf
    COMMON /zblpkx/ a(4),iis
    COMMON /zzzzzz/ z(1)
    DATA    NAME  / 4HAMGB,4HFS  /
    DATA    nhbfs , nhg,nha      / 4HBFS ,4HG   ,4HA   /
 
    nsb   = ntys + ntzs
    nzy2  = nsb*2
    nt02  = nt0*2
    ntp2  = ntp*2
    length= nt0 + nsb
    isl   = isk - 1
    ii    = isk + length
    nn    = nsk + nzy2 + ntp2
    ibuf2 = ecore
    IF (nsb == 0) GO TO 40
    ibuf2 = ecore - sysbuf
 
    !     CALL BFSMAT
    !     SCR1 HAS NTZS + NTYS ROWS WITH NTO*2 THEN NTZS+NTYS*2 TERMS
    !     ROWS ARE Z FOR Z , Y THEN Z FOR ZY , AND Y FOR Y
 
    CALL gopen (scr1,z(ibuf2),1)
    icorr = next
    IF (next+length*4 > ibuf2) GO TO 998
    CALL bfsmat (nd,NE,nb,np,ntp,length,nt0,scr1,jf,jl,z(inas),fmach,  &
        z(iyb),z(izb),z(iys),z(izs),z(ix),delx,ee,z(ixic),  &
        z(isg),z(icg),z(iarb),z(iria),z(inbea1),z(inbea2),  &
        z(inasb),z(inb),nc,z(icorr),z(iavr),refc,a0,xis1, xis2,rfk,nsbe,nt0)
    CALL WRITE  (scr1,0,0,1)
    CALL CLOSE  (scr1,1)
    CALL dmpfil (scr1,z(next),ibuf2-next)
    CALL gopen  (scr1,z(ibuf2),0)
    ncore = nt0*nzy2*2
    IF (ncore+next > ibuf2) GO TO 998
    CALL zeroc  (z(next),ncore)
    i    = next
    izbf = 1
    DO  j = 1,nsb
        CALL fread (scr1,z(i),nt02,0)
        CALL fread (scr1,z(i),-nzy2,0)
        i = i + nt02
        IF (jf == 0) GO TO 20
        IF (j < jf .OR. j > jl) GO TO 20
        izbf = -izbf
        IF (izbf < 0) CYCLE
        i = i + nt02
20      i = i + nt02
    END DO
    CALL bckrec (scr1)
 
    !     BUILD NT0 COLUMNS OF SKJ
 
40  IF (nt0 == 0) GO TO 100
    ibf  = next - 2
    k    = 1
    ks   = 1
    nbxr = nc(k)
    DO  i = 1,nt0
        CALL bldpk (3,3,skj,0,0)
        IF (i > ntp) GO TO 45
        a(1) = 2.0*ee(ks)*delx(i)
        a(2) = 0.0
        iis  = isl + (i-1)*2 + 1
        CALL zblpki
        a(1) = (ee(ks)*delx(i)**2)/2.0
        iis  = iis + 1
        CALL zblpki
        IF (i ==    ntp) GO TO 45
        IF (i == nba(k)) k = k + 1
        IF (i ==   nbxr) GO TO 44
        GO TO 45
44      ks   = ks + 1
        nbxr = nbxr + nc(k)
45      IF (nsb == 0) GO TO 60
        ibf  = ibf + 2
        DO  j = 1,nzy2
            l    = (j-1)*nt02
            a(1) = z(ibf+l  )
            a(2) = z(ibf+l+1)
            iis  = isl + ntp2 + j
            CALL zblpki
        END DO
60      CALL bldpkn (skj,0,tskj)
    END DO
 
    !     SLENDER BODY ONLY PART OF SKJ  BFS * G
 
100 IF (nsb == 0) GO TO 900
    ncore = nzy2*nsb*4 + nsb*nsb*2
    IF (ncore+next > ibuf2) GO TO 998
    CALL zeroc (z(next),ncore)
    i    = next
    izbf = 1
    DO  j = 1,nsb
        CALL fread (scr1,z(i),-nt02,0)
        CALL fread (scr1,z(i), nzy2,0)
        i = i + nzy2
        IF (jf == 0) GO TO 120
        IF (j < jf .OR. j > jl) GO TO 120
        izbf = -izbf
        IF (izbf < 0) CYCLE
        i = i + nzy2
120     i = i + nzy2
    END DO
 
    !     BFS AT NEXT  G AT IG
 
    ig   = i
    ia   = ig + nsb*nsb*2
    nfyb = nb + 1 - nby
    irow = ig
    rfkoc= 2.0*rfk/refc
    ibzy = 0
    p5   = .5
    IF (ntzs == 0) GO TO 170
    nfse = 1
    nlse = 0
    nfb  = 1
    nbx  = nbz
    141 DO  ib = nfb,nbx
        nlse = nlse + nsbe(ib)
        DO  it = nfse,nlse
            dx   = xis2(it) - xis1(it)
            a02p = 2.0/a0(it)*a0p(it)
            IF (nfse == nlse) GO TO 148
            IF (it   /= nfse) GO TO 142
            x2 = p5*(xis2(it+1) + xis1(it+1))
            x1 = p5*(xis2(it  ) + xis1(it  ))
            z(irow  ) = (-1.0/(x2-x1))*dx
            z(irow+2) = -z(irow)*(a0(it)/a0(it+1))**2
            GO TO 148
142         IF (it == nlse) GO TO 145
            x1 = p5*(xis2(it-1) + xis1(it-1))
            x2 = p5*(xis2(it  ) + xis1(it  ))
            x3 = p5*(xis2(it+1) + xis1(it+1))
            z(irow-2) = (1.0/(x3-x1) - 1.0/(x2-x1))*dx*(a0(it)/a0(it-1))**2
            z(irow)   = (1.0/(x2-x1) - 1.0/(x3-x2))*dx
            z(irow+2) = (1.0/(x3-x2) - 1.0/(x3-x1))*dx*(a0(it)/a0(it+1))**2
            GO TO 148
145         x1 = p5*(xis2(it-1) + xis1(it-1))
            x2 = p5*(xis2(it  ) + xis1(it  ))
            z(irow  ) = (1.0/(x2-x1))*dx
            z(irow-2) =-z(irow)*(a0(it)/a0(it-1))**2
148         z(irow  ) = z(irow) + dx*a02p
            z(irow+1) = dx*rfkoc
            irow = irow + nzy2 + 2
        END DO
        nfse = nfse + nsbe(ib)
    END DO
170 IF (ibzy == 1) GO TO 200
    ibzy = 1
    IF (ntys == 0) GO TO 200
    nfb  = nfyb
    nbx  = nb
    nfse = 1
    nlse = 0
    nl   = nfyb - 1
    IF (nl == 0) GO TO 141
    DO  j = 1,nl
        nlse = nlse + nsbe(j)
        nfse = nfse + nsbe(j)
    END DO
    GO TO 141
 
    !     MULTIPLY BFS * G
 
200 CALL bug (nhbfs ,200,z(next),nzy2*nzy2)
    CALL bug (nhg   ,200,z(ig),nsb*nsb*2)
    CALL gmmatc (z(next),nzy2,nsb,0,z(ig),nsb,nsb,0,z(ia))
    CALL bug (nha   ,200,z(ia),nzy2*nsb*2)
    irow = ia - 2
    DO  i = 1,nsb
        CALL bldpk (3,3,skj,0,0)
        irow = irow + 2
        k    = irow
        DO  j = 1,nzy2
            a(1) = z(k  )
            a(2) = z(k+1)
            iis  = isl + ntp2 + j
            CALL zblpki
            k    = k + nzy2
        END DO
        CALL bldpkn (skj,0,tskj)
    END DO
900 isk  = ii
    nsk  = nn
    CALL CLOSE (scr1,1)
1000 RETURN
 
    !     ERROR MESSAGES
 
998 CALL mesage (-8,0,NAME)
    GO TO 1000

END SUBROUTINE amgbfs
