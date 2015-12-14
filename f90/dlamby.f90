SUBROUTINE dlamby(INPUT,matout,skj)
     
!     DRIVER FOR DOUBLET LATTICE WITH BODIES
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: matout
 INTEGER, INTENT(IN OUT)                  :: skj
 INTEGER :: ecore,sysbuf,iz(1),tskj
 INTEGER :: scr1,scr2,scr3,scr4,scr5
 DIMENSION       NAME(2)
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / nk,nj
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /dlbdy / nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,  &
     izb,iyb,iavr,iarb,infl,ixle,ixte,int121,int122,  &
     izs,iys,ics,iee,isg,icg,ixij,ix,idelx,ixic,ixlam,  &
     ia0,ixis1,ixis2,ia0p,iria,inasb,ifla1,ifla2,  &
     ith1a,ith2a,ecore,next,scr1,scr2,scr3,scr4,scr5, ntbe
 EQUIVALENCE     (iz(1),z(1))
 DATA   NAME /4HDLAM,4HBY  /
 DATA nhaero,nhpoin,nhcore / 4HAERO,4HPOIN,4HCORE/
 
 scr1 = 301
 scr2 = 302
 scr3 = 303
 scr4 = 304
 scr5 = 305
 
!     GET CORE THEN SET POINTERS TO ACPT TABLE ARRAYS
 
 ecore = korsz(iz) - 4*sysbuf
 
!     READ LENGTHS OF ARRAYS
 
 CALL fread(INPUT,nj1,13,0)
 
!     COMPUTE POINTERS TO OPEN CORE
 
 lns = inc
 inc = 1
 ins = inc
 inb = ins + np
 inas = inb + np
 izin = inas
 iyin = izin
 inbea1 = iyin + np
 inbea2 = inbea1 + nb
 insbea = inbea2 + nb
 izb = insbea + nb
 iyb = izb + nb
 iavr = iyb + nb
 iarb = iavr + nb
 infl = iarb + nb
 ixle = infl + nb
 ixte = ixle + nb
 int121 = ixte + nb
 int122 = int121 + nb
 izs = int122 + nb
 n = 3*np + 12 * nb
 
!     READ FIXED ARRAYS
 
 IF(n > ecore) GO TO 998
 CALL fread(INPUT,iz,n,0)
 
!     GET LENGTHS OF VARIABLE ARRAYS, PANELS THEN BODIES
 
 lnas = 0
 IF(np == 0) GO TO 20
 DO  i=1,np
   lnas = lnas + iz(inas+i-1)
 END DO
 20 lnb = 0
 lnsb = 0
 lnfl = 0
 lt1 = 0
 lt2 = 0
 DO  i=1,nb
   k = i-1
   lnb = lnb + iz(inbea1+k)
   lnsb = lnsb + iz(insbea+k)
   lnfl = lnfl + iz(infl+k)
   lt1 = lt1 + iz(int121+k)
   lt2 = lt2 + iz(int122+k)
 END DO
 ntbe = ntp+lnb
 
!     READ VARIABLE  ARRAYS AND SET POINTERS TO CORE
 
 next = n+1
 n = 2*nb + 5*lns + 4*ntp + 3*lnb + 4*lnsb + lnas + 2*lnfl + lt1 + lt2
 IF(next+n+4*nj >= ecore) GO TO 998
 CALL fread(INPUT,iz(next),n,1)
 next = next + n + 1
 iys = izs + nb + lns
 ics = iys
 iee = ics + nb + lns
 isg = iee + lns
 icg = isg + lns
 ixij = icg
 ix = ixij + lns
 idelx = ix + ntp + lnb
 ixic = idelx + ntp + lnb
 ixlam = ixic + ntp
 ia0 = ixlam + ntp
 ixis1 = ia0 + lnsb
 ixis2 = ixis1 + lnsb
 ia0p = ixis2 + lnsb
 iria = ia0p + lnsb
 inasb = iria + lnb
 ifla1 = inasb + lnas
 ifla2 = ifla1 + lnfl
 ith1a = ifla2 + lnfl
 ith2a = ith1a + lt1
 
!     BUILD A MATRIX
 
 CALL bug(nhaero,100,nd,5)
 CALL bug(nhpoin,100,nj1,59)
 CALL bug(nhcore,100,z,next)
 n1 = next
 n = next + 2*ntbe
 next = next + 4*ntbe
 IF(nt0 /= 0) CALL gendsb(z(inc),z(inb),z(isg),z(icg),z(infl),  &
     z(inbea1),z(inbea2),z(ifla1),z(ifla2),z(n1),z(n1),z(n))
 n = ntzs + ntys
 next = n1
 beta = SQRT(1.0-fmach**2)
 IF( nt0 /= 0 .AND. n /= 0) CALL amgrod(z(n1),beta)
 CALL amgsba(matout,z(ia0),z(iarb),z(insbea),z(n1),z(iyb),z(izb))
 nrow = nrow + nj1
 
!     BUILD SKJ MATRIX BE SURE TO BUMP ISK NSK
 
 CALL amgbfs(skj,z(iee),z(idelx),z(inc),z(inb),z(ixis2),z(ixis1),  &
     z(ia0),z(ia0p),z(insbea))
 1000 RETURN
 
!     ERROR MESSAGES
 
 998 CALL mesage(-8,0,NAME)
 GO TO 1000
END SUBROUTINE dlamby
