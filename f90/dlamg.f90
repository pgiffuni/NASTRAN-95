SUBROUTINE dlamg(INPUT,matout,skj)
     
!     DRIVER FOR THE DOUBLET LATTICE METHOD
!     COMPUTATIONS ARE FOR THE AJJL MATRIX
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: matout
 INTEGER, INTENT(IN OUT)                  :: skj
 INTEGER :: ecore,sysbuf,iz(1)
 DIMENSION       NAME(2),a(2)
 COMPLEX :: dt(1)
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /dlcom / np,nstrip,ntp,f,njj,next,length,  &
     inc,inb,iys,izs,iee,isg,icg, ixic,idelx,ixlam,idt,ecore
 COMMON /zzzzzz/ work(1)
 COMMON /system/ sysbuf
 COMMON /BLANK / nk,nj
 EQUIVALENCE     (work(1),iz(1),dt(1))
 DATA    NAME  / 4HDLAM,4HG   /
 
 njj = nj
 
!     READ IN NP,NSIZE,NTP,F
 
 CALL READ(*999,*999,INPUT,np,4,0,n)
 
!     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
 
 ecore = korsz(iz)
 ecore =  ecore - 4 * sysbuf
 inc = 1
 inb = inc + np
 iys = inb + np
 izs = iys + nstrip
 iee = izs + nstrip
 isg = iee + nstrip
 icg = isg + nstrip
 ixic = icg + nstrip
 idelx = ixic + ntp
 ixlam = idelx + ntp
 nread = ixlam + ntp
!     IDT IS A COMPLEX POINTER
!     THE MATRIX PACKED OUT IS NJ LONG STARTING AT DT
 idt = (nread +2) / 2
 next = idt*2 + 2*nj + 1
 
!     FILL IN DATA
 
 IF(next > ecore) GO TO 998
 nread = nread -1
 CALL READ(*999,*999,INPUT,work,nread,1,n)
 
!     CHECK FOR ENOUGH SCRATCH STORAGE
 
 n = inc + np -1
 length = 1
 
!     PUT OUT SKJ
 
 iti = 1
 ito = 3
 ii = isk
 nsk = nsk + 2
 nn = nsk
 k = 0
 ks = 0
 nbxr=iz(inc+k)
 DO  i=1,ntp
   a(1) = 2.0 * work(iee+ks) * work(idelx+i-1)
   a(2) = (work(iee+ks) * work(idelx+i-1)**2)/ 2.0
   CALL pack( a,skj,tskj)
   ii = ii +2
   IF(i == ntp) CYCLE
   nn = nn +2
   IF(i == iz(inb+k)) k = k+1
   IF(i == nbxr) GO TO 4
   CYCLE
   4 ks = ks +1
   nbxr = nbxr + iz(inc+k)
 END DO
 isk = ii
 nsk = nn
 iti = 3
 ito = 3
 ii = 1
 nn = nj
 CALL gend(work(inc),work(inb),work(iys),work(izs),  &
     work(isg),work(icg),dt(idt),work(1),matout)
 nrow = nrow + ntp
 RETURN
 
!     ERROR MESSAGES
 
!     NOT ENOUGH CORE
 998 CALL mesage(-8,0,NAME)
!     INPUT NOT POSITIONED PROPERLY OR INCORRECTLY WRITTEN
 999 CALL mesage(-7,0,NAME)
 RETURN
END SUBROUTINE dlamg
