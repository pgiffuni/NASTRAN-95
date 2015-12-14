SUBROUTINE dlpt2 (INPUT,w1jk,w2jk)
     
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: w1jk
 INTEGER, INTENT(IN OUT)                  :: w2jk
 INTEGER :: sysbuf,ecore,tw1jk,tw2jk,NAME(2)
 DIMENSION       a(2),np(4),iz(1)
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /amgp2 / tw1jk(7),tw2jk(7)
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ work(1)
 EQUIVALENCE     (np(2),nstrip),(np(3),ntp)
 EQUIVALENCE     (work(1),iz(1))
 DATA     NAME / 4HDLPT,4H2   /
 
!     READ IN NP,NSIZE,NTP,F
 
 CALL READ(*999,*999,INPUT,np,4,0,n)
 
!     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
 
 ecore = korsz(iz)
 ecore = ecore - 4*sysbuf
 nn = ii +1
 inc = 0
 inb = inc + np(1)
 iys = inb + np(1)
 izs = iys + nstrip
 iee = izs + nstrip
 isg = iee + nstrip
 icg = isg + nstrip
 ixic = icg + nstrip
 idelx = ixic + ntp
 ixlam = idelx + ntp
 nread = ixlam + ntp
 
!     FILL IN DATA
 
 IF(nread > ecore) GO TO 998
 CALL READ(*999,*999,INPUT,work,nread,1,n)
 
!     COMPUTE TERMS AND PACK
 
 DO  i = 1,ntp
   a(1) = 0.0
   a(2) = 1.0
   CALL pack(a,w1jk,tw1jk)
   a(1) = -(2.0/refc)
   a(2)=work(idelx+i) / (2.0*refc)
   CALL pack(a,w2jk,tw2jk)
   
!     BUMP PACK INDEXES
   
   ii = ii +2
   IF(i == ntp) CYCLE
   nn = nn + 2
 END DO
 RETURN
 
!     ERROR MESSAGES
 
!     NOT ENOUGH CORE
 998 CALL mesage(-8,0,NAME)
!     FILE NOT POSITIONED PROPERLY
 999 CALL mesage(-7,0,NAME)
 RETURN
END SUBROUTINE dlpt2
