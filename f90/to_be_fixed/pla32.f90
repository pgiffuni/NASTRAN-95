SUBROUTINE pla32
     
!     THIS ROUTINE READS THE ESTNLS DATA BLOCK CREATED IN SUBROUTINE
!     PLA31, AND CALLS THE PROPER ELEMENT ROUTINE TO COMPUTE ELEMENT
!     STRESSES.
!     ELEMENT STRESS INFORMATION IS STORED BY THE ELEMENT ROUTINE IN
!     /STROUT/.  THE ELEMENT ROUTINE ALSO UPDATES THE EST ENTRY WHICH
!     HAS BEEN COMMUNICATED TO IT VIA /PLA32E/.  NOTE THAT THIS UPDATED
!     EST ENTRY DOES NOT CONTAIN DISPLACEMENT VECTOR INFORMATION.
 
 INTEGER :: bufsz,bufr1,bufr2,bufr3,cstm,dit,estnls,casecc,  &
     onles,estnl1,eor,clsrw,FILE,iz(1),iestbk(100),  &
     estwds(40),eltype,outrw,placnt,plsetn,planos(2), ostrt(7),estt(7),setno
 DIMENSION       NAME(2),nstwds(40),nwdsp2(40),p(4),ip(4),dum2(2),  &
     tubsav(20),ICHAR(9),ititle(3),iy(30)
 COMMON /BLANK / placnt,plsetn
 COMMON /system/ bufsz
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /zzzzzz/ z(1)
 COMMON /pla32c/ gamma,gammas,ipass
 COMMON /pla32e/ estbk(100)
 
!     SCRATCH BLOCK USED BY ELEMENT ROUTINES (325 SINGLE PRECISION
!     CELLS)  AND OUTPUT BLOCK FOR ELEMENT STRESSES
 
 COMMON /pla32s/ xxxxxx(325)
 COMMON /sout  / yyyyyy(30)
 EQUIVALENCE     (z(1),iz(1)),(estbk(1),iestbk(1)),(p(1),ip(1)),  &
     (yyyyyy(1),iy(1))
 DATA    NAME  / 4HPLA3, 4H2          /
 DATA    ititle/ 4HLOAD, 4H fac,4HTOR /
 DATA    cstm  , mpt,dit,estnls,casecc/ 101,102,103,301,106 /
 DATA    onles , estnl1  / 201,202    /
 DATA    inrw  , outrw,eor,neor,clsrw / 0,1,1,0,1 /
 DATA    planos/ 1103,11 /
 
!    1        ROD       BEAM      TUBE      SHEAR     TWIST
!    2        TRIA1     TRBSC     TRPLT     TRMEM     CONROD
!    3        ELAS1     ELAS2     ELAS3     ELAS4     QDPLT
!    4        QDMEM     TRIA2     QUAD2     QUAD1     DAMP1
!    5        DAMP2     DAMP3     DAMP4     VISC      MASS1
!    6        MASS2     MASS3     MASS4     CONM1     CONM2
!    7        PLOTEL    REACT     QUAD3     BAR       CONE
!    8        TRIARG    TRAPRG    TORDRG    CORE      CAP
 
 DATA    estwds/ 21,        0,       20,        0,        0,  &
     38,        0,        0,       27,       21,  &
     0,        0,        0,        0,        0,  &
     32,       32,       37,       43,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,       50,        0,  &
     0,        0,        0,        0,        0 /
 DATA    nstwds/ 5,        0,        5,        0,        0,  &
     17,        0,        0,        8,        5,  &
     0,        0,        0,        0,        0,  &
     8,       17,       17,       17,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,       16,        0,  &
     0,        0,        0,        0,        0 /
 DATA    nwdsp2/ 33,        0,       32,        0,        0,  &
     56,        0,        0,       36,       33,  &
     0,        0,        0,        0,        0,  &
     44,       50,       61,       67,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,        0,        0,  &
     0,        0,        0,       62,        0,  &
     0,        0,        0,        0,        0 /
 
!     DEFINE POSITION IN CASECC RECORD OF DESTINATION (PRINTER, PUNCH,
!     ETC.) OF ELEMENT STRESSES.
 
 DATA          idest /24/
 
 
!     DETERMINE SIZE OF CORE, DEFINE BUFFERS AND INITIALIZE CORE POINTER
!     AND COUNTERS
 
 izmax = korsz(z)
 bufr1 = izmax - bufsz
 bufr2 = bufr1 - bufsz
 bufr3 = bufr2 - bufsz
 left  = bufr3 - 1
 ipass = placnt- 1
 icstm = 0
 ncstm = 0
 DO  i = 1,7
   ostrt(i) = 0
   estt(i)  = 0
 END DO
 
!     ATTEMPT TO READ CSTM INTO CORE
 
 FILE = cstm
 CALL OPEN (*20,cstm,z(bufr1),inrw)
 CALL fwdrec (*9020,cstm)
 CALL READ (*9020,*10,cstm,z(icstm+1),left,eor,ncstm)
 CALL mesage (-8,0,NAME)
 10 left = left - ncstm
 CALL CLOSE (cstm,clsrw)
 CALL pretrs (z(icstm+1),ncstm)
 20 imat = ncstm
 
!     COMPUTE GAMMA AND GAMMAS FROM THE PROPER PLFACT CARD
 
 FILE = mpt
 CALL preloc (*9010,z(bufr1-3),mpt)
 CALL locate (*9040,z(bufr1-3),planos,iflag)
 30 CALL READ (*9020,*9030,mpt,setno,1,neor,iflag)
 IF (setno == plsetn) GO TO 40
 35 CALL READ (*9020,*9030,mpt,nn,1,neor,iflag)
 IF (nn == -1) GO TO 30
 GO TO 35
 40 IF (placnt <= 4) GO TO 45
 CALL READ (*9020,*9030,mpt,0,-(placnt-4),neor,iflag)
 45 nwdsrd = 4
 IF (placnt < 4) nwdsrd = placnt
 CALL READ (*9020,*9030,mpt,p,nwdsrd,neor,iflag)
 IF (ip(nwdsrd) /= -1) GO TO 48
 IF (placnt-3 < 0.0) THEN
   GO TO    42
 ELSE IF (placnt-3 == 0.0) THEN
   GO TO    43
 ELSE
   GO TO    44
 END IF
 42  gammas = 1.0
 GO TO 46
 43  gammas = (p(2) - p(1))/p(1)
 GO TO 46
 44  gammas = (p(3) - p(2))/(p(2) - p(1))
 46  gamma = 1.0
 GO TO 65
 48  a = p(2) - p(1)
 IF (placnt-3 < 0.0) THEN
   GO TO    50
 ELSE IF (placnt-3 == 0.0) THEN
   GO TO    55
 ELSE
   GO TO    60
 END IF
 50 gammas = 0.0
 gamma  = a/p(1)
 GO TO 65
 55  gammas = a/p(1)
 gamma = (p(3) - p(2))/a
 GO TO 65
 60 word = p(3) - p(2)
 gammas = word/a
 gamma  = (p(4) - p(3))/word
 65 CALL CLOSE (mpt,clsrw)
 
!     READ MPT AND DIT FILES.  NOTE MINUS SIGN ON DIT TO TRIGGER PLA
!     FLAG.
 
 CALL premat (iz(imat+1),z(imat+1),z(bufr1-3),left,mused,mpt,-dit)
 left = left  - mused
 icc  = ncstm + mused
 
!     READ CASECC INTO OPEN CORE
 
 FILE = casecc
 CALL OPEN (*9010,casecc,z(bufr1),inrw)
 CALL fwdrec (*9020,casecc)
 CALL READ (*9020,*68,casecc,z(icc+1),left,eor,ncc)
 CALL mesage (-8,0,NAME)
 68 left = left - ncc
 CALL CLOSE (casecc,clsrw)
 
! OPEN INPUT FILE
 
 FILE = estnls
 CALL OPEN (*9010,estnls,z(bufr1),inrw)
 CALL fwdrec (*9020,estnls)
 
!     OPEN THE ELEMENT STRESS FILE FOR OUTPUT AND BUILD HEADER WHICH IS
!     NON-CHANGING.
 
 FILE = onles
 CALL OPEN (*9010,onles,z(bufr2),outrw)
 CALL fname (onles,dum2)
 CALL WRITE (onles,dum2,2,eor)
 
!     THE FOLLOWING INDICES HAVE TO CHANGE  WHEN THERE ARE CHANGES IN
!     THE FORMAT OF THE CASECC DATA BLOCK
 
 ionles = icc + ncc
 iz(ionles+1) = iz(icc+18) + 100
 iz(ionles+2) = 5
 iz(ionles+4) = iz(icc+1)
 iz(ionles+5) = iz(icc+4)
 iz(ionles+6) = 0
 iz(ionles+7) = 0
 iz(ionles+8) = 0
 iz(ionles+9) = 0
 ilow  = ionles + 51
 ihigh = ionles + 146
 left  = left - 146
 IF (left < 0) CALL mesage (-8,0,NAME)
 j = icc + 38
 DO  i = ilow,ihigh
   j = j + 1
   iz(i) = iz(j)
 END DO
 
!     STORE LOAD FACTOR AND INTEGER IN LABEL PORTION OF OUTPUT
 
 iz(ionles+135) = ititle(1)
 iz(ionles+136) = ititle(2)
 iz(ionles+137) = ititle(3)
 iii = placnt - 1
 CALL int2al (iii,iz(ionles+138),ICHAR)
 
!     DEFINE DESTINATION OF OUTPUT
 
 i = icc + idest
 jdest = iz(i)
 
!     OPEN THE ESTNL1 FILE FOR OUTPUT.
 
 FILE = estnl1
 CALL OPEN (*9010,estnl1,z(bufr3),outrw)
 CALL fname (estnl1,dum2)
 CALL WRITE (estnl1,dum2,2,eor)
 FILE = estnls
 
!     READ ELEMENT TYPE
 
 80 CALL READ (*220,*9030,estnls,eltype,1,neor,iflag)
 
!     FILL IN REMAINDER OF ID RECORD FOR THE ONLES FILE
 
 iz(ionles+3)  = eltype
 iz(ionles+10) = nstwds(eltype)
 IF (nstwds(eltype) <= 0) CALL mesage (-30,91,eltype)
 
!     WRITE ID RECORD FOR ONLES FILE
 
 CALL WRITE (onles,iz(ionles+1),146,eor)
 CALL WRITE (estnl1,eltype,1,neor)
 
!     READ AN ENTRY FROM THE APPENDED ESTNL FILE AND CALL THE PROPER
!     ROUTINE
 
 90 CALL READ (*9020,*210,estnls,estbk,nwdsp2(eltype),neor,iflag)
 
!               1,ROD    2,BEAM    3,TUBE   4,SHEAR   5,TWIST
 GO TO (     110,      999,      120,      999,      999,
!             6,TRIA1   7,TRBSC   8,TRPLT   9,TRMEM 10,CONROD  &
 130,      999,      999,      140,      110,
!            11,ELAS1  12,ELAS2  13,ELAS3  14,ELAS4  15,QDPLT  &
 999,      999,      999,      999,      999,
!            16,QDMEM  17,TRIA2  18,QUAD2  19,QUAD1  20,DAMP1  &
 150,      160,      170,      180,      999,
!            21,DAMP2  22,DAMP3  23,DAMP4   24,VISC  25,MASS1  &
 999,      999,      999,      999,      999,
!            26,MASS2  27,MASS3  28,MASS4  29,CONM1  30,CONM2  &
 999,      999,      999,      999,      999,
!           31,PLOTEL  32,REACT  33,QUAD3    34,BAR   35,CONE  &
 999,      999,      999,      190,      999,
!           36,TRIARG 37,TRAPRG 38,TORDRG   39,CORE?   40,CAP?  &
 999,      999,      999,      999,      999), eltype
 
!     ROD, CONROD
 
 110 CALL psrod
 
!     IF ELEMENT IS A TUBE, RESTORE SAVED EST ENTRY AND STORE UPDATED
!     STRESS VARIABLES IN PROPER SLOTS.
 
 IF (eltype /= 3) GO TO 200
 DO  i = 1,16
   estbk(i) = tubsav(i)
 END DO
 estbk(17) = estbk(18)
 estbk(18) = estbk(19)
 estbk(19) = estbk(20)
 estbk(20) = estbk(21)
 GO TO 200
 
 
!     TUBE - REARRANGE ESTBK FOR THE TUBE SO THAT IT IS IDENTICAL TO THE
!            ONE FOR THE ROD
 
!     SAVE THE EST ENTRY FOR THE TUBE EXCEPT THE 4 WORDS WHICH WILL BE
!     UPDATED BY THE THE ROD ROUTINE AND THE DISPLACEMENT VECTORS
 
 120 DO  i = 1,16
   tubsav(i) = estbk(i)
 END DO
 
!     COMPUTE AREA, TORSIONAL INERTIA TERM AND STRESS COEFFICIENT
 
 d  = estbk(5)
 t  = estbk(6)
 dmt= d - t
 a  = dmt*t* pi
 fj = .25*a*(dmt**2 + t**2)
 c  = d/2.0
 
!     MOVE THE END OF THE ESTBK ARRAY DOWN ONE SLOT SO THAT ENTRIES 7
!     THRU 32 WILL BE MOVED TO POSITIONS 8 THRU 33.
 
 m = 33
 DO  i = 1,26
   estbk(m) = estbk(m-1)
   m = m - 1
 END DO
 estbk(5) = a
 estbk(6) = fj
 estbk(7) = c
 GO TO 110
 
!     TRIA1
 
 130 CALL pstri1
 GO TO 200
 
!     TRMEM
 
 140 CALL pstrm
 GO TO 200
 
!     QDMEM
 
 150 CALL psqdm
 GO TO 200
 
!     TRIA2
 
 160 CALL pstri2
 GO TO 200
 
!     QUAD2
 
 170 CALL psqad2
 GO TO 200
 
!     QUAD1
 
 180 CALL psqad1
 GO TO 200
 
!     BAR
 
 190 CALL psbar
 GO TO 200
 
!     ALTER ELEMENT IDENTIFICATION FROM EXTERNAL (USER) IDENTIFICATION
!     TO INTERNAL ID., AND WRITE OUTPUT FILES.
 
 200 iy(1) = 10*iy(1) + jdest
 CALL WRITE (onles,iy,nstwds(eltype),neor)
 CALL WRITE (estnl1,estbk,estwds(eltype),neor)
 ostrt(2) = ostrt(2) + 1
 estt(2)  = estt(2)  + 1
 GO TO 90
 
!     WRITE EORS
 
 210 CALL WRITE (onles,0,0,eor)
 CALL WRITE (estnl1,0,0,eor)
 GO TO 80
 
!     CLOSE FILES AND WRITE TRAILERS
 
 220 CALL CLOSE (onles,clsrw)
 CALL CLOSE (estnl1,clsrw)
 CALL CLOSE (estnls,clsrw)
 ostrt(1) = onles
 estt(1)  = estnl1
 CALL wrttrl (ostrt)
 CALL wrttrl (estt)
 RETURN
 
!     FATAL ERRORS
 
 999 CALL mesage (-30,92,eltype)
 9010 j = -1
 GO TO 9050
 9020 j = -2
 GO TO 9050
 9030 j = -3
 GO TO 9050
 9040 j = -5
 9050 CALL mesage (j,FILE,NAME)
 RETURN
END SUBROUTINE pla32
