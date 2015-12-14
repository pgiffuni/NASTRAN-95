SUBROUTINE plttra
     
!     PLTTRA MODIFIES THE SIL AND BGPDT TABLES FOR THE PURPOSE OF
!     PLOTTING SPECIAL SCALAR GRID POINTS
 
!     INPUT  SIL  BGPDT  LUSET
!     OUTPUT SIP  BGPDP  LUSEP
 
!     SPECIAL SCALAR GRID POINTS
!     BGPDT(I,1)= 1  SIL(I+1)-SIL(I)=1
!     BGPDP(I,1)=-2  SIP(I+1)-SIP(I)=6
 
!     LUSET IS THE VALUE OF SIL(LAST+1) IF IT EXISTED
!     LUSEP IS THE VALUE OF SIP(LAST+1) IF IT EXISTED
 
 LOGICAL :: leof
 INTEGER :: sysbuf,buf1,buf2,buf3,buf4,FILE,sil,bgpdt,sip,  &
     bgpdp,NAME(2),z,plt(2),flag,a,b,s1,s2,delta
 DIMENSION       a(4),b(2),mcb(7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / luset,lusep
 COMMON /system/ sysbuf,NOT
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (a(3),b(1)),(FILE,mcb(1))
 DATA    bgpdt , sil,bgpdp,sip/ 101,102,201,202 /
 DATA    plt   / 4HPLTT,4HRA  /,   mcb / 7*0    /
 DATA    leof  / .false./
 
 nadd = 0
 ns   = 0
 
!     LOCATE STORAGE AREA FOR FILE BUFFERS
 
 nz   = korsz(z)
 buf1 = nz   - sysbuf + 1
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 IF (buf4 <= 0) CALL mesage (-8,nz,plt)
 
!     READ TRAILER RECORDS OF INPUT FILES AND CHECK COMPATABILITY
!     OPEN AND FOREWARD SPACE LABEL RECORD OF INPUT FILES
!     OPEN AND WRITE LABEL RECORD OF OUTPUT FILES
 
 FILE = bgpdt
 CALL rdtrl (mcb)
 CALL fname (FILE,NAME)
 IF (FILE <= 0) GO TO 900
 CALL OPEN (*900, bgpdt, z(buf2), rdrew)
 CALL fwdrec (*1010,bgpdt)
 
 FILE = sil
 CALL rdtrl (mcb)
 CALL fname (FILE,NAME)
 IF (FILE <= 0) GO TO 900
 IF (mcb(3) /= luset) GO TO 1130
 CALL OPEN (*900,sil,z(buf1),rdrew)
 CALL fwdrec (*1010,sil)
 
 FILE = sip
 CALL fname (sip,a)
 CALL OPEN  (*1000,sip,z(buf3),wrtrew)
 CALL WRITE (sip,a,2,1)
 
 FILE = bgpdp
 CALL OPEN  (*1000,bgpdp,z(buf4),wrtrew)
 CALL fname (bgpdp,b)
 CALL WRITE (bgpdp,b,2,1)
 
!     READ SIL(I)
 
 FILE = sil
 CALL READ (*1010,*1020,sil,s1,1,0,flag)
 
!     READ SIL(I+1)
 
 10 FILE = sil
 CALL READ (*1010,*30,sil,s2,1,0,flag)
 
!     READ BGPDT(I,J)
 
 15 FILE = bgpdt
 CALL READ (*1010,*1020,bgpdt,a,4,0,flag)
 delta = 0
 ns = ns + 1
 
!     CHECK IF SPECIAL SCALAR GRID POINT
 
 IF (a(1) < 0 .OR. s2-s1 == 6) GO TO 20
 IF (s2-s1 /= 1) GO TO 1110
 
!     SPECIAL SCALAR GRID POINT
 
 delta = 5
 a(1)  =-2
 20 s1 = s1 + nadd
 
!     WRITE SIP AND BGPDP TABLE ENTRIES
 
 CALL WRITE (sip,s1,1,0)
 CALL WRITE (bgpdp,a,4,0)
 nadd = nadd + delta
 IF (leof) GO TO 40
 s1 = s2
 GO TO 10
 
!     SIL(I) IS SIL(LAST)
 
 30 leof = .true.
 s2   = luset + 1
 GO TO 15
 40 lusep = luset + nadd
 
!     CLOSE OUTPUT FILES AND WRITE TRAILER RECORDS
 
 CALL CLOSE (sil  ,clsrew)
 CALL CLOSE (bgpdt,clsrew)
 CALL CLOSE (sip  ,clsrew)
 CALL CLOSE (bgpdp,clsrew)
 mcb(1) = bgpdp
 mcb(3) = 0
 CALL wrttrl (mcb)
 mcb(1) = sip
 mcb(3) = lusep
 CALL wrttrl (mcb)
 RETURN
 
 900 lusep = luset
 RETURN
 
!     ERROR DIAGNOSTICS
 
 1000 ndx = -1
 GO TO 1100
 1010 ndx = -2
 GO TO 1100
 1020 ndx = -3
 1100 CALL mesage (ndx,FILE,plt)
 GO TO 1150
 1130 WRITE  (NOT,2001) ufm,luset,mcb(3)
 2001 FORMAT (a23,' 5011, FIRST PARAMETER',i6,' NE TRAILER RECORD ',  &
     'PARAMETER',i6)
 GO TO 1150
 1110 WRITE  (NOT,2002) ufm,ns
 2002 FORMAT (a23,' 5012, ENTRY',i6,' OF SIL TABLE INCOMPATIBLE WITH ',  &
     'NEXT ENTRY')
 1150 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE plttra
