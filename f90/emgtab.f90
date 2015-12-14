SUBROUTINE emgtab
!*****
!     THIS ROUTINE OF THE -EMG- MODULE PREPARES OPEN CORE WITH SOME
!     VARIOUS TABLES.  CSTM, MAT, ETC.
 
!     UTILITY ROUTINES ARE USED FOR THE MOST PART.
!*****
 LOGICAL :: anycon, error, heat
 INTEGER :: rdrew, wrt, wrtrew, cls, clsrew, buf1, subr(2),  &
     precis, sysbuf, est, cstm, dit, geom2, z, FILE, eor, rd, flags, ditfil
 COMMON /system/ ksystm(65)
 COMMON /names / rd, rdrew, wrt, wrtrew, clsrew, cls
 COMMON /emgprm/ icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat ,icmbar, lcstm, lmat, lhmat
 COMMON /emgfil/ est, cstm, mpt, dit, geom2
 COMMON /hmatdd/ iihmat, nnhmat, mptfil, ditfil
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm(1), sysbuf)
 DATA    subr  / 4HEMGT,  4HAB  /,   eor/ 1 /
!*****
!     READ -CSTM- INTO CORE.
!*****
 buf1 = ncore - sysbuf - 2
 icrq = jcore - buf1
 IF (buf1 <= jcore) GO TO 10
 icstm = jcore
 ncstm = jcore - 1
 FILE  = cstm
 CALL OPEN (*30,cstm,z(buf1),rdrew)
 CALL fwdrec (*30,cstm)
 CALL READ (*60,*20,cstm,z(icstm),buf1-jcore,eor,lcstm)
 icrq = buf1 - jcore
 10 CALL mesage (-8,icrq,subr)
 20 CALL CLOSE (cstm,clsrew)
 ncstm = icstm + lcstm - 1
 CALL pretrs (z(icstm),lcstm)
 CALL pretrd (z(icstm),lcstm)
!*****
!     HAMT AND PREMAT
!*****
 30 IF (.NOT.heat) GO TO 40
 
!     HEAT PROBLEM THUS USE -HMAT-
 
 imat   = ncstm + 1
 nmat   = ncstm
 iihmat = nmat
 nnhmat = ncore
 mptfil = mpt
 ditfil = dit
 CALL prehma (z)
 ihmat = iihmat
 nhmat = nnhmat
 lhmat = nhmat - ihmat
 jcore = nhmat + 1
 GO TO 50
 
!     NON-HEAT PROBLEM THUS USE -MAT-
 
 40 imat = ncstm + 1
 CALL premat (z(imat),z(imat),z(buf1),buf1-imat,lmat,mpt,dit)
 nmat  = imat + lmat - 1
 ihmat = nmat + 1
 nhmat = nmat
 jcore = nhmat + 1
 
 50 CONTINUE
 RETURN
 
 60 CALL mesage (-2,FILE,subr)
 RETURN
END SUBROUTINE emgtab
