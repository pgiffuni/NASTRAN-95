SUBROUTINE pull (bcd,out,icol,nchar,flag)
     
!     THIS ROUTINE EXTRACTS BCD DATA (OUT) FROM A STRING,(BCD)
!     STARTING AT POSITION ICOL
 
 
 INTEGER, INTENT(IN OUT)                  :: bcd(1)
 INTEGER, INTENT(OUT)                     :: out(1)
 INTEGER, INTENT(IN OUT)                  :: icol
 INTEGER, INTENT(IN)                      :: nchar
 INTEGER, INTENT(IN OUT)                  :: flag
 EXTERNAL        orf
 LOGICAL :: first
 INTEGER :: cperwd,BLANK,orf
 COMMON /system/ idum(38),nbpc,nbpw,ncpw
 DATA    cperwd/ 4 /, BLANK / 4H     /, first / .true. /
 
 nwds  = (nchar-1)/cperwd + 1
 IF (.NOT.first) GO TO 5
 first = .false.
 nx    = ncpw - cperwd
 ixtra = nx*nbpc
 ibl   = 0
 ib1   = krshft(BLANK,ncpw-1)
 IF (nx == 0) GO TO 5
 DO  i = 1,nx
   ibl = orf(ibl,klshft(ib1,i-1))
 END DO
 5 DO  i = 1,nwds
   out(i) = ibl
 END DO
 
 iwd = (icol-1)/cperwd + 1
 m1  = (icol-(iwd-1)*cperwd-1)*nbpc
 m2  = cperwd*nbpc - m1
 
 DO  i = 1,nwds
   ibcd   = iwd + i - 1
   itemp  = krshft(bcd(ibcd),ixtra/nbpc)
   out(i) = orf(out(i),klshft(itemp,(m1+ixtra)/nbpc))
   itemp  = krshft(bcd(ibcd+1),(m2+ixtra)/nbpc)
   out(i) = orf(out(i),klshft(itemp,ixtra/nbpc))
 END DO
 IF (nwds*cperwd == nchar) RETURN
 
!     REMOVE EXTRA CHARACTERS FROM LAST OUT WORD
 
 nbl   = (nwds-1)*cperwd + ncpw - nchar
 lword = krshft(out(nwds),nbl)
 out(nwds) = klshft(lword,nbl)
 
 itemp = 0
 DO  i = 1,nbl
   itemp = orf(itemp,klshft(ib1,i-1))
 END DO
 out(nwds) = orf(out(nwds),itemp)
 
 RETURN
END SUBROUTINE pull
