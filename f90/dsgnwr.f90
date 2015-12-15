SUBROUTINE dsgnwr
     COMMON /system/ isybuf, iwr
 INCLUDE 'XNSTRN.COM'
 INCLUDE 'GINOX.COM'
 INCLUDE 'DSIOF.COM'
 CHARACTER (LEN=4) :: cbuff(3)
 EQUIVALENCE     (cbuff,ibase)
 idsn   = mdsfcb( 2,ifilex )
 idsnr  = idsn
 10   istrb  = fcb( 5,idsnr  )
 IF ( nblock >= istrb ) GO TO 20
 inext  = mdsfcb( 3,idsnr  ) / mulq2
 GO TO 40
 20   ialloc = fcb( 7, idsnr )
 IF ( nblock <= ( ialloc+istrb-1 ) ) GO TO 50
 IF ( idsn == 8 ) CALL dsmsg( 9 )
 inext  = IAND( mdsfcb( 3,idsnr ), maskh2 )
 IF ( inext /= 0 ) GO TO 40
 maxpr1 = maxpri + 1
 DO  i = maxpr1, maxfcb
   iavail = mdsfcb( 3,i )
   IF ( iavail /= 0 ) CYCLE
   ifirst = ialloc + istrb
   ialloc = 20000000
   fcb( 5,i ) = ifirst
   fcb( 6,i ) = ifirst-1
   mdsfcb( 3,i ) = idsnr * mulq2
   inext  = i
   mdsfcb( 3,idsnr ) = ior( mdsfcb( 3,idsnr ), i )
   EXIT
 END DO
 40   idsnr  = inext
 IF ( idsnr >= 1 .AND. idsnr <= maxdsn ) GO TO 10
 CALL dsmsg ( 122 )
 50   IF ( idsn == idsnr ) GO TO 60
 CALL dsclos( idsn )
 mdsfcb( 1,idsn ) = IAND( mdsfcb( 1,idsn ), maskh1 )
 idsn   = idsnr
 mdsfcb( 1,idsn )   =  ior( mdsfcb( 1,idsn ), maskh2 )
 mdsfcb( 2,ifilex ) = idsn
 CALL dsmsg( 8 )
 idevic = 0
 DO  kk = 1, numdev
   mdsnam(idsn)(1:2) = dev(kk)
   isave = iop
   iop = 0
   CALL dsopen( mdsnam( idsn ), idsn, iop )
   iop = isave
   cbuff( indbas ) = mdsnam( idsn )
   CALL dswrit( idsn, ibase( indbas+3 ), nbuff, ioblk, iccer )
   IF ( iccer == 0 ) GO TO 60
   CALL dsclos (idsn)
 END DO
 57   WRITE ( iwr, 901 )
 901   FORMAT(///,' NO MORE DISK SPACE AVAILABLE, JOB ABORTED.')
 CALL dsmsg( 122 )
 60   ioblk  = nblock - istrb + 1
 CALL dswrit( idsn, ibase( indbas+3 ), nbuff, ioblk, iccer )
 IF ( iccer /= 0 ) GO TO 70
 lasblk = fcb( 6,idsn )
 IF ( lasblk >= nblock ) GO TO 7000
 fcb( 6,idsn ) = fcb( 6,idsn ) + 1
 GO TO 7000
 70   IF ( iccer /= 28 ) CALL dsmsg( 101 )
 IF ( idsn <= 21 .AND. idsn /= 8 .AND. idsn /= 9) GO TO 80
! ALLOW XPDT TO EXTEND (IDSN=9)---NOTE IDSN=8 IS THE NPTP
 itest = INDEX( mdsnam(8), 'ZAP' )
 IF ( idsn == 8  .AND. itest == 0 ) GO TO 80
 fcb( 7,ifilex ) = fcb( 6,ifilex )
 idsnr = idsn
 GO TO 10
 80   WRITE ( iwr, 902 )
 902   FORMAT(///,' NO MORE DISK SPACE AVAILABLE IN DEFAULT DIRECTORY',  &
     ' FOR PERMANENT FILES',/,' JOB ABORTED')
 CALL dsmsg( 122 )
 7000   RETURN
END SUBROUTINE dsgnwr
