SUBROUTINE dsgnrd
     INCLUDE 'XNSTRN.COM'
 INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 CHARACTER (LEN=4) :: cbuff(2)
 EQUIVALENCE     (cbuff,ibase)
 idsn   =  mdsfcb( 2, ifilex )
 idsnr  = idsn
 10   istrb  = fcb( 5, idsnr )
 IF ( nblock >= istrb ) GO TO 20
 idsnr  = mdsfcb( 3, idsnr ) / mulq2
 GO TO 30
 20   iend   = fcb( 6, idsnr )
 IF ( nblock <= iend ) GO TO 40
 idsnr  = IAND( mdsfcb( 3, idsnr  ), maskh2 )
 30   IF ( idsnr >= 1 .AND. idsnr <= maxdsn ) GO TO 10
 CALL dsmsg( 121 )
 40   IF ( idsn == idsnr ) GO TO 50
 CALL dsclos( idsn )
 mdsfcb( 1,idsn ) = IAND( mdsfcb( 1,idsn ), maskh1 )
 idsn   = idsnr
 mdsfcb( 1, idsn )   = ior ( mdsfcb( 1,idsn ), maskh2 )
 mdsfcb( 2, ifilex ) = idsn
 isave = iop
 iop = 0
 CALL dsopen( mdsnam(idsn), idsn, iop )
 iop = isave
 cbuff( indbas ) = mdsnam( idsn )
 50   ioblk  = nblock - istrb + 1
 CALL dsread( idsn, ibase(indbas+3), nbuff, ioblk )
 RETURN
END SUBROUTINE dsgnrd
