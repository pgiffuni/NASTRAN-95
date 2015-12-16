SUBROUTINE dsiodd
 
 INCLUDE 'GINOX.COM'
 INCLUDE 'DSIOF.COM'
 lginox = 5*numfcb + numsof + 2
 lhalf  = 16
 lendsp = 0
 lenwpb = 0
 maskh1 = 'FFFF0000'X
 maskh2 = '0000FFFF'X
 maske1 = 'FF000000'X
 maske2 = '00FF0000'X
 maske3 = '0000FF00'X
 maske4 = '000000FF'X
 mcbmas = '40000000'X
 maxdsn = numfcb
 maskq1 = maske1
 maskq2 = maske2
 maskq3 = maske3
 maskq4 = maske4
 mulq1  = 2**24
 mulq2  = 2**16
 mulq3  = 2**8
 idsx   = '00EE0000'X
 idsp   = '000E0000'X
 idsc   = '000C0000'X
 idsrh  = '11000000'X
 idsrt  = '77000000'X
 idssb  = '22000000'X
 idsse  = '7F000000'X
 idsch  = '3B000000'X
 idsct  = '3F000000'X
 idssh  = '4B000000'X
 idsst  = '4F000000'X
 idssd  = 'DD000000'X
 idseb  = 'EB000000'X
 idsef  = 'EF000000'X
 nwrdel( 1 ) = 1
 nwrdel( 2 ) = 2
 nwrdel( 3 ) = 2
 nwrdel( 4 ) = 4
 RETURN
END SUBROUTINE dsiodd
