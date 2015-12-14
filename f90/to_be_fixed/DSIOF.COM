PARAMETER  ( maxpri = 80, maxfcb = 89 )
COMMON / dsio /  ieor  , ioerr , iprvop, iretrn, irword, idatad, &
                 idsn  , lcw   , lwords, maskh1, maskh2,         &
                 maske1, maske2, maske3, maske4, maxdsn, nwords, &
                 nbuff , ioblk , nbfz  , nlr   ,                 &
                 maskq1, maskq2, maskq3, maskq4, idsx  , idsp  , &
                 idsc  , idsrh , idsrt , idssb , idsse , idsch , &
                 idsct , idssh , idsst , idssd , idseb , idsef , &
                 iblock, lasnam, mcbmas, mulq1 , mulq2 , mulq3 , &
                 lhalf,  lendsp, lenwpb, nwrdel(4)
 COMMON /dsname/ mdsnam(maxfcb)
 COMMON /dsdevc/ numdev, dev(10)
 CHARACTER (LEN=2) :: dev
 CHARACTER (LEN=80) :: mdsnam
 COMMON /dbm   / idbbas, idbfre, idbdir, indbas, indclr, indcbp, &
                 nblock, lenalc, iocode, ifilex, name,   maxalc, &
                 maxblk, maxdsk, idblen, idbadr, ibasbf, inddir, &
                 numopn, numcls, numwri, numrea, lenopc
 INTEGER :: fcb
 COMMON / fcb  / fcb(17,maxfcb)
