PARAMETER  ( numfcb=89, numsof=10 )
COMMON / ginox / lginox, idslim, mdsfcb( 3,numfcb ), lensof( numsof )
 
