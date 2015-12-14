COMMON / dosnam / dirtry, rfdir, INPUT, output, LOG , punch,       &
                  plot,   nptp , dic  , optp  , rdic, in12, out11, &
                  inp1,   inp2
CHARACTER (LEN=72) :: dirtry, rfdir, INPUT, output, LOG , punch
CHARACTER (LEN=72) :: plot  , nptp , dic  , optp  , rdic, in12, out11
CHARACTER (LEN=72) :: inp1,   inp2
COMMON / dsname / dsnames(89)
CHARACTER (LEN=80) :: dsnames
