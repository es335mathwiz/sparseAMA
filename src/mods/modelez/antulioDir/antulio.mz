MODEL > antulio

ENDOG >
         YA     0  _NOTD
         KA     0  _NOTD
         YB     0  _NOTD
         KB     0  _NOTD
         LAT    0  _NOTD
         AT     0  _NOTD
         SMALLA 0  _NOTD
         ONE    0  _DTRM


EQUATION > YA
EQTYPE   > IMPOSED
EQ       > YA = YAK* LAG(KA,1) +  YAE* AT +  YAB* YB + YAL* LAT

EQUATION > KA
EQTYPE   > IMPOSED
EQ       > KA = MAK* LAG(KA,1) + QA11* AT + QA12* YB + MAL* LAT

EQUATION > YB
EQTYPE   > IMPOSED
EQ       > YB = YBK* LAG(KB,1) + YBE* AT + YBA* LAG(YA,1)

EQUATION > KB
EQTYPE   > IMPOSED
EQ       > KB = MBK* LAG(KB,1) + MBE* AT + MBA* LAG(YA,1)

EQUATION > LAT
EQTYPE   > IMPOSED
EQ       > LAT - 1/MUA2* LEAD(LAT,1) = P11* LEAD(AT,1) + P21* AT 
                                     + P12* LEAD(YB,1) + P22* YB

EQUATION > AT
EQTYPE   > IMPOSED
EQ       > AT = RHO* LAG(AT,1) + SMALLA

EQUATION > SMALLA
EQTYPE >   IMPOSED
EQ >       SMALLA = 0*ONE

EQUATION > ONE
EQTYPE >   IMPOSED
EQ >       ONE = LAG(ONE,1)


END

