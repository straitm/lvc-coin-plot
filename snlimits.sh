#!/bin/bash

cd ligopass2
if ! make; then
  exit 1
fi
cd -

ligopass2/ligopass2.sh # no args means compile

printf ' GW150914~\\cite{ligocat} & Untriggered & Bad         & ---\\phantom{0}    & ---\\phantom{0}  \\\\\n'
printf ' GW151012~\\cite{ligocat} & Untriggered & No data     & ---\\phantom{0}    & ---\\phantom{0}  \\\\\n'

parallel -k -j 4 << EOF
ligopass2/justt02   GW151226
ligopass2/justt02   GW170104
ligopass2/justt02   GW170608
ligopass2/justt02   GW170729
ligopass2/justt02   GW170809
ligopass2/justt02   GW170814
ligopass2/justt02   GW170817
ligopass2/justt02   GW170818
ligopass2/justt02   GW170823
EOF

printf '     \\citegcn{S190408an} & No data     & No data     &  ---\\phantom{0}   & ---\\phantom{0}  \\\\\n'

parallel -k -j 3 << EOF
ligopass2/justt02   S190412m
ligopass2/justt02   S190421ar
ligopass2/justt02   S190425z
ligopass2/t02ndligo S190426c   '44.1\,s'
ligopass2/justt02   S190503bf
ligopass2/justt02   S190510g
ligopass2/justt02   S190512at
ligopass2/t02ndligo S190513bm  '24.7\,s'
ligopass2/justt02   S190517h
ligopass2/justt02   S190519bj 
ligopass2/fdndligo  S190521g   '45.0\,s' '45.0\,s'
ligopass2/justt02   S190521r
EOF

parallel -k -j 2 << EOF
ligopass2/fdndligo  S190602aq  '45.0\,s' '45.0\,s'
ligopass2/fdndligo  S190630ag  '45.0\,s' '45.0\,s'
ligopass2/fdndligo  S190701ah  '45.0\,s' '45.0\,s'
ligopass2/fdndligo  S190706ai  '45.0\,s' '17.5\,s'
ligopass2/justt02   S190707q
EOF
