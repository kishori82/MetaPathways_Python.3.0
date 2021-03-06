#!/bin/bash


function extract  { 
~/pathway-tools/pathway-tools -api &
sleep 10
a=1 
for sample in ${samples[@]}
do
      perl extract_pathway_table_from_pgdb.pl -f $sample -out $sample.long.pwy.txt -t long

       a=$(( $a+1 ))
       f=$(( $a%10 ))
       if [ "$f" -eq 0 ]; then
          perl extract_pathway_table_from_pgdb.pl -exit 'a'
          sleep 10
          ~/pathway-tools/pathway-tools -api &
          sleep 10
       fi  
     
done
perl extract_pathway_table_from_pgdb.pl -exit 'a'

}




# script to extract long-pathway tables from each TARA ePGDB
samples=(
cenf
ceng
cenh
ceni
cenj
cenk
cenl
cenm
cenn
ceno
cenp
cenq
cenr
cens
cent
cenu
cenv
cenw
cenx
ceny
cenz
ceoa
ceoc
ceod
ceoe
ceof
ceog
ceoh
ceoi
ceoj
ceok
ceol
ceom
ceon
ceoo
ceop
ceoq
ceor
ceos
ceot
ceou
ceov
ceow
ceox
ceoy
ceoz
cepa
cepb
cepc
cepd
cepe
cepf
cepg
ceph
cepi
cepj
cepk
cepl
cepm
cepn
cepo
cepp
cepq
cepr
ceps
cept
cepu
cepv
cepw
cepx
cepy
cepz
ceqa
ceqb
ceqc
ceqd
ceqe
ceqf
ceqg
ceqh
ceqi
ceqj
ceqk
ceql
ceqm
ceqn
ceqo
ceqp
ceqq
ceqr
ceqs
ceqt
cequ
ceqv
ceqw
ceqx
ceqy
ceqz
cera
cerb
cerc
cerd
cere
cerf
cerg
cerh
ceri
cerj
cerk
cerl
cerm
cern
cero
cerp
cerq
cerr
cers
cert
ceru
cerv
cerw
cerx
cery
cerz
cesa
cesb
cesc
cesd
cese
cesf
cesg
cesh
cesi
cesj
cesk
cesl
cesm
cesn
ceso
cesp
cesq
cesr
cess
cest
cesu
cesv
cesw
cesx
cesy
cesz
ceta
cetb
cetc
cetd
cete
cetf
cetg
ceth
ceti
cetj
cetk
cetl
cetm
cetn
ceto
cetp
cetq
cetr
cets
cett
cetu
cetv
cetw
cetx
cety
cetz
ceua
ceub
ceuc
ceud
ceue
ceuf
ceug
ceuh
ceui
ceuj
ceuk
ceul
ceum
ceun
ceuo
ceup
ceuq
ceur
ceus
ceut
ceuu
ceuv
ceuw
ceux
ceuy
ceuz
ceva
cevb
cevc
cevd
ceve
cevf
cevg
cevh
cevi
cevj
cevk
cevl
cevm
cevn
cevo
cevp
cevq
cevr
cevs
cevt
cevu
cevv
cevw
cevx
cevy
cevz
cewa
cewb
cewc
cewe
cewf
cewg
cewh
cewi
cewj
cewk
cewo
cewp
cewq
cewr
cxwf
gov_assembly
tcf_18_dcm
tcf_18_srf
tcf_23_dcm
tcf_25_dcm
tcf_25_srf
tcf_30_dcm
tcf_31_srf
tcf_32_dcm
tcf_32_srf
tcf_34_dcm
tcf_34_srf
tcf_36_dcm
tcf_36_srf
tcf_37_mes
tcf_38_dcm
tcf_38_mes
tcf_38_srf
tcf_39_dcm
tcf_39_mes
tcf_39_srf
tcf_41_dcm
tcf_41_srf
tcf_42_dcm
tcf_42_srf
tcf_48_srf
tcf_52_dcm
tcf_64_dcm
tcf_64_srf
tcf_65_dcm
tcf_65_srf
tcf_66_dcm
tcf_66_srf
tcf_67_srf
tcf_68_dcm
tcf_68_mes
tcf_68_srf
tcf_70_mes
tcf_70_srf
tcf_72_dcm
tcf_72_mes
tcf_72_srf
tcf_76_dcm
tcf_76_mes
tcf_76_srf
tcf_78_dcm
tcf_78_mes
tcf_78_srf
tcf_82_dcm
tcf_85_dcm
tcf_109_dcm
tcf_109_srf
tcf_122_dcm
tcf_122_mes
tcf_122_srf
tcf_123_mix
tcf_123_srf
tcf_124_mix
tcf_124_srf
tcf_125_mix
tcf_125_srf
)

extract  
