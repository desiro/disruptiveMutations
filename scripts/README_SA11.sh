#/bin/bash

################################################################
## paths
# tools
vsite=~/projects/rotavirus_SA11/vRNAsite.py
deopt=~/projects/rotavirus_SA11/disruptiveMutations.py
vcomp=~/projects/rotavirus_SA11/compare_vRNAsite.py

# data
data=~/projects/rotavirus_SA11/fasta
resVsite=~/projects/rotavirus_SA11/vRNAsite
resDeopt=~/projects/rotavirus_SA11/deopt
resVcomp=~/projects/rotavirus_SA11/comp


################################################################
## functions
# call vRNAsite
function vRNAsite() {
    settings=("-10.0 -clp -10.0 -vrh -vrs -vre -10.0 -plh -plv -pls 2.0 -bkh -bks 0.5" "-13.0 -clp -13.0 -mat" "-15.0 -clp -15.0 -mat" "-17.0 -clp -17.0 -mat")
    mutation=("$@")
    for m in "${mutation[@]}"
    do
        for s in "${settings[@]}"
        do
            echo nice -n 10 $vsite -pfx $resVsite/${n} -fsa $data/${m}.fa -thr 42 -ovr -nex peak -cit inter -cih -cis -cie $s
        done
    done
}

# call disruptiveMutations
function disruptiveMutations() {
    mutation=("$@")
    for m in "${mutation[@]}"
    do
        echo nice -n 10 $deopt -pfx $resDeopt/${m} -pss 52 -pse 354 -fsa $data/${m}.fa -thr 32 -slc 20 -ovr -cds -10.0
    done
}

# call compare_vRNAsite
function compare_vRNAsite() {
    mutation=$1
    peaks=(-10 -13 -15 -17)
    for p in "${peaks[@]}"
    do
        echo nice -n 10 $vcomp -pfx $resVcomp/${mutation}_comp_${p} -vst $resVsite/${mutation}/${mutation}_peak-10.0.tsv -ovr -cdp ${p}.0
    done
}


################################################################
## workflow
echo "# Step 1: Predicting interactions of SA11 wild-type with vRNAsite"
mut=(SA11_WT)
vRNAsite ${mut[@]}

echo "# Step 2: Creating mutants mutA,mutB,mutC,mutD,mutE and mutX"
mut=(mutA_seg1_seg10 mutB_seg2_seg9 mutC_seg3_seg8 mutD_seg4_seg7 mutE_seg5_seg6 mutX_seg1to10)
disruptiveMutations ${mut[@]}

echo "# Step 3: Predicting interactions of SA11 mutants with vRNAsite"
mut=(SA11_mutA SA11_mutB SA11_mutC SA11_mutD SA11_mutE SA11_mutX SA11_mutAtoE-mutX)
vRNAsite ${mut[@]}

echo "# Step 4: Compare interactions of SA11 mutants A to E, X and wild-type"
compare_vRNAsite SA11_mutAtoE-mutX

echo "# Step 5: Laboratory experiments on mutants A to E and X"

echo "# Step 6; Creating reverse complement of SA11 segment 4 and insert it into SA11 segment 11 mutant A"
#WT-S4:       862-875:GGAGGGCTGGGATA
#WT-S11:      231-244:ACTCACCAGTTTTT
#mutA-S11:    231-244:AaaCcCCttaTcTT
#rev-comp-S4:         tatcccagccctcc
#GGAGGGCUGGGAUA&ACUCACCAGUUUUU WT: -10.90 kcal/mol
#.(((((((((((..&..)).)))))))))
#GGAGGGCUGGGAUA&AaaCcCCuuaUcUU mutA: -7.70 kcal/mol
#.((((...(((...&....)))...))))
#GGAGGGCUGGGAUA&uaucccagcccucc mutA7: -28.80 kcal/mol
#((((((((((((((&))))))))))))))

echo "# Step 7: Predicting interactions of SA11 mutant A3 to A7 with vRNAsite"
mut=(SA11_mutAr SA11_mutA-mutAr-mutX)
vRNAsite ${mut[@]}

echo "# Step 8: Compare interactions of SA11 mutants A, Ar, X and wild-type"
compare_vRNAsite SA11_mutA-mutAr-mutX

echo "# Step 9: Laboratory experiments on mutant Ar"

