#!/bin/bash
#$ -S /bin/sh

for j in `cat lista_pdb`
  do
    echo $j
    echo $j > kk
    cp /home/julia/Documentos/Evo-idp/Evo-ordenadas/SCPE/Estructuras/$j.pdb .

    uno=$(cut -f1 -d"_" kk)
    dos=$(cut -f2 -d"_" kk)

    echo $uno
    echo $dos
    # sed "/PROTEIN CHAIN=/c\PROTEIN CHAIN= $dos" inputfile1.in > inputfile2.in
    # sed "s/var/"$dos"/g" inputfile1.in > inputfile2.in
    sed "s/XXXXXXXXXXX/"$j"/g" inputfile1.in > inputfile3.in
    sed "s/var/"$dos"/g" inputfile3.in > inputfile2.in
    rm inputfile3.in
    ./scpe.exe  inputfile2.in
    # #inputfile2.in
    #
    #
    mv Contact-information.dat /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/Contactos/Contact-$j.dat

    mkdir /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/OUT/$j/
    mv  /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/$j*.dat  /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/OUT/$j/
    mv  /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/$j*.Freq  /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/OUT/$j/
    mv  /home/julia/Documentos/Evo-idp/Evo_UNC/ORDENADAS/SCPE/$j*.Prob  /home/julia/Documentos/Evo-idp/Evo-ordenadas/SCPE/OUT/$j/

    rm $j.pdb



  done
