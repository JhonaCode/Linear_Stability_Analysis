#!/bin/bash
#help#
################################

current=`pwd`

name_run='teste'

para="Parameters_${name_run}"

DataIn=$current/src 

Dataout=$current/${name_run}

N=150

#(1)'spatial'    
#(2)'temporal'  
analysis="spatial"

#Derivative Matrix
#(1)'Fnitediff'
#(2)'Chebyshev'

diffM="Fnitediff"
#diffM="Chebyshev"

#(1)   = 'Morris_mxl'  
#(2)   = 'Hu_mixingl'  
#(3)   = 'Joncas_jet'  

base="Morris_mxl"
L=10.0

#Mapping
#(1)'__not'    
#(2)'point'    
#Only work with Chevisheb
#(3) 'tanmp'    
#(4) 'cylyn'  	 
#(5) 'sqrtm'        	

#mapp='point'  
#mapp='__not'  
##############
#mapp="tanmp" 
#mapp='cylyn' 
mapp='sqrtm' 

label="Icomp"


if [ ! -d "${Dataout}" ]; then
    mkdir -p "${Dataout}"
fi

sed  -e "s;%_%_%path;${Dataout};g" \
     -e "s;%_%_%analysis;${analysis};g" \
     -e "s;%_%_%diffM;${diffM};g" \
     -e "s;%_%_%base;${base};g" \
     -e "s;'%_%_%L';${L};g" \
     -e "s;'%_%_%N';${N};g" \
     -e "s;%_%_%mapp;${mapp};g" \
     -e "s;%_%_%label;${label};g" \
     ${DataIn}/Parameters.py > ${Dataout}/${para}.py

sed  -e "s;%_%_%Parameters;${para};g" \
     ${DataIn}/main.py > ${Dataout}/main.py

sed  -e "s;%_%_%Parameters;${para};g" \
     ${DataIn}/type_analysis.py > ${Dataout}/type_analysis.py

sed  -e "s;%_%_%Parameters;${para};g" \
     ${DataIn}/spectral.py > ${Dataout}/spectral.py

cp ${DataIn}/chebPy.py         ${Dataout}
cp ${DataIn}/mapping.py        ${Dataout}
cp ${DataIn}/Baseflow.py       ${Dataout}
cp ${DataIn}/Eigenvalues.py    ${Dataout}
cp ${DataIn}/diffusion.f90     ${Dataout}
cp ${DataIn}/Boundary.py       ${Dataout}
cp ${DataIn}/plotbaseflow.py   ${Dataout}
cp ${DataIn}/plotparameters.py ${Dataout}
cp ${DataIn}/diffusion.f90 ${Dataout}
cp ${DataIn}/diffusion.so ${Dataout}

cd ${Dataout}
#rm *.so
#rm *.o
#f2py -m diffusion -c  global.f90 diffusion.f90

python main.py
