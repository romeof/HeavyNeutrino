#!/bin/bash
#Specify needed variables
#int
#varType=int
#varTYPE=INT
#capLetter=I
#varList=(mu_charge ele_charge lep_charge)
#varLast=lep_charge
#varNum=mu_num
#varCount=mu_num
#double
varType=double
varTYPE=DOUBLE
capLetter=D
varList=(mu_3dip)
varLast=mu_3dip
varNum=mu_num
varCount=mu_num
#Print info
echo " "
#Declare variables
echo -e " $varType \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}[DEF_SIZE1D], \c"
  else
   echo "${varList[$pos]}[DEF_SIZE1D];"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  INIT_1DARRAY(${varList[$pos]},DEF_SIZE1D,DEF_VAL_$varTYPE);"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->Branch(\"${varList[$pos]}\", &${varList[$pos]}, \"${varList[$pos]}[$varNum]/$capLetter\");"
  let pos=pos+1
done
echo " "
#Set branch address
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->SetBranchAddress(\"${varList[$pos]}\", &${varList[$pos]});"
  let pos=pos+1
done
echo " "
#Analyzer
pos=0
for count in ${varList[@]}; 
do
  echo "tree->${varList[$pos]}[$varCount] = ;"
  let pos=pos+1
done
echo " "
