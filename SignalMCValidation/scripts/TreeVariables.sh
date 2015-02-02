#!/bin/bash
#Specify needed variables
#int
varType=int
varTYPE=INT
capLetter=I
varList=(ele_charge)
varLast=ele_charge
varNum=ele_num
varCount=ele_num
#double
#varType=double
#varTYPE=DOUBLE
#capLetter=D
#varList=(ele_pt ele_eta ele_phi ele_en)
#varLast=ele_en
#varNum=ele_num
#varCount=gl
#Print info
echo " "
#Declare variables
echo -e " $varType \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}, \c"
  else
   echo "${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  ${varList[$pos]} = DEF_VAL_$varTYPE;"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo "  tree->Branch(\"${varList[$pos]}\", &${varList[$pos]}, \"${varList[$pos]}/$capLetter\");"
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
  echo "tree->${varList[$pos]} = ;"
  let pos=pos+1
done
echo " "
