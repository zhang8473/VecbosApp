#!/bin/bash
#######################################################
# VecbosApp submitter based on Samples.csv            #
# Author: Zhang Jinzhong, zhangjin@cern.ch            #
#######################################################
#usage ./SubmitAll.sh [local,test]
j=0
while read line
do
    (( j+=1 ))
    declare -a arr=($(echo $line | tr "," " "))
    List[$j]=${arr[0]}
    XSec[$j]=${arr[1]}
    if [ `expr match "${XSec[$j]}" "[0-9./*()]*"` == ${#XSec[$j]} ]; then
	XSec[$j]=`bc<<< "scale=1;${XSec[$j]}"`
    fi
    NEvt[$j]=${arr[2]}
    echo "$j:${List[$j]} ---- ${XSec[$j]} / ${NEvt[$j]}"
done < "Samples.csv"
#Choice
until [ `expr match "$Choice" '[:blank: \.0-9]*'` = "${#Choice}" ]&&[[ "$Choice" == [0-9]* ]]
do
  echo -n "Which ones do you want to produce for? (Seperate by space, continuous ones by "..", like 1 2 4..17)"
  read Choice
done
#expand choices
for i in $Choice
do
  posdots=`expr index "$i" ..`
  if [ $posdots -gt 0 ]; then
      posdots=`expr $posdots - 1`
      numBegin=${i:0:$posdots}
      posdots=`expr $posdots + 2`
      numEnd=${i:$posdots:${#i}}
      Choice=${Choice/$i/$(seq $numBegin 1 $numEnd)}
  fi
done

for j in $Choice
do
    nsplitfiles=`eval echo ${List[$j]/"."/"_*."}|wc|awk '{print $2}'`
#    nsplitfiles=1
    for i in $(seq 1 1 $nsplitfiles)
    do
	if [ $nsplitfiles -gt 1 ]; then
	    inputname=${List[$j]/"."/"_$i."}
	else
	    inputname=${List[$j]}
	fi
	outputname=${inputname#${inputname%/*}}
	outputname=/afs/cern.ch/work/z/zhangjin${outputname%.list}
	echo "#!/bin/bash">$outputname.sh
	echo "cd $CMSSW_BASE/src/">>$outputname.sh
	echo "eval \`scramv1 runtime -sh\`">>$outputname.sh
	echo "cd VecbosApp">>$outputname.sh
	if [ `expr match "${XSec[$j]}" ".*JSON.*"` != "0" ]; then
            echo "./VecbosApp $inputname $outputname.root --isData -json=${XSec[$j]}">>$outputname.sh
	    echo "$j:$inputname --- X:${XSec[$j]} ===>"
	else
	    if [ $nsplitfiles -gt 1 ]; then
		NEvent=`bc<<< "scale=1;${NEvt[$j]}"`
    		echo "./VecbosApp $inputname $outputname.root -xsec=${XSec[$j]} -nevt=$NEvent">>$outputname.sh
		echo "$j:$inputname ---X:${XSec[$j]} ---N:$NEvent ---njobs:$nsplitfiles ===>"
	    else
    		echo "./VecbosApp $inputname $outputname.root -xsec=${XSec[$j]}">>$outputname.sh
		echo "$j:$inputname ---X:${XSec[$j]} ===>"
	    fi
	fi
	echo "rm $outputname.sh">>$outputname.sh
	chmod +x $outputname.sh
	case $1 in
	    local)
		nohup $outputname.sh>$outputname.out&;;
            test)
		;;
	    "")
    		bsub -q 2nd $outputname.sh;;
	esac
    done
done
