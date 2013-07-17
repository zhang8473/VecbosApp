#!/bin/bash
j=0
while read line
do
    (( j+=1 ))
    declare -a arr=($(echo $line | tr "," " "))
    List[$j]=`echo ${arr[0]}`
    XSec[$j]=${arr[1]}
    echo $j":"${List[$j]}"----"${XSec[$j]}
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
    outputname=${List[$j]#${List[$j]%/*}}
    outputname=/afs/cern.ch/work/z/zhangjin${outputname%.list}
    echo "#!/bin/bash">$outputname.sh
    echo "cd $CMSSW_BASE/src/">>$outputname.sh
    echo "eval \`scramv1 runtime -sh\`">>$outputname.sh
    echo "cd VecbosApp">>$outputname.sh
    if [ `expr match "${XSec[$j]}" ".*JSON.*"` != "0" ]; then
        echo "./VecbosApp ${List[$j]} $outputname.root --isData -json=${XSec[$j]}">>$outputname.sh		
    else
    	echo "./VecbosApp ${List[$j]} $outputname.root -weight=${XSec[$j]}">>$outputname.sh
    fi
    echo "rm $outputname.sh">>$outputname.sh
    chmod +x $outputname.sh
    echo $j":"${List[$j]}"----"${XSec[$j]}" ===>"
    case $1 in
	local)
		nohup $outputname.sh>$outputname.out&;;
        test)
		;;
	"")
    		bsub -q 1nw $outputname.sh;;
    esac
done
