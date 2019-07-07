gsql -g state_estimation loadflow_state_estimation.gsql
flg="loop"
iter=1
while [ "$flg" != '"No Bad Data Detected":"No Bad Data Detected"},' ]
do
	gsql -g state_estimation 'run query -dm state_estimation_weight_gain()'
	# gsql -g state_estimation 'run query -dm state_estimation_weight_solve("/home/tigergraph/output/",1,1,0,0.0001, 20)'   # flat start
	a=$(curl -X GET 'http://192.168.2.231:9000/query/state_estimation_weight_solve?outputfile="/home/tigergraph/output/"&flatstart=1&initial_Vm=1&initial_Va=0&tol=0.0001&IterLim=20')
	# flg=$(cut -d '{' -f4<<<$a)
	flg=$(echo $a|rev|cut -d '{' -f 2|rev)
	gsql -g state_estimation 'run query -dm state_estimation_BDI()'
	gsql -g state_estimation 'run query -dm SITS_R()'
	gsql -g state_estimation 'run query -dm delt_MNodes()'
	echo "$flg"
  iter+=1
  echo $iter
done