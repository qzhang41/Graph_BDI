gsql -g state_estimation loadflow_state_estimation.gsql
flg="loop"
iter=0
while [ "$flg" != '"flg_BDI":0}]}' ]
do
	gsql -g state_estimation 'run query -dm state_estimation_weight_gain()'
	gsql -g state_estimation 'run query -dm state_estimation_weight_solve("/home/tigergraph/output/",1,1,0,0.0001, 20)'   # flat start
	a=$(curl -X GET 'http://192.168.2.231:9000/query/state_estimation_BDI')
	flg=$(echo $a|rev|cut -d '{' -f 1|rev)
	#gsql -g state_estimation 'run query -dm state_estimation_BDI()'
	gsql -g state_estimation 'run query -dm SITS_R()'
	gsql -g state_estimation 'run query -dm delt_MNodes()'
	echo "$flg"
  iter+=1
  echo $iter
  echo $flg
done