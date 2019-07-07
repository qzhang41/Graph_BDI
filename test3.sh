gsql -g state_estimation loadflow_state_estimation.gsql
gsql -g state_estimation 'run query -dm state_estimation_weight_gain()'
gsql -g state_estimation 'run query -dm state_estimation_weight_solve("/home/tigergraph/output/",1,1,0,0.0001, 20)'   # flat start
gsql -g state_estimation 'run query -dm state_estimation_BDI()'
gsql -g state_estimation 'run query -dm SITS_R()'
gsql -g state_estimation 'run query -dm state_estimation_weight_solve("/home/tigergraph/output/",1,1,0,0.0001, 20)'   # flat start