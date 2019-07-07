#gsql -g state_estimation drop query SE_BDI
#gsql -g state_estimation drop query add_MNodes
#gsql -g state_estimation drop query add_Medges
#gsql -g state_estimation drop query delt_MNodes

#gsql -g state_estimation add_MNodes.gsql
#gsql -g state_estimation delt_MNodes.gsql
#gsql -g state_estimation add_Medges.gsql
gsql -g state_estimation state_estimation_BDI.gsql
#gsql -g state_estimation SE_BDI.gsql