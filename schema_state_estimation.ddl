// 06/16/2017 [Chen Yuan] This schema is for state estimation
// 06/26/2017 [Chen Yuan] add an edge attribute "neighbors" to indicate if it is an edge of a triangle. 

drop all

set segsize_in_bits=8

// sort_rpi_cpi_matrix:for the L'U' and L"U" matrix HeapAccums (key based on rpi and cpi)
//                     - contains factorized LU values of B' and B" (both edge and node) 
typedef tuple <key1 int(8), cpi int(8), value1 double> sort_rpi_cpi_matrix
 
// sort_rpi_vertex:    for the L'U' and L"U" matrix vertex HeapAccums (key based on rpi) 
//                     - contains pointers, permutation and scaling info  
typedef tuple <key1 int(8), Lp int(8), Up int(8), rp int(8), cpi int(8), row_scaling double, col_scaling double> sort_rpi_vertex

// sort_id_vertex:     for the h array (key based on exId) 
//                     - contains information on the nodes and also for the pointer arrays 
typedef tuple <index1 int(8), value1 double> sort_id

//create vertex GNode (primary_id cid string, exId uint, flag uint, GenP double, GenQ double, LdP double, LdQ double, G double, B double, Vm double, Va double, Ri_V double, Ri_vP double, Ri_vQ double) with stats="OUTDEGREE_BY_EDGETYPE"
//the voltage, Vm, and angle, Va, here are state variables, not measurements

create vertex GNode (primary_id cid string, exId int, flag uint, M_P double, M_Q double, G double, B double, Vm double, Va double, P double, M_Vm double, M_Va double, Ri_V double, Ri_vP double, Ri_vQ double, busname string, M_id int) with stats="OUTDEGREE_BY_EDGETYPE"
//the voltage, Vm, and angle, Va, here are state variables, not measurements

create directed edge connected (from GNode, to GNode, G double, B double, hB double, K double, Kcount int, BIJ double, M_P_TLPF double, M_Q_TLPF double, neighbors int, deltaP double, deltaQ double, Ri_eP double, Ri_eQ double, reverse int, M_id int)

create vertex MNode (primary_id cid string, Mid int, residual double, residual_org double) with stats="OUTDEGREE_BY_EDGETYPE"

create directed edge Medge (from MNode, to MNode)

create graph state_estimation (GNode, connected, MNode, Medge)

EXPORT SCHEMA state_estimation

//clear graph store -HARD
//init graph store for graph company_graph
