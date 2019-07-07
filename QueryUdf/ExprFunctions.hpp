/******************************************************************************
 * Copyright (c) 2015-2016, GraphSQL Inc.
 * All rights reserved.
 * Project: GraphSQL Query Language
 * udf.hpp: a library of user defined functions used in queries.
 *
 * - This library should only define functions that will be used in
 *   GraphSQL Query scripts. Other logics, such as structs and helper
 *   functions that will not be directly called in the GQuery scripts,
 *   must be put into "ExprUtil.hpp" under the same directory where
 *   this file is located.
 *
 * - Supported type of return value and parameters
 *     - int
 *     - float
 *     - double
 *     - bool
 *     - string (don't use std::string)
 *     - accumulators
 *
 * - Function names are case sensitive, unique, and can't be conflict with
 *   built-in math functions and reserve keywords.
 *
 * - Please don't remove necessary codes in this file
 *
 * - A backup of this file can be retrieved at
 *     <graphsql_root_path>/dev_<backup_time>/gdk/gsql/src/QueryUdf/ExprFunctions.hpp
 *   after upgrading the system.
 *
 ******************************************************************************/

#ifndef EXPRFUNCTIONS_HPP_
#define EXPRFUNCTIONS_HPP_

#define MAX_VERTICES 15000

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <algorithm> // Header file needed for ListAccum sorting
#include <gle/engine/cpplib/headers.hpp>
#include <tuple> // for accessing elements in tuple
#include <list>
#include <vector>

/**     XXX Warning!! Put self-defined struct in ExprUtil.hpp **
 *  No user defined struct, helper functions (that will not be directly called
 *  in the GQuery scripts) etc. are allowed in this file. This file only
 *  contains user-defined expression function's signature and body.
 *  Please put user defined structs, helper functions etc. in ExprUtil.hpp
 */
#include "ExprUtil.hpp"

//// ============= Header files needed for NISCLU =============
#include "graphlu.h"
#include "graphlu_util.h"
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include <unordered_set> // use for ILU filter
#include <iostream>
// #include <chrono>
// typedef std::chrono::high_resolution_clock Clock;
//// ==========================================================
//
//[Chen Yuan]
//// ============= Header files needed for self-defined functions =============
#include "get_chi.h"
// ==

typedef std::string string; //XXX DON'T REMOVE

/****** BIULT-IN FUNCTIONS **************/
/****** XXX DON'T REMOVE ****************/
inline int str_to_int (string str) {
  return atoi(str.c_str());
}

inline int float_to_int (float val) {
  return (int) val;
}

inline string ToString (double val) {
  char result[200];
  sprintf(result, "%g", val);
  return string(result);
}
/****************************************/

// [Chen Yuan] Achieve Chi Square Distribution Index
inline double chi_square_index(const double p, const double freedom_degree){
	return get_chi(p, freedom_degree);
}

// // [Chen Yuan] Testing
// inline void static_C_testing1(){
	// //list<LU_data> testing = new list<LU_data>();
	// testing.push_back(1,1,0.03);
	// std::cout << "teting " << testing << std::endl;
// }


// inline void static_C_testing2(){
	// std::cout << "testing " << testing << std::endl;
// }


//**********************************************************************************************************************************
// Created by: Chen Yuan, chen.yuan@geirina.net
// Date: 10/24/2017
// This code performs part of state estimation Gain matrix formulation. This version assigns gain matrix off-diagonal elements from mapaccum to listaccum for later system gain matrix built-up
// History: 
// 10/24/2017 [Chen Yuan] 
// **********************************************************************************************************************************
template <typename key, typename value,
          typename T_P,
		  typename T_Q>
inline void Map2List (MapAccum<key, value>& gMap_P_1step, MapAccum<key, value>& gMap_Q_1step, MapAccum<key, value>& gMap_P_2step, MapAccum<key, value>& gMap_Q_2step, ListAccum<T_P>& Gip, ListAccum<T_Q>& Giq) {
	
	int i;
	Gip.data_.reserve(Gip.data_.size()+gMap_P_1step.size()+gMap_P_2step.size());
	Giq.data_.reserve(Giq.data_.size()+gMap_Q_1step.size()+gMap_Q_2step.size());
	
	// struct timeval t2_st, t2_end, t3_st, t3_end; 
	// long seconds, useconds;
	// gettimeofday(&t2_st, 0);
	// gettimeofday(&t3_st, 0);
	//for(i=0; i<gMap.data_.size(); i++){
	for (const auto &map1: gMap_P_1step) {
		//Gip.data_.push_back(T(map1.first, map1.second));
		Gip.data_.emplace_back(T_P(map1.first, map1.second));
	}
	// gettimeofday(&t3_end, 0);
	// seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	// useconds = t3_end.tv_usec  - t3_st.tv_usec;
	// printf("\n\n============================================================================================== ");
	// std::cout << "Time1:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	// printf("\n================================================================================================ "); 
	
	// gettimeofday(&t3_st, 0);
	for (const auto &map2: gMap_Q_1step) {
		//Giq.data_.push_back(T(map2.first, map2.second));
		Giq.data_.emplace_back(T_Q(map2.first, map2.second));
	}
	// gettimeofday(&t3_end, 0);
	// seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	// useconds = t3_end.tv_usec  - t3_st.tv_usec;
	// printf("\n\n============================================================================================== ");
	// std::cout << "Time2:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	// printf("\n================================================================================================ "); 
	
	// gettimeofday(&t3_st, 0);
	
	for (const auto &map3: gMap_P_2step) {
		//Gip.data_.push_back(T(map3.first, map3.second));
		Gip.data_.emplace_back(T_P(map3.first, map3.second));
	}
	// gettimeofday(&t3_end, 0);
	// seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	// useconds = t3_end.tv_usec  - t3_st.tv_usec;
	// printf("\n\n============================================================================================== ");
	// std::cout << "Time3:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	// printf("\n================================================================================================ "); 
	
	// gettimeofday(&t3_st, 0);
		
	for (const auto &map4: gMap_Q_2step) {
		//Giq.data_.push_back(T(map4.first, map4.second));
		Giq.data_.emplace_back(T_Q(map4.first, map4.second));
	}	
	// gettimeofday(&t3_end, 0);
	// seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	// useconds = t3_end.tv_usec  - t3_st.tv_usec;
	// printf("\n\n============================================================================================== ");
	// std::cout << "Time3:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	// printf("\n================================================================================================ "); 
	
	
	// gettimeofday(&t2_end, 0);
	// seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	// useconds = t2_end.tv_usec  - t2_st.tv_usec;
	// printf("\n\n============================================================================================== ");
	// std::cout << "Total Time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	// printf("\n================================================================================================ "); 
	//}

	
}



// [Chen Yuan] gain matrix sorting
template <typename T_gainmatrix> //typename vertex_pp_comp>
inline string gain_matrix (ArrayAccum<ListAccum<T_gainmatrix>>& gGainMatrix) {
   struct timeval t1_st, t1_end; long seconds, useconds;
   
   gettimeofday(&t1_st, 0);
   
   
	// extern size_t g_si, g_sd, g_sp;
	printf("\n\n------------------------------------------------------------------------------------------- ");
	printf("\nStart Running gain_matrix function!\n");
    std::cout << "Gain Matrix Number of Rows:" << gGainMatrix.data_.size() << std::endl;
    printf("-------------------------------------------------------------------------------------------- \n\n");	
	// ------------------------------------------------------------------------------------------------
	// 				Initialize variables and arrays
	// ------------------------------------------------------------------------------------------------
	
 
	// Initialize arrays and variables
	//uint__t n, nnz, n_pp, nnz_pp, n_e, nnz_e;
	//int ret, i, j, p, iter;
	//real__t maxDeltaP=0, maxDeltaQ=0, max_change_ex=maxchange;
	//SGraphLU *graphlu, *graphlu_pp;
	uint__t n;
	int i;
  
    string result = "FAILED";
  
	//const double pi_value = 3.141592653589793;
    //const int64_t LUflag = 1; 
	// real__t *x, err;	
	// uint__t fctnnz; // number of nonzeros after the matrix is factorized

	// Get the dimension and the nnz of the matrix B' and B"
	n = gGainMatrix.data_.size(); // number of the row of gain matrix	//nnz=gBp; n_pp=gVertex.data_.size(); nnz_pp=gBpp;
    //n_e=gVertex.data_.size();  nnz_e/ = gYbus;	//nnz_e=gMatrix_all.data_.size(); // get the size of the Y bus matrix
  
	
	//// Convert vectors to GRAPHLU Data Structure (pointers)

	std::cout << " ======================== Initialization of ararys used to store the factorized matrix and LU ========================"<< std::endl;

	struct timeval t2_st, t2_end, t3_st, t3_end;

  //initialization to 0
  // [tc] use memset here
  //memset(Vm, 0, sizeof(real__t)*(n));
  //memset(Va, 0, sizeof(real__t)*(n));
  //memset(deltaP, 0, sizeof(real__t)*(n));
  //memset(deltaQ, 0, sizeof(real__t)*(n));
	
  // =================================== Sort all input HeapAccum =================================
  gettimeofday(&t3_st, 0);
/*  gettimeofday(&t2_st, 0);
  // Input HeapAccums are assumed to be unsorted. Sort them before assigning values to
  // local arrays
  gVertex.sort_heap();
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort gVertex HeapAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");   */
  
  struct ListAccum_sort{
            inline bool operator() (const T_gainmatrix& tuple1, const T_gainmatrix& tuple2)
            {
                return (tuple1.index < tuple2.index);
            }
	};
  
  
  gettimeofday(&t2_st, 0);  
  
  for(i=0;i<gGainMatrix.data_.size(); ++i){
	std::sort(gGainMatrix.data_[i].begin(), gGainMatrix.data_[i].end(), ListAccum_sort());
  }
  
  //std::sort(gMatrix_all.data_.begin(), gMatrix_all.data_.end(), ListAccum_sort());
 	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort gMatrix_all ListAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	std::cout << "gMatrix_all Size:: " << gGainMatrix.data_.size() << std::endl;
	printf("\n================================================================================================ ");   
    
  // get the time for sorting all HeapAccums 
	gettimeofday(&t3_end, 0);
	seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	useconds = t3_end.tv_usec  - t3_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort all HeapAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ "); 
	

EXIT:
	
	
  printf("\n\n----------------------------------------------------------------------------------------");
	printf("\t\t End of Running Gain Matrix function!");
	printf("\t\t ----------------------------------------------------------------------------------------\n\n");	
	
	gettimeofday(&t1_end, 0);
	seconds=t1_end.tv_sec  - t1_st.tv_sec; 
	useconds = t1_end.tv_usec  - t1_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Total Time of ExpressFunction " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");  
 
  return result;
  
  

  
}




template <typename T>
inline void check_symetric(ListAccum<T>& gGainP, ListAccum<T>& gGainQ, int64_t& count_P, int64_t& count_Q)
{
	//int count_P = 0;
	//int count_Q = 0;
	int i = 0;
	struct Gain_sort
	{
            inline bool operator() (const T& tuple1, const T& tuple2)
            {
                return (tuple1.value1 <= tuple2.value1);
            }
	};
	
    for(i=0;i<gGainP.data_.size(); ++i){
	   std::sort(gGainP.data_.begin(), gGainP.data_.end(), Gain_sort());
    }
	
	for(i=0;i<gGainQ.data_.size(); ++i){
	   std::sort(gGainQ.data_.begin(), gGainQ.data_.end(), Gain_sort());
    }
	
	for(i=0;i<gGainP.data_.size(); ++i){
		if (gGainP.data_[i].value1 != gGainP.data_[i+1].value1){
			count_P = count_P + 1;
			++i;
		}
	}
	
	for(i=0;i<gGainQ.data_.size(); ++i){
		if (gGainQ.data_[i].value1 != gGainQ.data_[i+1].value1){
			count_Q = count_Q + 1;
			++i;
		}
	}
	
	std::cout << "count_P: " << count_P << std::endl;
	std::cout << "count_Q: " << count_Q << std::endl;
}
// **********************************************************************************************************************************
//Created by Qiwei
//2019 summer intern 
inline double BDI_getER (uint64_t& index, const ArrayAccum<SumAccum<double>> &gVa) {
  return gVa.data_[index];
} 

inline double BDI_getK (const int64_t index, const ArrayAccum<SumAccum<double>> &gVa) {
  return gVa.data_[index];
} 

inline double store_S (ArrayAccum<SumAccum<double>>& S_x, int m_c)
{
real__t *ax;
uint__t *ai;
uint__t *ap;
ax = (real__t *)malloc(sizeof(real__t)*(m_c*m_c)); // values in B' 
ai = (uint__t *)malloc(sizeof(uint__t)*(m_c*m_c)); // column indices of B'
ap = (uint__t *)malloc(sizeof(uint__t)*(m_c+1)); // initial row pointers
for (int i= 0; i < m_c*m_c; i++)
	{
		ai[i] = i % m_c;
		ax[i] = S_x.data_[i];
	}
for (int i= 0; i < m_c; i++)
	{
		ap[i] = i *m_c;
	}
	ap[m_c] = m_c*m_c;
	int m_size = m_c*m_c;
	SingletonMatrixInterface::deleteS();
	CSRMatrix *S = SingletonMatrixInterface::getS();
	// S->loadMatrix(n, nnz, CSR_p, CSR_i, CSR_x);
	S->loadMatrix(m_c, m_size, ap, ai, ax);
free(ax);
free(ai);
free(ap);
	return 1;
}



// **********************************************************************************************************************************
//Created by Qiwei
//2019 summer intern 
inline double cal_SR(ListAccum<ListAccum<double>> R_c, ListAccum<ListAccum<double>> R_i, ArrayAccum<SumAccum<double>> &cor_residual, int m_c)
{
	
double threshhold = 3;
//obtain S matrix
	vector<vector<double>> S_N_all;
	S_N_all.resize(m_c, vector<double>(m_c, 0));
	CSRMatrix *S_N = SingletonMatrixInterface::getS();	
	double *lx; 
	lx = S_N->getX();
	for(int i=0;i<m_c;i++)
	{
		for(int j=0;j<m_c;j++)
		{
			S_N_all[i][j] = lx[i*m_c+j];
		}
	}
//cal error	
	vector<double> residual;
	vector<double> residual_idx;
	vector<double> org_residual;
	vector<double> org_residual_idx;
	vector<double> residual_copy;
	vector<double> residual_idx_copy;
	vector<int> error_idx_set;
	vector<int> error_set;
	int inte_idx = 0;
	int s_c = 0;
	for(int i=0;i<R_c.size();i++)
	{
		s_c+=1;
		std::cout<<"s_C  "<<s_c<<std::endl;
		residual.clear();
		residual_idx.clear();
		org_residual.clear();
		org_residual_idx.clear();
		residual_copy.clear();
		residual_idx_copy.clear();
		
		residual.resize(R_c.data_[i].size(), 0);
		org_residual.resize(R_c.data_[i].size(), 0);
		org_residual_idx.resize(R_c.data_[i].size(), 0);
		for(int j=0;j<R_c.data_[i].size();j++)
		{
			if(R_i.data_[i].get(j)!=0)
			{
				residual[j] = abs(R_c.data_[i].get(j));
				org_residual[j] = abs(R_c.data_[i].get(j));
				org_residual_idx[j] = R_i.data_[i].get(j);			
			}
		}
		sort( residual.begin(), residual.end() );
		residual.erase( unique( residual.begin(), residual.end() ), residual.end() ); // sort the residual from small to large, delete the same M in group 
		residual_idx.resize(residual.size(), 0);
		for(int k=0;k<residual.size();k++)
		{
			// std::cout << "residual_value2: "<< residual[k]<< std::endl;	 
			for(int l=0;l<org_residual.size();l++)
			{
				if (residual[k]==org_residual[l])
				{
					residual_idx[k] = org_residual_idx[l];                   // find the corresponding M idx according to residual
					// std::cout << "residual_idx3: "<< residual_idx[k]<< std::endl;	 
					break;
				}
			} 
		}
		for(int k=0;k<residual.size();k++)
		{
			std::cout << "residual_value1  : "<< residual[k]<< std::endl;	 
			std::cout << "residual_idx1  : "<< residual_idx[k]<< std::endl;	 
		}
		// copy the value of residual and idx
		residual_copy.resize(R_c.data_[i].size(), 0);
		residual_idx_copy.resize(R_c.data_[i].size(), 0);
		for(int k=0;k<residual.size();k++)
		{
			residual_copy[k] = residual[k];
			residual_idx_copy[k] = residual_idx[k];
		}
		
		
		// cal error 
		error_idx_set.clear();
		error_idx_set.reserve(50);
		bool flg = false;
		int max_M_idx = residual_idx[residual_idx.size()-1];
		double max_error = residual[residual.size()-1]/S_N_all[max_M_idx-1][max_M_idx-1];
		error_idx_set.push_back(max_M_idx);
		for(int k = 0; k<error_idx_set.size();k++)
		{
			std::cout << "error_idx_set[i]1   "<<error_idx_set[k] <<std::endl;	
		}
		residual.pop_back();
		residual_idx.pop_back();
		for(int p = 0; p<residual.size(); p++)//except the last one(max one)
		{
			std::cout << "residual22222 "<<residual[p]<< std::endl;
			// std::cout << "p "<<p<< std::endl;
			inte_idx = residual_idx[p];
			// std::cout << "inte_idx "<<inte_idx<< std::endl;
			residual[p] = residual[p] - abs(S_N_all[inte_idx-1][max_M_idx-1] * max_error);
			std::cout << "residual33333 "<<residual[p]<< std::endl;
			// std::cout << "S_N_all "<<S_N_all[inte_idx-1][max_M_idx-1]<< std::endl;
			// std::cout << "max_error "<<max_error<< std::endl;	
		}
		sort(residual.begin(), residual.end());
		if(residual.size()==0)
		{
		    cor_residual.data_[max_M_idx-1] = max_error;
			continue;
		}
		if (residual[residual.size()-1] < threshhold)
		{
			 std::cout << "single error "<< std::endl;
			 cor_residual.data_[max_M_idx-1] = max_error;
			 std::cout << "cor_residual.data_[max_M_idx-1]"<<cor_residual.data_[max_M_idx-1]<< std::endl;	
			 std::cout << "max_M_idx"<<max_M_idx<< std::endl;	
			 // return max_error, max_M_idx;
		} else{
			 std::cout << "multiple error "<< std::endl;
			 while(flg!= true)
			 {
				for(int k=0;k<residual.size();k++)
				{
					residual_idx.clear();
					for(int l=0;l<org_residual.size();l++)
					{
						if (residual[k]==org_residual[l])
						{
							residual_idx[k] = org_residual_idx[l];                   // find the corresponding M idx according to residual
							break;
						}
					} 
				}
				// for(int k=0;k<residual.size();k++)
				// {
					// std::cout << "residual_idx[k]   "<<residual_idx[k] <<std::endl;	
				// }
				max_M_idx = residual_idx[residual.size()-1];
				// std::cout << "max_M_idx2   "<<max_M_idx <<std::endl;	
				// std::cout << "residual_idx.size()   "<<residual.size() <<std::endl;	
				error_idx_set.push_back(max_M_idx);
				
				for(int k = 0; k<error_idx_set.size();k++)
				{
					std::cout << "error_idx_set[i]2   "<<error_idx_set[k] <<std::endl;	
				}
//////////////////////////////////////////////////////////////////////////////////////////////solve the LP accord to error_idx_set (copy from solve and factorization)
				uint__t n_s, nnz_s;
				int ret; // p, iter;
				int error_p = 0;
				SGraphLU *graphlu_s;
				real__t *ax_s;
				uint__t *ai_s;
				uint__t *ap_s;
				n_s = error_idx_set.size();
				nnz_s = n_s *n_s;
				ax_s = (real__t *)malloc(sizeof(real__t)*(nnz_s)); // values in B' 
				ai_s = (uint__t *)malloc(sizeof(uint__t)*(nnz_s)); // column indices of B'
				ap_s = (uint__t *)malloc(sizeof(uint__t)*(n_s+1)); // initial row pointers

				for(int q = 0; q<n_s; q++)
				{
					for(int j = 0; j<n_s; j++)
					{
						int idx_1 = error_idx_set[q];
						int idx_2 = error_idx_set[j];
						ax_s[q*n_s+j] = abs(S_N_all[idx_1-1][idx_2-1]);
						std::cout << "idx_1"<<idx_1 <<std::endl;		
						std::cout << "idx_2"<<idx_2 <<std::endl;								
						std::cout << "ax_s[i]"<<ax_s[q*n_s+j] <<std::endl;							
					}
				}
				for(int q = 0; q<nnz_s;q++)
				{
						ai_s[q] = q % n_s;
						std::cout << "a1_s[q]"<<ai_s[q] <<std::endl;			
				}
				for (int q= 0; q < n_s; q++)
				{
					ap_s[q] = q * n_s;
					std::cout << "ap_s[q]"<<ap_s[q] <<std::endl;				
				}
				ap_s[n_s] = nnz_s;
				graphlu_s = (SGraphLU *)malloc(sizeof(SGraphLU));
				ret = GraphLU_Initialize(graphlu_s);
				ret = GraphLU_CreateMatrix(graphlu_s, n_s, nnz_s, ax_s, ai_s, ap_s);
				ret = GraphLU_Analyze(graphlu_s);
				ret = GraphLU_CreateScheduler(graphlu_s);
				if (!ret){ 
				  GraphLU_CreateThreads(graphlu_s, 2, TRUE);
				  GraphLU_BindThreads(graphlu_s, FALSE);
				  error_p = GraphLU_Factorize_MT(graphlu_s);
				if (error_p < 0)
				  std::cout << "GraphLU_Factorize error code:" << error_p <<std::endl;  
				}   
				else{
				  error_p = GraphLU_Factorize(graphlu_s);
				  printf("Factorization time: %.8g\n", graphlu_s->stat[1]);
				  				std::cout << "stop where5"<< std::endl;
				if (error_p < 0) //there is an error, print the code 
				  std::cout << "GraphLU_Factorize error code:" << error_p <<std::endl; 
				}
				real__t *lx_s, *ux_s; 
				uint__t *li_s, *ui_s; 
				size_t *lp_s, *up_s; 
				uint__t *rp_s, *cp_s, *rpi_s, *cpi_s, nnz_Ls, nnz_Us; //row (column) permutation, row (column) permutation inverse
				real__t *rows_s, *cols_s; //*ldiag, *cscale,
			    lx_s = ux_s = NULL; 
				li_s = ui_s = NULL; 
				lp_s = up_s = NULL; 

				int row, col, k;
				int LUrow, LUcol;
    
				GraphLU_DumpLU(graphlu_s, &lx_s, &li_s, &lp_s, &ux_s, &ui_s, &up_s); 
				rp_s = graphlu_s->row_perm; // row permutation, permuted row # --> original row #  (original row)
				rpi_s = graphlu_s->row_perm_inv;  // row permutation inverse, original row # --> permuted row # (permuted row)
				cp_s = graphlu_s->col_perm;
				cpi_s = graphlu_s->col_perm_inv;
				rows_s = graphlu_s->row_scale; 
				cols_s = graphlu_s->col_scale_perm;
				nnz_Ls = graphlu_s->l_nnz;
				nnz_Us = graphlu_s->u_nnz;
											
				uint__t n, nnz_p_L, nnz_p_U, n_e, nnz_e;  // nnz_pp_L, nnz_pp_U,
				const double pi_value = 3.141592653589793;

				n =n_s;	
				nnz_p_L=nnz_Ls; nnz_p_U=nnz_Us;
				
				double *r_set, *e_set; 

				r_set = (double *)malloc(sizeof(double)*(n+1)); 

				e_set = (double *)malloc(sizeof(double)*(n+1)); 

				double sum, diag;
				double *b;
				int q, j, p; // ret, iter
				b = (double *)malloc(sizeof(double)*(n));

				memset(r_set, 0, sizeof(double)*(n));
				memset(e_set, 0, sizeof(double)*(n));

				for(int q = 0; q<n_s; q++)
				{
					std::cout<<"error_idx_set[q]2 "<<error_idx_set[q]<<std::endl;
					for(int k=0;k<org_residual.size();k++)
					{
						if (org_residual_idx[k]==error_idx_set[q])
						{
							r_set[q] = org_residual[k];                   // find the corresponding M idx according to residual
							std::cout<<"r_set[q]"<<r_set[q]<<std::endl;
							break;
						}
					}
				}
				for (q=0; q<n; ++q){
				  b[q]=r_set[rp_s[q]]*rows_s[rp_s[q]];
				}
				// Full forward substitution 
				for (q=0; q<n; ++q)
				{
					sum=0.;
					diag=0.;
					for(p=lp_s[q]; p< lp_s[q+1]; ++p){
						if (q!=li_s[p]){ // if it is not diagnal element
							j=li_s[p];
							sum += lx_s[p] * b[j];
						}
						else{
							diag = lx_s[p];
						}
					}
					b[q] -= sum; 
					b[q] /= diag; 
				} 
				
				// Full backward substitution 
				for (q=n-1; q>=0; --q)
				{
					sum=0.;
					for(p=up_s[q]; p< up_s[q+1]; ++p){
						if (q!=ui_s[p]){ // if it is not diagnal element
							sum += ux_s[p] * b[ui_s[p]]; 
						}
					}
			   
					b[q] -= sum; 
				}      
				// Update V angle (Va)
				for (q=0; q<n; ++q)
				{ 
					e_set[q] += b[cpi_s[q]]*cols_s[cpi_s[q]];    // x(k+1) = x(k) + delta_x 
					std::cout << "e_set"<<e_set[q]<< std::endl;		
				}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// error_set = ;           
				//update cor_residual accord to error_set
				for(int q = 0; q<error_idx_set.size();q++)
				{
					cor_residual.data_[error_idx_set[q]-1] = e_set[q];
				}
				residual.pop_back();
				residual_idx.pop_back();
				if(residual.size()==0)
				{
					break;
				}
				for(int p = 0; p<residual.size(); p++)//except the last one(max one)
				{
					inte_idx = residual_idx[p];
					residual[p] = residual[p] - abs(S_N_all[inte_idx-1][max_M_idx-1] * max_error);
					std::cout << "S_N_all[inte_idx-1][max_M_idx-1] "<<S_N_all[inte_idx-1][max_M_idx-1] << std::endl;
					std::cout << "max_error"<<max_error<< std::endl;	
					std::cout << "residual[p]_after"<<residual[p]<< std::endl;	
				}
				sort(residual.begin(), residual.end());
				if (residual[residual.size()-1] < threshhold)
				{
					flg = true;
					std::cout << "flg = true;"<< std::endl;	
				}else{
					std::cout << "flg = false;"<< std::endl;	
				}
				
				
				EXIT:
					GraphLU_Destroy(graphlu_s);
					free(ax_s);
					free(ai_s);
					free(ap_s);
					free(graphlu_s);
					free(lx_s);
					free(li_s);
					free(lp_s);
					free(ux_s);
					free(ui_s);
					free(up_s);
					
					free(e_set);
					free(r_set);
					free(b);
			 }
		}

		
	}

	return 1;
}
// **********************************************************************************************************************************
//Created by Qiwei
//2019 summer intern 
// perform matrix multiplication
inline double matrix_KH_cal (int m_c , int n_b,ArrayAccum<SumAccum<double>>K_all, ArrayAccum<SumAccum<double>>H_t_all, ArrayAccum<SumAccum<double>>& S_all)
{
	vector<vector<double>> K_all_2;
	vector<vector<double>> H_t_all_2;
	vector<vector<double>> S_all_2;
	K_all_2.resize(n_b, vector<double>(m_c, 0));
	H_t_all_2.resize(n_b, vector<double>(m_c, 0));
	S_all_2.resize(m_c, vector<double>(m_c, 0));
	double inte_s = 0;
	// for(int i=0;i<n_b;i++)
	// {
		// for(int j=0;j<m_c;j++)
		// {
			// K_all_2[i][j] = K_all.data_[i*(m_c)+j];
			// H_t_all_2[i][j] = H_t_all.data_[i*(m_c)+j];
		// }
	// }
	for(int i=0;i<m_c;i++)
	{
		for(int j=0;j<m_c;j++)
		{
			inte_s = 0;
			std::cout<<"inte_s"<<std::endl;
			for(int k = 0; k<n_b;k++)
			{	
				// inte_s = inte_s + K_all_2[k][j]*H_t_all_2[k][i];
				inte_s = inte_s + K_all.data_[k*(m_c)+j]*H_t_all.data_[k*(m_c)+i];
			}
			if(i == j)
			{
				S_all.data_[i*m_c+j] = 1 - inte_s;
			}else{
				S_all.data_[i*m_c+j] = (-1) * inte_s;
			}
		}
	}
}

//**********************************************************************************************************************************
// Created by: Chen Yuan, chen.yuan@geirina.net
// Date: 08/01/2017
// This code performs state estimation Gain matrix factorization. This version performs LU factorization on the Gain_P and Gain_Q matrix
// After performing LU factorization, the resulting B' and B" LU factors are stored in ListAccums. 
// History: 
// 08/01/2017 [Chen Yuan] Added the GainMatrix_factorize C/C++ function to implement LU decomposition for matrices Gain_P and Gain_Q
// **********************************************************************************************************************************
template <typename T> 
          //typename T_GainQ, 
          //typename T_GainP_L_matrix, typename GainP_L_matrix_comp,
          //typename T_GainP_U_matrix, typename GainP_U_matrix_comp, 
          //typename T_GainQ_L_matrix, typename GainQ_L_matrix_comp,
          //typename T_GainQ_U_matrix, typename GainQ_U_matrix_comp,
          //typename T_vertex_GainP, typename vertex_GainP_comp,
          //typename T_vertex_GainQ, typename vertex_GainQ_comp>
inline string GainMatrix_factorize_static (int64_t &gslackbus, int64_t &gnnz_GainP, int64_t &gnnz_GainQ, ArrayAccum<SumAccum<int64_t>>& gGainP_p, ArrayAccum<SumAccum<int64_t>>& gGainQ_p, ArrayAccum<ListAccum<T>>& gGainP, ArrayAccum<ListAccum<T>>& gGainQ)
              //HeapAccum<T_GainP_L_matrix, GainP_L_matrix_comp>& gMatrix_GainP_L, HeapAccum<T_GainP_U_matrix, GainP_U_matrix_comp>& gMatrix_GainP_U, HeapAccum<T_vertex_GainP, vertex_GainP_comp>& gVertex_GainP,
              //HeapAccum<T_GainQ_L_matrix, GainQ_L_matrix_comp>& gMatrix_GainQ_L, HeapAccum<T_GainQ_U_matrix, GainQ_U_matrix_comp>& gMatrix_GainQ_U, HeapAccum<T_vertex_GainQ, vertex_GainQ_comp>& gVertex_GainQ) 
			  {

	// extern size_t g_si, g_sd, g_sp;
	printf("\n\n------------------------------------------------------------------------------------------- ");
	printf("\n Start Running GainMatrix_factorize_static function!\n");
	std::cout << "GainP Number of rows:" << gGainP.data_.size() << ",\tNumber of nonzeros:" << gnnz_GainP << std::endl;
	std::cout << "GainQ Number of rows:" << gGainQ.data_.size() << ",\tNumber of nonzeros:" << gnnz_GainQ << std::endl;
	//std::cout << "Y bus Number of rows:" << gVertex.data_.size() << ",\tNumber of nonzeros:" << gMatrix_all.data_.size() <<------------------------------------------------------------ \n\n");	
	// ------------------------------------------------------------------------------------------------
	// 				Initialize variables and arrays
	// ------------------------------------------------------------------------------------------------
	
 
	// Initialize arrays and variables
	uint__t n_GainP, nnz_GainP, n_GainQ, nnz_GainQ, slackbus; //n_e, nnz_e;
	int ret, i, j; // p, iter;
	SGraphLU *graphlu_GainP, *graphlu_GainQ;
  
	string result = "FAILED";

	// Get the dimension, n, and the nnz of the matrix GainP and GainQ
	n_GainP = gGainP.data_.size();	nnz_GainP = gnnz_GainP; n_GainQ = gGainQ.data_.size(); nnz_GainQ = gnnz_GainQ;
	slackbus = gslackbus;
	//n_e=gVertex.data_.size(); 	nnz_e=gMatrix_all.data_.size(); // get the size of the Y bus matrix
	
	std::cout << "GainP size: " << n_GainP << std::endl;
	std::cout << "GainQ size: " << n_GainQ << std::endl;
	
	//// Convert vectors to GRAPHLU Data Structure (pointers)
	real__t *ax_GainP, *ax_GainQ; //*deltaP, *deltaQ, *Vm, *Va; *Pn, *Qn; *eG, *eB,
	uint__t *ai_GainP, *ai_GainQ; // *ei, *piv; *btype; 
	//uint__t *ap, *ap_pp, *ep;
    uint__t *ap_GainP, *ap_GainQ; 	
 	
	ax_GainP = (real__t *)malloc(sizeof(real__t)*(nnz_GainP)); // values in B' 
	ai_GainP = (uint__t *)malloc(sizeof(uint__t)*(nnz_GainP)); // column indices of B'
	ap_GainP = (uint__t *)malloc(sizeof(uint__t)*(n_GainP+1)); // initial row pointers
	// b = (real__t *)malloc(sizeof(real__t)*(n)); // b in the Ax=b

	ax_GainQ = (real__t *)malloc(sizeof(real__t)*(nnz_GainQ)); // values in B"
	ai_GainQ = (uint__t *)malloc(sizeof(uint__t)*(nnz_GainQ)); // column indices of B"
	ap_GainQ = (uint__t *)malloc(sizeof(uint__t)*(n_GainQ+1)); // initial row pointers

  
	std::cout << " ======================== Initialization of ararys used to store the factorized matrix and LU ========================"<< std::endl;

	struct timeval t2_st, t2_end, t3_st, t3_end; long seconds, useconds;
	
	// =================================== Sort both GainP and GainQ =================================
	gettimeofday(&t3_st, 0);
	gettimeofday(&t2_st, 0);
	// Input ArrayAccums' ListAccums are assumed to be unsorted. Sort them before assigning values to local arrays
    struct Gain_sort
	{
            inline bool operator() (const T& tuple1, const T& tuple2)
            {
                return (tuple1.index1 < tuple2.index1);
            }
	};
	
    for(i=0;i<gGainP.data_.size(); ++i){
	   std::sort(gGainP.data_[i].begin(), gGainP.data_[i].end(), Gain_sort());
    }
	
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort GainP ArrayAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ "); 
	
	
	gettimeofday(&t2_st, 0); 
	
    // struct GainQ_sort
	// {
            // inline bool operator() (const T_GainQ& tuple3, const T_GainQ& tuple4)
            // {
                // return (tuple3.index < tuple4.index);
            // }
	// };
	
     for(i=0;i<gGainQ.data_.size(); ++i){
	    std::sort(gGainQ.data_[i].begin(), gGainQ.data_[i].end(), Gain_sort());
     }
	
 	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort GainQ ArrayAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ "); 
    
	// get the time for sorting all HeapAccums 
	gettimeofday(&t3_end, 0);
	seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	useconds = t3_end.tv_usec  - t3_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to sort both GainP and GainQ ArrayAccum:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ "); 
	
	// ========================================= Convert ap_GainP, ap_GainQ =========================================
	gettimeofday(&t2_st, 0); 
	int ii = 1; // counter for the array  
	// Assign value to the first position of the pointer array 
	ap_GainP[0] = 0;
	ap_GainQ[0] = 0;
	//ep[0] = 0;
	for (i=0;i<gGainP_p.size();i++){
		ap_GainP[ii] = gGainP_p.data_[i] + ap_GainP[ii-1];
		//std::cout<< ii << ",ap_GainP,"<< ap_GainP[ii] <<  std::endl;
		ii++;
	}
	
	ii = 1;
	for (i=0;i<gGainQ_p.size();i++){
		ap_GainQ[ii] = gGainQ_p.data_[i] + ap_GainQ[ii-1];
		//std::cout<< ii << ",ap_GainQ,"<< ap_GainQ[ii] <<  std::endl;
		ii++;
	}
 
 	// ========================================= Convert GainP (ai_GainP, ax_GainP) and GainQ (ai_GainQ, ax_GainQ) =========================================
	int i_GainP = 0, i_GainQ = 0;
	for (i=0; i<gGainP.data_.size(); ++i) {
		for (ii=0; ii<gGainP.data_[i].size(); ++ii) {
	//for (ii=0; ii<gGainP.size(); ++ii) {
    //ei[ii] = gMatrix_all.data_[ii].ei;
    //eG[ii] = gMatrix_all.data_[ii].eG;
    //eB[ii] = gMatrix_all.data_[ii].eB; 
		//if(gMatrix_all.data_[ii].Bp_x != 0)
		//{
			if (gGainP.data_[i].get(ii).index1 < slackbus){
				ai_GainP[i_GainP] = gGainP.data_[i].get(ii).index1 - 1;
			}
			else if (gGainP.data_[i].get(ii).index1 > slackbus){
				ai_GainP[i_GainP] = gGainP.data_[i].get(ii).index1 - 2;
			}
			ax_GainP[i_GainP] = gGainP.data_[i].get(ii).value1;
      //std::cout<< i_GainP << ",ax,"<< ax_GainP[i_GainP] << ",ai,"<< ai_GainP[i_GainP] << std::endl;
			i_GainP++;
		}
	}
    //}
    //if(gMatrix_all.data_[ii].Bpp_x != 0)
    //{
	for (i=0; i<gGainQ.data_.size(); ++i) {
		for (ii=0; ii<gGainQ.data_[i].size(); ++ii) {
			ai_GainQ[i_GainQ] = gGainQ.data_[i].get(ii).index1 - 1;
			ax_GainQ[i_GainQ] = gGainQ.data_[i].get(ii).value1;
			//std::cout<< i_GainQ << ",ax,"<< ax_GainQ[i_GainQ] << ",ai,"<< ai_GainQ[i_GainQ] << std::endl;
			i_GainQ++;
		}
	} 	
	// Done converting all input HeapAccums to arrays
  
	// Get the time to convert the data between containers and arrays
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	printf("\n\n============================================================================================== ");
	std::cout << "Time to convert data to GRAPHLU time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");		
	
	

	
// ------------------------------------------ [Debug] print out gain matrix and check --------------------------------------------
// -------------------------------------------------------------------------------------------------------------------

	//// Convert vectors to GRAPHLU Data Structure (pointers)
	// real__t *GainP, *GainQ; //*deltaP, *deltaQ, *Vm, *Va; *Pn, *Qn; *eG, *eB,
	
	// GainP = (real__t *)malloc(sizeof(real__t)*(n_GainP*n_GainP)); 
	// GainQ = (real__t *)malloc(sizeof(real__t)*(n_GainQ*n_GainQ));

	// int temp_ap;
	// for (i=0; i<n_GainP; ){
		// temp_ap = ap_GainP[i];
		// //std::cout << i << ", ap: " << ap_GainP[i] << std::endl;
		// for(ii=0; ii<n_GainP; ) {
			// if (ii == ai_GainP[temp_ap]) {
				// GainP[i*n_GainP + ii] = ax_GainP[temp_ap];
				// temp_ap++;
			// } 
			// else {
				// GainP[i*n_GainP + ii] = 0;
			// }
			// std::cout << GainP[i*n_GainP + ii] << ",";
			// ii++;
		// }
		// std::cout << "" << std::endl;
		// i++;
	// }
	
	
	
	// ----------------------------------------------------------------------------------------------------
	// 								Call GRAPHLU and Factorize GainP Matirx
	// ----------------------------------------------------------------------------------------------------
	int error_p = 0, error_pp = 0;
	gettimeofday(&t2_st, 0);
	std::cout << "\n ======================== Start factorizing GainP Matirx ======================== \n"<< std::endl;
	// Initialize the structure for GRAPHLU
	graphlu_GainP = (SGraphLU *)malloc(sizeof(SGraphLU));
	ret = GraphLU_Initialize(graphlu_GainP);
	std::cout<<"Initialize: " << ret << std::endl; 
	
	// Create Matrix
	ret = GraphLU_CreateMatrix(graphlu_GainP, n_GainP, nnz_GainP, ax_GainP, ai_GainP, ap_GainP);
	// graphlu->cfgf[0] = 1.;
	// // Set control parameters
	//graphlu->cfgi[1] = 0; // 0 is no MC64 scaling
	// 					 // 1 is with MC64 scaling
	std::cout<<"CreateMatrix: " << ret << std::endl; 

	// Analyze (column/row ordering, scaling. For details please refer to GRAPHLU user's guide)
	GraphLU_Analyze(graphlu_GainP);
	printf("analysis time: %.8g\n", graphlu_GainP->stat[0]);
	
	// creates the task scheduler for parallel LU factorization. If you want to run 
	// parallel factorization or parallel re-factorization, it should be called after GraphLU Analyze.
	ret = GraphLU_CreateScheduler(graphlu_GainP);
	printf("time of creating scheduler: %.8g\n", graphlu_GainP->stat[4]);
	printf("suggestion: %s.\n", ret==0?"parallel":"sequential");
	//std::cout << "GraphLU MC64 scaling 4th:" << graphlu->cfgi[1] <<std::endl; 
	// This function creates threads for parallel computation. The second argument (thread)
	// specifies the number of threads, including the main thread. The last argument (check)
	// specifies whether to check the number of threads or not.
  if (!ret){ // parallel factorization 
	  GraphLU_CreateThreads(graphlu_GainP, 2, TRUE);
	  printf("total cores: %d, threads created: %d\n", (int)(graphlu_GainP->stat[9]), (int)(graphlu_GainP->cfgi[5]));
	
	  // This function binds threads to cores (unbind = FALSE) or unbinds threads from cores (unbind = TRUE).
	  GraphLU_BindThreads(graphlu_GainP, FALSE);
	
	  // Numerical LU factorization with partial pivoting, parallel
	  error_p = GraphLU_Factorize_MT(graphlu_GainP);
	  printf("factorization time: %.8g\n", graphlu_GainP->stat[1]);
    if (error_p < 0) //there is an error, print the code 
      std::cout << "GraphLU_Factorize error code:" << error_p <<std::endl;  
  }   
  else{  // Sequential factorization
	  error_p = GraphLU_Factorize(graphlu_GainP);
	  printf("Factorization time: %.8g\n", graphlu_GainP->stat[1]);
    if (error_p < 0) //there is an error, print the code 
      std::cout << "GraphLU_Factorize error code:" << error_p <<std::endl; 
	}
 
  real__t *lx_GainP, *ux_GainP; 
  uint__t *li_GainP, *ui_GainP; 
  size_t *lp_GainP, *up_GainP; 
  
 	//uint__t fctnnz_GainP; // number of nonzeros after the matrix is factorized
	//size_t l_nnz_GainP, u_nnz_GainP; // number of nonzeros in L and U factors
  // get the permutation arrays, please refer to graphlu.h for details
  uint__t *rp_GainP, *cp_GainP, *rpi_GainP, *cpi_GainP, nnz_LGainP, nnz_UGainP; //row (column) permutation, row (column) permutation inverse
  real__t *rows_GainP, *cols_GainP; //*ldiag, *cscale,
  //int__t *pivot, *pivot_inv;
  // rpi = (real__t *)malloc(sizeof(real__t)*(n));
  
  lx_GainP = ux_GainP = NULL; 
  li_GainP = ui_GainP = NULL; 
  lp_GainP = up_GainP = NULL; 

  int row, col, k;
  int LUrow, LUcol;
    
  GraphLU_DumpLU(graphlu_GainP, &lx_GainP, &li_GainP, &lp_GainP, &ux_GainP, &ui_GainP, &up_GainP); 
     
   // Get the number of nonzeros in the factorized matrix
 	//fctnnz = graphlu->lu_nnz;
 	//l_nnz = graphlu->l_nnz;
 	//u_nnz = graphlu->u_nnz;
 
 	std::cout << " ======================== Number of Total nonzeros after factorization is: "<< graphlu_GainP->lu_nnz << "========================"<<std::endl;
 	std::cout << " ======================== Number of nonzeros in L: "<< graphlu_GainP->l_nnz << "========================"<<std::endl;
 	std::cout << " ======================== Number of nonzeros in U: "<< graphlu_GainP->u_nnz << "========================"<<std::endl;
 	
   // get the permutation arrays, rp and cp may be different if the matrix is not symmetric
  rp_GainP = graphlu_GainP->row_perm; // row permutation, permuted row # --> original row #  (original row)
  rpi_GainP = graphlu_GainP->row_perm_inv;  // row permutation inverse, original row # --> permuted row # (permuted row)
  cp_GainP = graphlu_GainP->col_perm;
  cpi_GainP = graphlu_GainP->col_perm_inv;

  //ldiag = graphlu->ldiag; // diagnal elements
  //cscale = graphlu->cscale; // sclaing factors
  rows_GainP = graphlu_GainP->row_scale; 
  cols_GainP = graphlu_GainP->col_scale_perm;
  nnz_LGainP = graphlu_GainP->l_nnz;
  nnz_UGainP = graphlu_GainP->u_nnz;
  //pivot = graphlu->pivot; 
  //pivot_inv = graphlu->pivot_inv; // pivoting array  
	
	gettimeofday(&t3_st, 0);
	std::cout << "\n ======================== Save L_GainP and U_GainP ======================== \n"<< std::endl;

	
	// free MatrixGainP_L() and MatrixGainP_U();
	SingletonMatrixInterface::deleteMatrixGainP_U();
	SingletonMatrixInterface::deleteMatrixGainP_L();
	
	// static testing [Chen Yuan]  
	CSRMatrix *L_GainP = SingletonMatrixInterface::getMatrixGainP_L();
	CSRMatrix *U_GainP = SingletonMatrixInterface::getMatrixGainP_U();
	L_GainP->loadMatrix(n_GainP, nnz_LGainP, lp_GainP, li_GainP, lx_GainP, rp_GainP, cpi_GainP, rows_GainP, cols_GainP);
	U_GainP->loadMatrix(n_GainP, nnz_UGainP, up_GainP, ui_GainP, ux_GainP);
	// static testing end

	gettimeofday(&t3_end, 0);
	seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	useconds = t3_end.tv_usec  - t3_st.tv_usec;
	
	printf("\n\n============================================================================================== ");
	std::cout << "Time for Save L_GainP and U_GainP:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");


    
// ------------------------------------------ [Debug] print out and check --------------------------------------------
// -------------------------------------------------------------------------------------------------------------------

 	// for (i=0; i<n_GainP; ++i){ // check the permutation arrays
 	 	// std::cout << i << ",rp_GainP," << rp_GainP[i] << ",rpi_GainP," << rpi_GainP[i] << ",cp_GainP," << cp_GainP[i] << ",cpi_GainP," << cpi_GainP[i] << std::endl;
 	// } 
  
  //std::cout << " ======================== Get L Factor ========================"<<std::endl;
 	/*for (i=0; i<l_nnz; ++i){
     std::cout<< i << ",lx,"<< lx[i] << ",li,"<< li[i] << std::endl; 
  }
   
  for (i=0; i<n+1; ++i){
     std::cout<< i << ",lp,"<< lp[i] << std::endl;
  }
   
  //std::cout << "\n ======================== Get U Factor ========================"<<std::endl;
 	for (i=0; i<u_nnz; ++i){
     std::cout<< i << ",ux,"<< ux[i] << ",ui,"<< ui[i] << std::endl;
  }
   
  for (i=0; i<n+1; ++i){
 		std::cout<< i << ",up,"<< up[i] << std::endl;
 	}*/

 	// check arrays
  //for (i=0; i<n; ++i){
 	// 	std::cout<< i << ",\t diag:"<< ldiag[i] << ",\t cscale:"<< cscale[i] << ",\t rows:"<< rows[i] << ",\t cols:"<< cols[i] << std::endl;
  //}
  //for (i=0; i<n; ++i){
 	// 	std::cout<< i << ",rp," << rp[i]  << ",cpi," << cpi[i] << ",rows,"<< rows[i] << ",cols," << cols[i] <<std::endl;
  //}  
 
 	// // multiply U factor by cols (which is )
 	// for (i=0; i<n; ++i){
 	// 	for (j=up[i]; j<up[i+1]; ++j){
 	// 		p=rp[piv[ui[j]]]; // this may be rp!!!!!!!!!!!
 	// 		ux[j] = ux[j]/cols[p]; // include the scaling factor into L factor
 	// 	}
 	// }
 
 	// ================================= Get the pivoting vectors =================================
 	//for (i=0; i<n; ++i){
 	//	std::cout<< i << ",\t pivot:"<< pivot[i] << ",\t pivot_inv:"<< pivot_inv[i] << std::endl;
 	//}

  
    //get the factorized LU values and find the row and col before permutation
    // for (i = 0; i < n_GainP; ++i){
      // row = rp_GainP[i];
      // LUrow = row+1;
      
      // gVertex_GainP += T_vertex_GainP(i+1, lp_GainP[i+1]-lp_GainP[i], up_GainP[i+1]-up_GainP[i], rp_GainP[i], cpi_GainP[i], rows_GainP[i], cols_GainP[i]); 
      // // process all non-zeros in L matrix on row i
      // for (j = lp_GainP[i]; j < lp_GainP[i+1]; ++j){
        // col = cp_GainP[li_GainP[j]];             
        // LUcol = col+1;
        
        // gMatrix_GainP_L += T_GainP_L_matrix((i+1)*100000+(li_GainP[j]+1), li_GainP[j], lx_GainP[j]); 
      // } 
         
      // // process all non-zeros in U matrix on row i
      // for (j = up_GainP[i]; j < up_GainP[i+1]; ++j){
        // col = cp_GainP[ui_GainP[j]];
        // LUcol = col+1; 
        
        // gMatrix_GainP_U += T_GainP_U_matrix((i+1)*100000+(ui_GainP[j]+1), ui_GainP[j], ux_GainP[j]);          
      // }   
    // }
    // gVertex_GainP.sort_heap();
    // gMatrix_GainP_L.sort_heap();
    // gMatrix_GainP_U.sort_heap();

  
  //gVertex_p.sort_heap();
  //gMatrix_p_L.sort_heap();
  //gMatrix_p_U.sort_heap();
	// ----------------------------------------------------------------------------------------------------
	// 								Call GRAPHLU and Factorize GainQ Matirx
	// ----------------------------------------------------------------------------------------------------
	std::cout << "\n ======================== Start factorizing GainQ Matirx ======================== \n"<< std::endl;
	// Initialize the structure for GRAPHLU
	graphlu_GainQ = (SGraphLU *)malloc(sizeof(SGraphLU));
	GraphLU_Initialize(graphlu_GainQ);
 
	// Create Matrix
	GraphLU_CreateMatrix(graphlu_GainQ, n_GainQ, nnz_GainQ, ax_GainQ, ai_GainQ, ap_GainQ);
	// graphlu_pp->cfgf[0] = 1.;
	
	// // Set control parameters
	//graphlu_pp->cfgi[1] = 0; // 0 is no MC64 scaling
	// 					 // 1 is with MC64 scaling
 
	// Analyze (column/row ordering, scaling. For details please refer to GRAPHLU user's guide)
	GraphLU_Analyze(graphlu_GainQ);
	printf("analysis time: %.8g\n", graphlu_GainQ->stat[0]);
  
	// creates the task scheduler for parallel LU factorization. If you want to run 
	// parallel factorization or parallel re-factorization, it should be called after GraphLU Analyze.
	ret = GraphLU_CreateScheduler(graphlu_GainQ);
	printf("time of creating scheduler: %.8g\n", graphlu_GainQ->stat[4]);
	printf("suggestion: %s.\n", ret==0?"parallel":"sequential");
	
	// This function creates threads for parallel computation. The second argument (thread)
	// specifies the number of threads, including the main thread. The last argument (check)
	// specifies whether to check the number of threads or not.
  if (!ret){ // parallel factorization
	  GraphLU_CreateThreads(graphlu_GainQ, 2, TRUE);
	  printf("total cores: %d, threads created: %d\n", (int)(graphlu_GainQ->stat[9]), (int)(graphlu_GainQ->cfgi[5]));
	
	  // This function binds threads to cores (unbind = FALSE) or unbinds threads from cores (unbind = TRUE).
	  GraphLU_BindThreads(graphlu_GainQ, FALSE);
	
	  // Numerical LU factorization with partial pivoting, parallel
	  error_pp = GraphLU_Factorize_MT(graphlu_GainQ);
	  printf("factorization time: %.8g\n", graphlu_GainQ->stat[1]);
    if (error_pp < 0) //there is an error, print the code 
      std::cout << "GraphLU_Factorize error code:" << error_pp <<std::endl;      
  }
  else{  //Sequential factorization 
	  error_pp = GraphLU_Factorize(graphlu_GainQ);
    printf("Factorization time: %.8g\n", graphlu_GainQ->stat[1]);
    if (error_pp < 0) //there is an error, print the code 
      std::cout << "GraphLU_Factorize error code:" << error_pp <<std::endl; 
  }

  real__t *lx_GainQ, *ux_GainQ; 
  uint__t *li_GainQ, *ui_GainQ; 
  size_t *lp_GainQ, *up_GainQ; 
  //uint__t fctnnz_pp; // number of nonzeros after the matrix is factorized
  //size_t l_nnz_pp, u_nnz_pp; // number of nonzeros in L and U factors
  // get the permutation arrays, please refer to graphlu.h for details
  uint__t *rp_GainQ, *cp_GainQ, *rpi_GainQ, *cpi_GainQ, nnz_LGainQ, nnz_UGainQ; //row (column) permutation, row (column) permutation inverse
  real__t *rows_GainQ, *cols_GainQ;  //*ldiag_pp, *cscale_pp,
  //int__t *pivot_pp, *pivot_inv_pp; 
  
  lx_GainQ = ux_GainQ = NULL; 
  li_GainQ = ui_GainQ = NULL; 
  lp_GainQ = up_GainQ = NULL; 
  
  
  GraphLU_DumpLU(graphlu_GainQ, &lx_GainQ, &li_GainQ, &lp_GainQ, &ux_GainQ, &ui_GainQ, &up_GainQ); 
   
   // Get the number of nonzeros in the factorized matrix
 	//fctnnz_pp = graphlu_pp->lu_nnz;
 	//l_nnz_pp = graphlu_pp->l_nnz;
 	//u_nnz_pp = graphlu_pp->u_nnz;
 
 	std::cout << " ======================== Number of Total nonzeros after factorization is: "<< graphlu_GainQ->lu_nnz << "========================"<<std::endl;
 	std::cout << " ======================== Number of nonzeros in L: "<< graphlu_GainQ->l_nnz << "========================"<<std::endl;
 	std::cout << " ======================== Number of nonzeros in U: "<< graphlu_GainQ->u_nnz << "========================"<<std::endl;
 	
  // get the permutation arrays, rp and cp may be different if the matrix is not symmetric
  rp_GainQ = graphlu_GainQ->row_perm;
  rpi_GainQ = graphlu_GainQ->row_perm_inv;
  cp_GainQ = graphlu_GainQ->col_perm;
  cpi_GainQ = graphlu_GainQ->col_perm_inv;
  
  //ldiag_pp = graphlu_pp->ldiag; // diagnal elements
  //cscale_pp = graphlu_pp->cscale; // sclaing factors
  rows_GainQ = graphlu_GainQ->row_scale; 
  cols_GainQ = graphlu_GainQ->col_scale_perm;
  nnz_LGainQ = graphlu_GainQ->l_nnz;
  nnz_UGainQ = graphlu_GainQ->u_nnz;
  //pivot_pp = graphlu_pp->pivot; 
  //pivot_inv_pp = graphlu_pp->pivot_inv; // pivoting array   
  
	gettimeofday(&t3_st, 0);
	std::cout << "\n ======================== Save L_GainQ and U_GainQ ======================== \n"<< std::endl;
  	// static testing [Chen Yuan]  
	
	// free MatrixGainP_L() and MatrixGainP_U();
	SingletonMatrixInterface::deleteMatrixGainQ_U();
	SingletonMatrixInterface::deleteMatrixGainQ_L();
	
	
	CSRMatrix *L_GainQ = SingletonMatrixInterface::getMatrixGainQ_L();
	CSRMatrix *U_GainQ = SingletonMatrixInterface::getMatrixGainQ_U();
	L_GainQ->loadMatrix(n_GainQ, nnz_LGainQ, lp_GainQ, li_GainQ, lx_GainQ, rp_GainQ, cpi_GainQ, rows_GainQ, cols_GainQ);
	U_GainQ->loadMatrix(n_GainQ, nnz_UGainQ, up_GainQ, ui_GainQ, ux_GainQ);
	// static testing end
	
	gettimeofday(&t3_end, 0);
	seconds=t3_end.tv_sec  - t3_st.tv_sec; 
	useconds = t3_end.tv_usec  - t3_st.tv_usec;
	
	printf("\n\n============================================================================================== ");
	std::cout << "Time for Save L_GainQ and U_GainQ:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");
  
  
  
    
// ------------------------------------------ [Debug] print out and check --------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	
 	// for (i=0; i<n_GainQ; ++i){ // check the permutation arrays
 	 	// std::cout << "rp_GainQ," << rp_GainQ[i] << ",rpi_GainQ," << rpi_GainQ[i] << ",cp_GainQ,"<< cp_GainQ[i] << ",cpi_GainQ," << cpi_GainQ[i] << std::endl;
 	// } 
  
  //std::cout << " ======================== Get L Factor ========================"<<std::endl;
 	/*for (i=0; i<l_nnz_pp; ++i){
    std::cout<< i << ",lx_pp,"<< lx_pp[i] << ",li_pp,"<< li_pp[i] << std::endl; 
  }
   
  for (i=0; i<n+1; ++i){
    std::cout<< i << ",lp_pp,"<< lp_pp[i] << std::endl;
  }
   
   //std::cout << "\n ======================== Get U Factor ========================"<<std::endl;
 	for (i=0; i<u_nnz_pp; ++i){
    std::cout<< i << ",ux_pp,"<< ux_pp[i] << ",ui_pp,"<< ui_pp[i] << std::endl;
  }
   
 	for (i=0; i<n+1; ++i){
 		std::cout<< i << ",up_pp,"<< up_pp[i] << std::endl;
 	}*/
 
 	// check arrays
  //for (i=0; i<n; ++i){
 	// 	std::cout<< i << ",\t diag_pp:"<< ldiag_pp[i] << ",\t cscale_pp:"<< cscale_pp[i] << ",\t rows_pp:"<< rows_pp[i] << ",\t cols_pp:"<< cols_pp[i] << std::endl;
  //}
  //for (i=0; i<n; ++i){
 	// 	std::cout<< i << ",rp_pp," << rp_pp[i]  << ",cpi_pp," << cpi_pp[i] << ",rows_pp,"<< rows_pp[i] << ",cols_pp," << cols_pp[i] <<std::endl;
  //}   
 
 	// // multiply U factor by cols (which is )
 	// for (i=0; i<n; ++i){
 	// 	for (j=up[i]; j<up[i+1]; ++j){
 	// 		p=rp[piv[ui[j]]]; // this may be rp!!!!!!!!!!!
 	// 		ux[j] = ux[j]/cols[p]; // include the scaling factor into L factor
 	// 	}
 	// }
 
 	// ================================= Get the pivoting vectors =================================
 	//for (i=0; i<n; ++i){
 	//	std::cout<< i << ",\t pivot_pp:"<< pivot_pp[i] << ",\t pivot_inv_pp:"<< pivot_inv_pp[i] << std::endl;
 	//}
  
      	
    //get the factorized LU values and find the row and col before permutation
    // for (i = 0; i < n_GainQ; ++i){
      // row = rp_GainQ[i];
      // LUrow = row+1;
      
      // gVertex_GainQ += T_vertex_GainQ(i+1, lp_GainQ[i+1]-lp_GainQ[i], up_GainQ[i+1]-up_GainQ[i], rp_GainQ[i], cpi_GainQ[i], rows_GainQ[i], cols_GainQ[i]); 
      // // process all non-zeros in L matrix on row i
      // for (j = lp_GainQ[i]; j < lp_GainQ[i+1]; ++j){
        // col = cp_GainQ[li_GainQ[j]];
        // LUcol = col+1;
        
        // gMatrix_GainQ_L += T_GainQ_L_matrix((i+1)*100000+(li_GainQ[j]+1), li_GainQ[j], lx_GainQ[j]); 
      // } 
      
      // // process all non-zeros in U matrix on row i
      // for (j = up_GainQ[i]; j < up_GainQ[i+1]; ++j){
        // col = cp_GainQ[ui_GainQ[j]];
        // LUcol = col+1;
        
        // gMatrix_GainQ_U += T_GainQ_U_matrix((i+1)*100000+(ui_GainQ[j]+1), ui_GainQ[j], ux_GainQ[j]);
      // }    
    // }
    // gVertex_GainQ.sort_heap();
    // gMatrix_GainQ_L.sort_heap();
    // gMatrix_GainQ_U.sort_heap();          

  
  // Get the time for factorizing B' and B"
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;
	
	printf("\n\n============================================================================================== ");
	std::cout << "Time for GRAPHLU to factorize GainP and GainQ:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");
  
  if(error_p >= 0 && error_pp >= 0){
      result = "Factorization Completed";
      printf("\n\n============================================================================================== ");
  	  std::cout << "Factorization Completed" << std::endl;
  	  printf("\n================================================================================================ "); 
    //continue to do the iterations only if there is no error from factorization
//  	printf("\n\n============================================================================================== ");
//  	std::cout << "\n Start iteratively updating deltaP and deltaQ " << std::endl;
//  	printf("\n================================================================================================\n ");
//  


  } else {
    //factorization failed
 	  result = "Factorization FAILED";
    printf("\n\n============================================================================================== ");
  	std::cout << "Factorization FAILED" << std::endl;
  	printf("\n================================================================================================ ");  
  }
EXIT:
	// std::cout << "\n free1 " << std::endl;
	GraphLU_Destroy(graphlu_GainP);
	GraphLU_Destroy(graphlu_GainQ);
	// std::cout << "\n free2 " << std::endl;
	free(ax_GainP);
	free(ai_GainP);
	free(ap_GainP);
	free(ax_GainQ);
	free(ai_GainQ);
	free(ap_GainQ);
	free(graphlu_GainP);
	free(graphlu_GainQ);
	// std::cout << "\n free3 " << std::endl;
		
	//free(eG);
	//free(eB);
	//free(ei);
	//free(ep);

	//free(deltaP);
	//free(deltaQ);
	//free(Vm);
	//free(Va);		
	//free(Pn);
	//free(Qn);	
  
 	free(lx_GainP);
	free(li_GainP);
	free(lp_GainP);
    free(ux_GainP);
	free(ui_GainP);
	free(up_GainP);
	
 	//free(rp_GainP);
	//free(cp_GainP);
	//free(rpi_GainP);
 	//free(cpi_GainP);
	//free(rows_GainP);
	//free(cols_GainP); 
	

   	
 	free(lx_GainQ);
	free(li_GainQ);
	free(lp_GainQ);
 	free(ux_GainQ);
	free(ui_GainQ);
	free(up_GainQ);   
 	
	//free(rp_GainQ);
	//free(cp_GainQ);
	//free(rpi_GainQ);
 	//free(cpi_GainQ);
	//free(rows_GainQ);
	//free(cols_GainQ); 	

	
  printf("\n\n----------------------------------------------------------------------------------------");
	printf("\t\t End of Running GainMatrix_factorize_static C function!");
	printf("\t\t ----------------------------------------------------------------------------------------\n\n");	
 
  return result;
}





// **********************************************************************************************************************************
 // Created by: Chen Yuan, chen.yuan@geirina.net
 // Date: 11/14/2017
 // This code performs to solve state estimation GainP * deltaVa = H_P' * r_P. This version does NOT perform LU factorization. GainP is assumed to be previously factorized and
 // the resulting L and U matrices are taken as inputs. The equation is solved with forward/backward substitution.
 // History: 
 // 08/24/2017: [Chen Yuan] GainP matrix has been previously factorized. Built this function to solve GainP * deltaVa = H_P'*r_P, by employing forward/backward substitution, and update voltage angle, Va. 
 
//template <//typename T_deltaP_vertex, typename deltaP_vertex_comp,
          //typename T_vertex_Va,  
          //typename GainP_L_matrix, typename GainP_L_matrix_comp,  
          //typename GainP_U_matrix, typename GainP_U_matrix_comp,
          //typename GainP_LU_vertex, typename GainP_LU_vertex_comp>
//inline string SE_solve_GainP (ArrayAccum<SumAccum<double>>& gVa, ArrayAccum<SumAccum<double>>& gH_r_P,
              //HeapAccum<GainP_L_matrix, GainP_L_matrix_comp>& gGainP_L, HeapAccum<GainP_U_matrix, GainP_U_matrix_comp>& gGainP_U, HeapAccum<GainP_LU_vertex, GainP_LU_vertex_comp>& gGainP_LU_vertex, double& gmax_dVa) {
inline string SE_solve_GainP_static (ArrayAccum<SumAccum<double>>& gVa, ArrayAccum<SumAccum<double>>& gH_r_P, double& gmax_dVa) {

	struct timeval t1_st, t1_end;
	gettimeofday(&t1_st, 0); 
	//printf(xxx.c_str());
	
	// extern size_t g_si, g_sd, g_sp;
	printf("\n\n------------------------------------------------------------------------------------------- ");
	printf("\nStart Running SE_solve_GainP_static function!\n");	
	printf("-------------------------------------------------------------------------------------------- \n\n");	
	// ------------------------------------------------------------------------------------------------
	// 				Initialize variables and arrays
	// ------------------------------------------------------------------------------------------------


	CSRMatrix *L_GainP = SingletonMatrixInterface::getMatrixGainP_L();
	CSRMatrix *U_GainP = SingletonMatrixInterface::getMatrixGainP_U();

		
	
	// For L' and U'
	double *lx, *ux, *rows, *cols; 
	uint *li, *ui, *rp, *cpi; 
	size_t *lp, *up; 
	
	lx = L_GainP->getX();
	li = L_GainP->getI();
	lp = L_GainP->getP();
		
	ux = U_GainP->getX();
	ui = U_GainP->getI();
	up = U_GainP->getP();
	
	rp = L_GainP->getRp();
	cpi = L_GainP->getCpi();
	rows = L_GainP->getRows();
	cols = L_GainP->getCols();	// Initialize arrays and variables
	
	uint__t n, nnz_p_L, nnz_p_U, n_e, nnz_e;  // nnz_pp_L, nnz_pp_U,
	int i, j, p; // ret, iter
	//double maxDeltaP=0, maxDeltaQ=0, max_change_ex=maxchange;

	const double pi_value = 3.141592653589793;

	// Get the dimension and the nnz for the LU of B' and B"
	n = L_GainP->getn();	// get number of rows or vertices
	nnz_p_L=L_GainP->getnnz(); nnz_p_U=U_GainP->getnnz();
	//nnz_pp_L=gMatrix_pp_L.size(); nnz_pp_U=gMatrix_pp_U.size();
	//n_e=gVertex_Va.data_.size(); // nnz_e=gMatrix_Ybus.data_.size(); // get the size of the Y bus matrix


	
	printf("\n\n------------------------------------------------------------------------------------------- ");
	std::cout << "Lp Number of rows:" << n << ",\tNumber of nonzeros:" << nnz_p_L << std::endl;
	std::cout << "Up Number of rows:" << n << ",\tNumber of nonzeros:" << nnz_p_U << std::endl;
	printf("-------------------------------------------------------------------------------------------- \n\n");
	
	
	//// For Y-bus
	double *H_r_P, *Va; // *Pn, *Qn; *eG, *eB, *deltaQ, *Vm, 
	//uint *ei, *piv, *btype; 
	//uint *ep; 

	// For L' and U'
	//double *lx, *ux, *rows, *cols; 
	//uint *li, *ui, *rp, *cpi; 
	//size_t *lp, *up; 
	//uint__t fctnnz; // number of nonzeros after the matrix is factorized
	//size_t l_nnz = gLp_edge.data_.size(), u_nnz = gUp_edge.data_.size(); // number of nonzeros in L and U factors 
	
	// lx = (double *)malloc(sizeof(double)*(nnz_p_L));
	// li = (uint *)malloc(sizeof(uint)*(nnz_p_L)); // column indices 
	// lp = (size_t *)malloc(sizeof(size_t)*(n+1)); // initial pointer
	// ux = (double *)malloc(sizeof(double)*(nnz_p_U));
	// ui = (uint *)malloc(sizeof(uint)*(nnz_p_U)); // column indices 
	// up = (size_t *)malloc(sizeof(size_t)*(n+1)); // initial pointer 
	// rows = (double *)malloc(sizeof(double)*(n));
	// cols = (double *)malloc(sizeof(double)*(n));
	// rp = (uint *)malloc(sizeof(uint)*(n));
	// cpi = (uint *)malloc(sizeof(uint)*(n));


	H_r_P = (double *)malloc(sizeof(double)*(n)); // b in the Ax=b
	//deltaQ = (double *)malloc(sizeof(double)*(n)); // b in the Ax=b

	//Vm = (double *)malloc(sizeof(double)*(n)); 
	Va = (double *)malloc(sizeof(double)*(n)); 

	//Pn = (double *)malloc(sizeof(double)*(n)); 
	//Qn = (double *)malloc(sizeof(double)*(n));

	//btype = (uint *)malloc(sizeof(uint)*(n));

	double sum, diag;
	double *b;  // *delta, 
	//delta = (double *)malloc(sizeof(double)*(n));
	b = (double *)malloc(sizeof(double)*(n));

	std::cout << " ======================== Initialization of ararys used to store the factorized matrix and LU ========================"<< std::endl;

	struct timeval t2_st, t2_end, t3_st, t3_end; long seconds, useconds;

	//initialization to 0
	// [tc] use memset here
	//memset(Vm, 0, sizeof(double)*(n));
	memset(Va, 0, sizeof(double)*(n));
	memset(H_r_P, 0, sizeof(double)*(n));
	//memset(deltaQ, 0, sizeof(double)*(n));


	// ========================================= Convert pointers and the Vertex =========================================
	gettimeofday(&t2_st, 0); 
	//int i_p; // i_pp
	int ii = 1; // counter for the array  
	// Assign value to the first position of the pointer array 
	//lp[0] = 0;
	//up[0] = 0;
	//lp_pp[0] = 0;
	//up_pp[0] = 0;    
	//ep[0] = 0;
	for (i=0; i<n; ++i){
	//process pointers for the matrices
	
	//lp[ii] = gGainP_LU_vertex.data_[i].Lp + lp[ii-1];
	//up[ii] = gGainP_LU_vertex.data_[i].Up + up[ii-1];
	
	//lp_pp[ii] = gVertex_pp.data_[i].Lp + lp_pp[ii-1];
	//up_pp[ii] = gVertex_pp.data_[i].Up + up_pp[ii-1];    
	//ep[ii] = gVertex_Ybus.data_[i].ep + ep[ii-1];
	//ii++;
	
	Va[i] = gVa.data_[i];
	H_r_P[i] = gH_r_P.data_[i];
	  
	// debug print


	//std::cout<< i << ",deltaP,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  

	//process factorized information for B'
	//i_p = gVertex_p.data_[i].exId - 1;
	//rows[i]= gGainP_LU_vertex.data_[i].row_scaling;
	//cols[i]= gGainP_LU_vertex.data_[i].col_scaling;
	//rp[i] = gGainP_LU_vertex.data_[i].rp;
	//cpi[i] = gGainP_LU_vertex.data_[i].cpi;  

	}

	// ========================================= Convert L' and U' =========================================
	//int i_l = 0, i_u = 0;

	// for (ii=0; ii<gGainP_L.data_.size(); ++ii) {
	// li[ii] = gGainP_L.data_[ii].cpi;
	// lx[ii] = gGainP_L.data_[ii].value1; 
	// //std::cout<< i << ",li,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  	
	// } 

	// for (ii=0; ii<gGainP_U.data_.size(); ++ii) {  
	// ui[ii] = gGainP_U.data_[ii].cpi;
	// ux[ii] = gGainP_U.data_[ii].value1; 
	// //std::cout<< i << ",deltaP,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  	
	// } 


	// Done converting all input HeapAccums to arrays

	// Get the time to convert the data between containers and arrays
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;

	printf("\n\n============================================================================================== ");
	std::cout << "Time to convert data to GRAPHLU time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");		

	//continue to do the iterations only if there is no error from factorization
	printf("\n\n============================================================================================== ");
	std::cout << "\n Start iteratively updating deltaP and deltaQ " << std::endl;
	printf("\n================================================================================================\n ");

	gettimeofday(&t2_st, 0);

 		// solve for deltaVa
 		// std::cout << " ======================== Solve for V angle ========================"<<std::endl;
     
    // A*x = b
    // multiply deltaP with row scaling (rows) to get b
   	for (i=0; i<n; ++i){
      b[i]=H_r_P[rp[i]]*rows[rp[i]];
   	}
    
   	// Full forward substitution 
   	for (i=0; i<n; ++i)
   	{
   		sum=0.;
   		diag=0.;
   		for(p=lp[i]; p< lp[i+1]; ++p){
   			if (i!=li[p]){ // if it is not diagnal element
   				j=li[p];
   				sum += lx[p] * b[j];
   			}
   			else{
   				diag = lx[p];
   			}
   		}
   		b[i] -= sum; 
   		b[i] /= diag; 
   	} 
    
   	// Full backward substitution 
   	for (i=n-1; i>=0; --i)
   	{
   		sum=0.;
   		for(p=up[i]; p< up[i+1]; ++p){
   			if (i!=ui[p]){ // if it is not diagnal element
   				sum += ux[p] * b[ui[p]]; 
   			}
   		}
   
   		b[i] -= sum; 
		
		//if(gmax_dVa < abs(b[i]*cols[i]))
		//{
		//	gmax_dVa = abs(b[i]*cols[i]);
		//}
		//std::cout << "gmax_dVa: " << gmax_dVa << std::endl;
		
   	}      
   	// permute the array back to the original order
   	//for (i=0; i<n; ++i)
     //{
   	//	delta[i]=b[cpi[i]];
     //  //std::cout << "delta: " << i << ", " << delta[i] << std::endl;           
 		//}
       
    // Update V angle (Va)
 		for (i=0; i<n; ++i)
		{ 
      //Va[i] -= delta[i];  
			Va[i] += b[cpi[i]]*cols[cpi[i]];    // x(k+1) = x(k) + delta_x 
	  	
			if(gmax_dVa < abs(b[cpi[i]]*cols[cpi[i]]))
			{
				gmax_dVa = abs(b[cpi[i]]*cols[cpi[i]]);
			}
			//std::cout << "gmax_dVa: " << gmax_dVa << std::endl;
	  
	  
		}
    
 	gettimeofday(&t2_end, 0);
 	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
 	useconds = t2_end.tv_usec  - t2_st.tv_usec;
 	
 	printf("\n\n============================================================================================== ");
 	std::cout << "Time for solving V angle: " << (seconds*1000000 + useconds) << " us" << std::endl;
 	printf("\n================================================================================================ ");
 
  string result = "Finished Va Forward/Backward Substitution";
 	// ----------------------------------------------------------------------------------------------------
 	// 								Store the Result back to the array
 	// ----------------------------------------------------------------------------------------------------
 	// Get solved Vm and Va back into HeapAccum
 	for (ii=0; ii<n; ++ii) {
       gVa.data_[ii] = Va[ii];  // in radian
       //gVertex_Vm.data_[ii].Vm = Vm[ii];
       //std::cout << "Final, " << ii+1 << ", Vm, " << Vm[ii]<< ", Va" << ii+1 << ", " << Va[ii]/pi_value*180<< std::endl;
 	}
	gettimeofday(&t1_end, 0);
	seconds=t1_end.tv_sec  - t1_st.tv_sec; 
	useconds = t1_end.tv_usec  - t1_st.tv_usec;

	printf("\n\n============================================================================================== ");
	std::cout << "Time to solve Gain_P_static time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");		

  
EXIT:

	free(H_r_P);
	//free(deltaQ);
	//free(Vm);
	free(Va);
	free(b);
	//free(Pn);
	//free(Qn);	
	//free(L_GainP);
	//free(U_GainP);

	//free(lx);
	//free(li);
	//free(lp);
	//free(ux);
	//free(ui);
	//free(up); 
	//free(rows);
	//free(cols);
	//free(rp);
	//free(cpi); 
	



	printf("\n\n----------------------------------------------------------------------------------------");
	printf("\t\t End of Running SE_solve_GainP_static C function!");
	printf("\t\t ----------------------------------------------------------------------------------------\n\n");	
 
  return result;
}

// **********************************************************************************************************************************
// Created by: Qiwei
//2019 summer intern
// LU decoposition to solve K matrix
inline string SE_solve_BDIG_static (ArrayAccum<SumAccum<double>>& gK_m, ArrayAccum<SumAccum<double>>& gH_R_m, double& slackbus) {

	CSRMatrix *L_GainP = SingletonMatrixInterface::getMatrixGainP_L();
	CSRMatrix *U_GainP = SingletonMatrixInterface::getMatrixGainP_U();
	int slack_bus = (int)slackbus;
    std::cout << " slack_bus " << slack_bus << std::endl;
	// For L' and U'
	double *lx, *ux, *rows, *cols; 
	uint *li, *ui, *rp, *cpi; 
	size_t *lp, *up; 
	int i, j, p; // ret, iter
	lx = L_GainP->getX();
	li = L_GainP->getI();
	lp = L_GainP->getP();
		
	ux = U_GainP->getX();
	ui = U_GainP->getI();
	up = U_GainP->getP();
	
	rp = L_GainP->getRp();
	cpi = L_GainP->getCpi();
	rows = L_GainP->getRows();
	cols = L_GainP->getCols();	// Initialize arrays and variables
	
	uint__t n, nnz_p_L, nnz_p_U, n_e, nnz_e;  // nnz_pp_L, nnz_pp_U,
	//double maxDeltaP=0, maxDeltaQ=0, max_change_ex=maxchange;

	const double pi_value = 3.141592653589793;

	// Get the dimension and the nnz for the LU of B' and B"
	n = L_GainP->getn();	// get number of rows or vertices
	nnz_p_L=L_GainP->getnnz(); nnz_p_U=U_GainP->getnnz();
	
	//// For Y-bus
	double *H_R_m, *K_m; // *Pn, *Qn; *eG, *eB, *deltaQ, *Vm, 

	H_R_m = (double *)malloc(sizeof(double)*(n+1)); // b in the Ax=b

	K_m = (double *)malloc(sizeof(double)*(n+1)); 

	double sum, diag;
	double *b;  // *delta, 
	//delta = (double *)malloc(sizeof(double)*(n));
	b = (double *)malloc(sizeof(double)*(n));

	std::cout << " ======================== Initialization of ararys used to store the factorized matrix and LU ========================"<< std::endl;

	memset(K_m, 0, sizeof(double)*(n));
	memset(H_R_m, 0, sizeof(double)*(n));


	// ========================================= Convert pointers and the Vertex =========================================
	int ii = 1; // counter for the array  
	for (i=0; i<n+1; ++i){
		K_m[i] = gK_m.data_[i];
		H_R_m[i] = gH_R_m.data_[i];
	}
	for (i = slack_bus-1; i<n;i++)
	{
		K_m[i] = K_m[i+1];
	    H_R_m[i] = H_R_m[i+1];
	}
    // A*x = b
    // multiply deltaP with row scaling (rows) to get b
   	for (i=0; i<n; ++i){
      b[i]=H_R_m[rp[i]]*rows[rp[i]];
   	}
    
   	// Full forward substitution 
   	for (i=0; i<n; ++i)
   	{
   		sum=0.;
   		diag=0.;
   		for(p=lp[i]; p< lp[i+1]; ++p){
   			if (i!=li[p]){ // if it is not diagnal element
   				j=li[p];
   				sum += lx[p] * b[j];
   			}
   			else{
   				diag = lx[p];
   			}
   		}
   		b[i] -= sum; 
   		b[i] /= diag; 
   	} 
    
   	// Full backward substitution 
   	for (i=n-1; i>=0; --i)
   	{
   		sum=0.;
   		for(p=up[i]; p< up[i+1]; ++p){
   			if (i!=ui[p]){ // if it is not diagnal element
   				sum += ux[p] * b[ui[p]]; 
   			}
   		}
   
   		b[i] -= sum; 
   	}      
    // Update V angle (Va)
 		for (i=0; i<n; ++i)
		{ 
			K_m[i] += b[cpi[i]]*cols[cpi[i]];    // x(k+1) = x(k) + delta_x 
		}
    
 
  string result = "Finished K_m Forward/Backward Substitution";
 	// ----------------------------------------------------------------------------------------------------
 	// 								Store the Result back to the array
 	// ----------------------------------------------------------------------------------------------------
 	// Get solved Vm and Va back into HeapAccum
 	for (ii=0; ii<n; ++ii) {
       gK_m.data_[ii] = K_m[ii];  // in radian
 	}
  
EXIT:

	free(K_m);
	free(H_R_m);
	free(b);



	printf("\n\n----------------------------------------------------------------------------------------");
	printf("\t\t End of Running SE_solve_GainP_static C function!");
	printf("\t\t ----------------------------------------------------------------------------------------\n\n");	
 
  return result;
}

// **********************************************************************************************************************************
 // Created by: Chen Yuan, chen.yuan@geirina.net
 // Date: 11/15/2017
 // This code performs to solve state estimation GainQ * deltaVm = H_Q' * r_Q. This version does NOT perform LU factorization. GainQ is assumed to be previously factorized and
 // the resulting L and U matrices are taken as inputs. The equation is solved with forward/backward substitution.
 // History: 
 // 08/25/2017: [Chen Yuan] GainP matrix has been previously factorized. Built this function to solve GainQ * deltaVa = H_Q'*r_Q, by employing forward/backward substitution, and update voltage angle, Vm. 
 
//template <//typename T_deltaP_vertex, typename deltaP_vertex_comp,
          //typename T_vertex_Va,  
          //typename GainQ_L_matrix, typename GainQ_L_matrix_comp,  
          //typename GainQ_U_matrix, typename GainQ_U_matrix_comp,
          //typename GainQ_LU_vertex, typename GainQ_LU_vertex_comp>
//inline string SE_solve_GainQ (ArrayAccum<SumAccum<double>>& gVm, ArrayAccum<SumAccum<double>>& gH_r_Q,
              //HeapAccum<GainQ_L_matrix, GainQ_L_matrix_comp>& gGainQ_L, HeapAccum<GainQ_U_matrix, GainQ_U_matrix_comp>& gGainQ_U, HeapAccum<GainQ_LU_vertex, GainQ_LU_vertex_comp>& gGainQ_LU_vertex, double& gmax_dVm) {
inline string SE_solve_GainQ_static (ArrayAccum<SumAccum<double>>& gVm, ArrayAccum<SumAccum<double>>& gH_r_Q, double& gmax_dVm) {

	struct timeval t1_st, t1_end;
	gettimeofday(&t1_st, 0);
	// extern size_t g_si, g_sd, g_sp;
	printf("\n\n------------------------------------------------------------------------------------------- ");
	printf("\n Start Running SE_solve_GainQ_static function!\n");
	printf("-------------------------------------------------------------------------------------------- \n\n");	
	// ------------------------------------------------------------------------------------------------
	// 				Initialize variables and arrays
	// ------------------------------------------------------------------------------------------------
	
	CSRMatrix *L_GainQ = SingletonMatrixInterface::getMatrixGainQ_L();
	CSRMatrix *U_GainQ = SingletonMatrixInterface::getMatrixGainQ_U();
	
	
	// For L' and U'
	double *lx, *ux, *rows, *cols; 
	uint *li, *ui, *rp, *cpi; 
	size_t *lp, *up; 
	
	
	lx = L_GainQ->getX();
	li = L_GainQ->getI();
	lp = L_GainQ->getP();
		
	ux = U_GainQ->getX();
	ui = U_GainQ->getI();
	up = U_GainQ->getP();
	
	rp = L_GainQ->getRp();
	cpi = L_GainQ->getCpi();
	rows = L_GainQ->getRows();
	cols = L_GainQ->getCols();	// Initialize arrays and variables
	
	


	// Initialize arrays and variables
	uint__t n, nnz_p_L, nnz_p_U, n_e, nnz_e;  // nnz_pp_L, nnz_pp_U,
	int i, j, p; // ret, iter
	//double maxDeltaP=0, maxDeltaQ=0, max_change_ex=maxchange;

	const double pi_value = 3.141592653589793;

	// Get the dimension and the nnz for the LU of B' and B"
	n=L_GainQ->getn();	
	nnz_p_L=L_GainQ->getnnz(); nnz_p_U=U_GainQ->getnnz();
	//nnz_pp_L=gMatrix_pp_L.size(); nnz_pp_U=gMatrix_pp_U.size();
	//n_e=gVertex_Va.data_.size(); // nnz_e=gMatrix_Ybus.data_.size(); // get the size of the Y bus matrix

	
	printf("\n\n------------------------------------------------------------------------------------------- ");
	std::cout << "Lp Number of rows:" << n << ",\tNumber of nonzeros:" << nnz_p_L << std::endl;
	std::cout << "Up Number of rows:" << n << ",\tNumber of nonzeros:" << nnz_p_U << std::endl;
	
	printf("-------------------------------------------------------------------------------------------- \n\n");	

	//// For Y-bus
	double *H_r_Q, *Vm; // *Pn, *Qn; *eG, *eB, *deltaQ, *Vm, 
	//uint *ei, *piv, *btype; 
	//uint *ep; 

	// For L' and U'
	//double *lx, *ux, *rows, *cols; 
	//uint *li, *ui, *rp, *cpi; 
	//size_t *lp, *up; 
	//uint__t fctnnz; // number of nonzeros after the matrix is factorized
	//size_t l_nnz = gLp_edge.data_.size(), u_nnz = gUp_edge.data_.size(); // number of nonzeros in L and U factors 

	// lx = (double *)malloc(sizeof(double)*(nnz_p_L));
	// li = (uint *)malloc(sizeof(uint)*(nnz_p_L)); // column indices 
	// lp = (size_t *)malloc(sizeof(size_t)*(n+1)); // initial pointer
	// ux = (double *)malloc(sizeof(double)*(nnz_p_U));
	// ui = (uint *)malloc(sizeof(uint)*(nnz_p_U)); // column indices 
	// up = (size_t *)malloc(sizeof(size_t)*(n+1)); // initial pointer 
	// rows = (double *)malloc(sizeof(double)*(n));
	// cols = (double *)malloc(sizeof(double)*(n));
	// rp = (uint *)malloc(sizeof(uint)*(n));
	// cpi = (uint *)malloc(sizeof(uint)*(n));


	H_r_Q = (double *)malloc(sizeof(double)*(n)); // b in the Ax=b
	//deltaQ = (double *)malloc(sizeof(double)*(n)); // b in the Ax=b

	//Vm = (double *)malloc(sizeof(double)*(n)); 
	Vm = (double *)malloc(sizeof(double)*(n)); 

	//Pn = (double *)malloc(sizeof(double)*(n)); 
	//Qn = (double *)malloc(sizeof(double)*(n));

	//btype = (uint *)malloc(sizeof(uint)*(n));

	double sum, diag;
	double *b;  // *delta, 
	//delta = (double *)malloc(sizeof(double)*(n));
	b = (double *)malloc(sizeof(double)*(n));

	std::cout << " ======================== Initialization of ararys used to store the factorized matrix and LU ========================"<< std::endl;

	struct timeval t2_st, t2_end, t3_st, t3_end; long seconds, useconds;

	//initialization to 0
	// [tc] use memset here
	//memset(Vm, 0, sizeof(double)*(n));
	memset(Vm, 0, sizeof(double)*(n));
	memset(H_r_Q, 0, sizeof(double)*(n));
	//memset(deltaQ, 0, sizeof(double)*(n));


	// ========================================= Convert pointers and the Vertex =========================================
	gettimeofday(&t2_st, 0); 
	//int i_p; // i_pp
	int ii = 1; // counter for the array  
	// Assign value to the first position of the pointer array 
	//lp[0] = 0;
	//up[0] = 0;
	//lp_pp[0] = 0;
	//up_pp[0] = 0;    
	//ep[0] = 0;
	for (i=0; i<n; ++i){
	//process pointers for the matrices
	//lp[ii] = gGainQ_LU_vertex.data_[i].Lp + lp[ii-1];
	//up[ii] = gGainQ_LU_vertex.data_[i].Up + up[ii-1];
	//lp_pp[ii] = gVertex_pp.data_[i].Lp + lp_pp[ii-1];
	//up_pp[ii] = gVertex_pp.data_[i].Up + up_pp[ii-1];    
	//ep[ii] = gVertex_Ybus.data_[i].ep + ep[ii-1];
	//ii++;
	
	Vm[i] = gVm.data_[i];
	H_r_Q[i] = gH_r_Q.data_[i];
	  
	// debug print


	//std::cout<< i << ",deltaP,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  

	//process factorized information for B'
	//i_p = gVertex_p.data_[i].exId - 1;
	//rows[i]= gGainQ_LU_vertex.data_[i].row_scaling;
	//cols[i]= gGainQ_LU_vertex.data_[i].col_scaling;
	//rp[i] = gGainQ_LU_vertex.data_[i].rp;
	//cpi[i] = gGainQ_LU_vertex.data_[i].cpi;  

	}

	// ========================================= Convert L' and U' =========================================
	//int i_l = 0, i_u = 0;

	//for (ii=0; ii<gGainQ_L.data_.size(); ++ii) {
	//li[ii] = gGainQ_L.data_[ii].cpi;
	//lx[ii] = gGainQ_L.data_[ii].value1; 
	////std::cout<< i << ",li,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  	
	//} 

	//for (ii=0; ii<gGainQ_U.data_.size(); ++ii) {  
	//ui[ii] = gGainQ_U.data_[ii].cpi;
	//ux[ii] = gGainQ_U.data_[ii].value1; 
	////std::cout<< i << ",deltaP,"<< deltaP[i] << ",deltaQ,"<< deltaQ[i] << std::endl;  	
	//} 


	// Done converting all input HeapAccums to arrays

	// Get the time to convert the data between containers and arrays
	gettimeofday(&t2_end, 0);
	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
	useconds = t2_end.tv_usec  - t2_st.tv_usec;

	printf("\n\n============================================================================================== ");
	std::cout << "Time to convert data to GRAPHLU time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");		

	//continue to do the iterations only if there is no error from factorization
	printf("\n\n============================================================================================== ");
	std::cout << "\n Start iteratively updating deltaP and deltaQ " << std::endl;
	printf("\n================================================================================================\n ");

	gettimeofday(&t2_st, 0);

 		// solve for deltaVa
 		// std::cout << " ======================== Solve for V angle ========================"<<std::endl;
     
    // A*x = b
    // multiply deltaP with row scaling (rows) to get b
   	for (i=0; i<n; ++i){
      b[i]=H_r_Q[rp[i]]*rows[rp[i]];
   	}
    
   	// Full forward substitution 
   	for (i=0; i<n; ++i)
   	{
   		sum=0.;
   		diag=0.;
   		for(p=lp[i]; p< lp[i+1]; ++p){
   			if (i!=li[p]){ // if it is not diagnal element
   				j=li[p];
   				sum += lx[p] * b[j];
   			}
   			else{
   				diag = lx[p];
   			}
   		}
   		b[i] -= sum; 
   		b[i] /= diag; 
   	} 
    
   	// Full backward substitution 
   	for (i=n-1; i>=0; --i)
   	{
   		sum=0.;
   		for(p=up[i]; p< up[i+1]; ++p){
   			if (i!=ui[p]){ // if it is not diagnal element
   				sum += ux[p] * b[ui[p]]; 
   			}
   		}
   
   		b[i] -= sum; 
		
		if(gmax_dVm < abs(b[i]*cols[i]))
		{
			gmax_dVm = abs(b[i]*cols[i]);
		}
		//std::cout << "gmax_dVm: " << gmax_dVm << std::endl;
		
   	}      
   	// permute the array back to the original order
   	//for (i=0; i<n; ++i)
     //{
   	//	delta[i]=b[cpi[i]];
     //  //std::cout << "delta: " << i << ", " << delta[i] << std::endl;           
 		//}
       
    // Update V angle (Va)
 		for (i=0; i<n; ++i)
		{ 
		  //Va[i] -= delta[i];  
		  Vm[i] += b[cpi[i]]*cols[cpi[i]];    // x(k+1) = x(k) + delta_x
		  if(gmax_dVm < abs(b[cpi[i]]*cols[cpi[i]]))
			{
				gmax_dVm = abs(b[cpi[i]]*cols[cpi[i]]);
			}
			//std::cout << "gmax_dVm: " << gmax_dVm << std::endl;
		  
		  
		}
    
 	gettimeofday(&t2_end, 0);
 	seconds=t2_end.tv_sec  - t2_st.tv_sec; 
 	useconds = t2_end.tv_usec  - t2_st.tv_usec;
 	
 	printf("\n\n============================================================================================== ");
 	std::cout << "Time for solving V angle: " << (seconds*1000000 + useconds) << " us" << std::endl;
 	printf("\n================================================================================================ ");
 
  string result = "Finished Vm Forward/Backward Substitution";
 	// ----------------------------------------------------------------------------------------------------
 	// 								Store the Result back to the array
 	// ----------------------------------------------------------------------------------------------------
 	// Get solved Vm and Va back into HeapAccum
 	for (ii=0; ii<n; ++ii) {
       gVm.data_[ii] = Vm[ii];  // in radian
       //gVertex_Vm.data_[ii].Vm = Vm[ii];
       //std::cout << "Final, " << ii+1 << ", Vm, " << Vm[ii]<< ", Va" << ii+1 << ", " << Va[ii]/pi_value*180<< std::endl;
 	}

	gettimeofday(&t1_end, 0);
	seconds=t1_end.tv_sec  - t1_st.tv_sec; 
	useconds = t1_end.tv_usec  - t1_st.tv_usec;

	printf("\n\n============================================================================================== ");
	std::cout << "Time to solve Gain_Q static time:: " << (seconds*1000000 + useconds) << " us" << std::endl;
	printf("\n================================================================================================ ");		
  
EXIT:

	free(H_r_Q);
	//free(deltaQ);
	//free(Vm);
	free(Vm);
	free(b);
	//free(Pn);
	//free(Qn);	
	
	//free(L_GainQ);
	//free(U_GainQ);

	//free(lx);
	//free(li);
	//free(lp);
	//free(ux);
	//free(ui);
	//free(up); 
	//free(rows);
	//free(cols);
	//free(rp);
	//free(cpi);	


	printf("\n\n----------------------------------------------------------------------------------------");
	printf("\t\t End of Running SE_solve_GainQ_static C function!");
	printf("\t\t ----------------------------------------------------------------------------------------\n\n");	
 
  return result;
}






// **********************************************************************************************************************************
// Created by: Chen Yuan, chen.yuan@geirina.net
// Date: 08/24/2017
// This function takes an index and an ArrayAccum<SumAccum<>> as inputs and returns Va (voltage magnitude) as output. Specified for gVa ArrayAccum<SumAccum<>>
// History:
// **********************************************************************************************************************************

inline double SE_getVa (uint64_t& index, const ArrayAccum<SumAccum<double>> &gVa) {
  return gVa.data_[index];
} 

// **********************************************************************************************************************************
// Created by: Chen Yuan, chen.yuan@geirina.net
// Date: 08/25/2017
// This function takes an index and an ArrayAccum<SumAccum<>> as inputs and returns Vm (voltage magnitude) as output. Specified for gVm ArrayAccum<SumAccum<>>
// History:
// **********************************************************************************************************************************

inline double SE_getVm (uint64_t& index, const ArrayAccum<SumAccum<double>> &gVm) {
  return gVm.data_[index];
} 
		
#endif /* EXPRFUNCTIONS_HPP_ */

