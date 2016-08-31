/*
 * ungapped_extension.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "ungapped_extension.h"
#include "energy_par.h"
#include "intloops.h"
#include <iostream>


void UngappedExtension::Run(vector<Hit> &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility){
  for(int x = 0; x < candidate.size(); x++){
    double min_energy = candidate[x].GetEnergy();
    double energy = candidate[x].GetEnergy();
    int i = candidate[x].GetQSp(); int p = candidate[x].GetQSp();
    int j = candidate[x].GetDbSp(); int q = candidate[x].GetDbSp();
    int min_p = p; int min_q = q;
    int db_seq_id = candidate[x].GetDbSeqId();
    int db_seq_id_start = candidate[x].GetDbSeqIdStart(); int db_seq_id_end = db_seq_id_start + candidate[x].GetDbLength()-1;
    int min_db_seq_id_start = db_seq_id_start;
    while(true){
      i--; j--;
      db_seq_id_end++;
      if(i < 0 || j < 0 || query_seq[i] < 2 || db_seq[j] < 2){
	break;
      }
      energy += query_accessibility[i] - query_accessibility[i+1] + query_conditional_accessibility[i+_min_accessible_length];
      energy += db_conditional_accessibility[db_seq_id][db_seq_id_end];
      int q_char = query_seq[i] <=5 ? query_seq[i]-1 : query_seq[i]-5;
      int db_char = db_seq[j] <=5 ? db_seq[j]-1 : db_seq[j]-5;
      int type = BP_pair[q_char][db_char];
      if(type!=0){
	int q_char2 = query_seq[p] <=5 ? query_seq[p]-1 : query_seq[p]-5;
	int db_char2 = db_seq[q] <=5 ? db_seq[q]-1 : db_seq[q]-5;
	int type2 = BP_pair[q_char2][db_char2];
	type2 = rtype[type2];
	energy += LoopEnergy(type, type2, i, j, p, q, query_seq, db_seq);
	if(energy < min_energy){
	  min_energy = energy;
	  min_p = i;
	  min_q = j;
	}
	
	p = i;
	q = j;
      }
      
      if(min_p-i >= _drop_out_score){
	break;
      }
    }

    energy = min_energy;
    int k = candidate[x].GetQSp() + candidate[x].GetQLength()-1; int r = candidate[x].GetQSp() + candidate[x].GetQLength()-1;
    int l = candidate[x].GetDbSp() + candidate[x].GetQLength()-1; int s = candidate[x].GetDbSp() + candidate[x].GetQLength()-1;
    int min_r = r; int min_s = s;
    while(true){
      k++; l++;
      db_seq_id_start--;
      if(query_seq[k] < 2 || db_seq[l] < 2){
	break;
      }
      energy += query_conditional_accessibility[k];
      energy += db_accessibility[db_seq_id][db_seq_id_start] - db_accessibility[db_seq_id][db_seq_id_start+1] + db_conditional_accessibility[db_seq_id][db_seq_id_start+_min_accessible_length];
      int q_char2 = query_seq[k] <=5 ? query_seq[k]-1 : query_seq[k]-5;
      int db_char2 = db_seq[l] <=5 ? db_seq[l]-1 : db_seq[l]-5;
      int type2 = BP_pair[q_char2][db_char2];
      type2 = rtype[type2];
      if(type2!=0){
	int q_char = query_seq[r] <=5 ? query_seq[r]-1 : query_seq[r]-5;
	int db_char = db_seq[s] <=5 ? db_seq[s]-1 : db_seq[s]-5;
	int type = BP_pair[q_char][db_char];
	energy += LoopEnergy(type, type2, r, s, k, l, query_seq, db_seq);
	if(energy < min_energy){
	  min_energy = energy;
	  min_r = k;
	  min_s = l;
	  min_db_seq_id_start = db_seq_id_start;
	}
	r = k;
	s = l;
      }
      if(k-min_r >= _drop_out_score){
	break;
      }
    }
    candidate[x].SetDbSeqIdStart(min_db_seq_id_start);
    candidate[x].SetQSp(min_p);
    candidate[x].SetDbSp(min_q);
    candidate[x].SetQLength(min_r-min_p+1);
    candidate[x].SetDbLength(min_r-min_p+1);
    candidate[x].SetEnergy(min_energy);
  }
}

double UngappedExtension::LoopEnergy(int type, int type2,int i,int j,int p,int q, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq){
  double z = 0;
  int u1 = p-i-1;
  int u2 = q-j-1;
  
  if ((u1==0) && (u2==0)){
    z = stack37[type][type2];
  }else{
    int a = query_seq[i+1] <=5 ? query_seq[i+1]-1 : query_seq[i+1]-5;
    int b = db_seq[j+1] <=5 ? db_seq[j+1]-1 : db_seq[j+1]-5;
    int c = query_seq[p-1] <=5 ? query_seq[p-1]-1 : query_seq[p-1]-5;
    int d = db_seq[q-1] <=5 ? db_seq[q-1]-1 : db_seq[q-1]-5;
    if (u1+u2==2) {
      z = int11_37[type][type2][a][b];
    }else if ((u1==1) && (u2==2)){
      z = int21_37[type][type2][a][d][b];
    }else if ((u1==2) && (u2==1)){
      z = int21_37[type2][type][d][a][c];
    }else if ((u1==2) && (u2==2)){
    z = int22_37[type][type2][a][c][d][b];
    }else{
      z = internal_loop37[u1+u2]+mismatchI37[type][a][b]+mismatchI37[type2][d][c];
    }
  }
  return double(z)/100;
}
