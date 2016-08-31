#include "raccess.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>


void Raccess::Run(string &sequence){
  clock_t start, end;
  Initiallize(sequence);
  CalcInsideVariable();
  CalcOutsideVariable();
  start = clock();
  CalcAccessibility();
  end = clock();
  Clear();
}

void Raccess::Run(string &sequence, vector<float> &accessibility, vector<float> &conditional_accessibility){
  Initiallize(sequence);
  CalcInsideVariable();
  CalcOutsideVariable();
  CalcAccessibility(accessibility, conditional_accessibility);
  Clear();
}

void Raccess::Initiallize(string &sequence){
  _seq_length = sequence.length();
  _int_sequence.resize(_seq_length+1);
  for(int i = 0;i < _seq_length;i++){
    if(sequence[i] == 'A' || sequence[i] == 'a'){
      _int_sequence[i+1] = 1;
    }else if(sequence[i] == 'C' || sequence[i] == 'c'){
      _int_sequence[i+1] = 2;
    }else if(sequence[i] == 'G' || sequence[i] == 'g'){
      _int_sequence[i+1] = 3;
    }else if(sequence[i] == 'T' || sequence[i] == 't' || sequence[i] == 'U' || sequence[i] == 'u'){
      _int_sequence[i+1] = 4;
    }else{
      _int_sequence[i+1] = 0;
    }
  }
  
  _Alpha_outer.resize(_seq_length+1, 0);
  _Alpha_stem.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_stemend.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multibif.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi1.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Alpha_multi2.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  

  _Beta_outer.resize(_seq_length+1, 0);
  _Beta_stem.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_stemend.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multibif.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi1.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
  _Beta_multi2.resize(_seq_length+1,vector<double>(_maximal_span+2, -INF));
}

void Raccess::CalcInsideVariable(){
  for (int j =TURN+1; j <= _seq_length; j++){
    for (int i=j-TURN; i >= max(0,j-_maximal_span-1); i--){
      //Alpha_stem
      int type = BP_pair[_int_sequence[i+1]][_int_sequence[j]];
      int type2 = BP_pair[_int_sequence[i+2]][_int_sequence[j-1]];
      
      double temp = 0; bool flag = 0;
      if (type != 0) {
	type2 = rtype[type2];
	if(_Alpha_stem[i+1][j-i-2] != -INF){
	  //Stem¨Stem
	  if(type2 != 0){
	    temp = _Alpha_stem[i+1][j-i-2]+ LoopEnergy(type, type2,i+1,j,i+2,j-1);
	  }
	  flag = 1;
	}
    
	if(_Alpha_stemend[i+1][j-i-2] != -INF){
	  //Stem¨StemEnd
	  temp = flag == 1 ? logsumexp(temp,_Alpha_stemend[i+1][j-i-2]) : _Alpha_stemend[i+1][j-i-2];
	  flag = 1;
	}

	_Alpha_stem[i][j-i] = flag == 0 ? -INF : temp;
      }else{
	_Alpha_stem[i][j-i] = -INF;
      }
      
      //Alpha_multiBif
      temp = 0; flag = 0;
      for (int k=i+1; k<=j-1; k++){
	if(_Alpha_multi1[i][k-i] != -INF && _Alpha_multi2[k][j-k] != -INF){
	  temp = flag == 0 ? _Alpha_multi1[i][k-i]+_Alpha_multi2[k][j-k] : logsumexp(temp,_Alpha_multi1[i][k-i]+_Alpha_multi2[k][j-k]);
	  flag = 1;
	}
      }
      _Alpha_multibif[i][j-i] = flag == 0 ? -INF : temp;
      
      //Alpha_multi2
      temp = 0; flag = 0; 
      if (type != 0) {
	if(_Alpha_stem[i][j-i] != -INF){
	  temp = _Alpha_stem[i][j-i]+MLintern+CalcDangleEnergy(type,i,j);
	  flag = 1;
	}
      }
      if(_Alpha_multi2[i][j-i-1] != -INF){
	_Alpha_multi2[i][j-i] = _Alpha_multi2[i][j-i-1]+MLbase;
	if(flag == 1){
	  _Alpha_multi2[i][j-i] = logsumexp(temp,_Alpha_multi2[i][j-i]);
	}
      }else{
	_Alpha_multi2[i][j-i] = flag == 0 ? -INF : temp;
      }
      
      //Alpha_multi1
      if(_Alpha_multi2[i][j-i] != -INF && _Alpha_multibif[i][j-i] != -INF){
	_Alpha_multi1[i][j-i] = logsumexp(_Alpha_multi2[i][j-i],_Alpha_multibif[i][j-i]);
      }else if(_Alpha_multi2[i][j-i] == -INF){
	_Alpha_multi1[i][j-i] = _Alpha_multibif[i][j-i];
      }else if(_Alpha_multibif[i][j-i] == -INF){
	_Alpha_multi1[i][j-i] = _Alpha_multi2[i][j-i];
      }else{
	_Alpha_multi1[i][j-i] = -INF;
      }
      
      //Alpha_multi
      flag = 0;
      if(_Alpha_multi[i+1][j-i-1] != -INF){
	_Alpha_multi[i][j-i] = _Alpha_multi[i+1][j-i-1]+MLbase;
	flag = 1;
      }
      
      if(flag == 1){
	if(_Alpha_multibif[i][j-i] != -INF){
	  _Alpha_multi[i][j-i] = logsumexp(_Alpha_multi[i][j-i],_Alpha_multibif[i][j-i]);
	}
      }else{
	_Alpha_multi[i][j-i] = _Alpha_multibif[i][j-i];
      }
      
      //Alpha_stemend
      if(j != _seq_length){
	temp = 0;
	type = BP_pair[_int_sequence[i]][_int_sequence[j+1]];
	if (type!=0) {
	  //StemEnd¨sn
	  temp = HairpinEnergy(type, i,j+1);
	  
	  //StemEnd¨sm_Stem_sn
	  for (int p =i; p <= min(i+MAXLOOP,j-TURN-2); p++) {
	    int u1 = p-i;
	    for (int q=max(p+TURN+2,j-MAXLOOP+u1); q<=j; q++) {
	      type2 = BP_pair[_int_sequence[p+1]][_int_sequence[q]];
	      if(_Alpha_stem[p][q-p] != -INF){
		if (type2 != 0 && !(p == i && q == j)) {
		  type2 = rtype[type2];
		  temp = logsumexp(temp,_Alpha_stem[p][q-p]+LoopEnergy(type, type2,i,j+1,p+1,q)); 
		}
	      }
	    }
	  }
	  
	  //StemEnd¨Multi
	  int tt = rtype[type];
	  temp = logsumexp(temp,_Alpha_multi[i][j-i]+MLclosing+MLintern+dangle3[tt][_int_sequence[i+1]]+dangle5[tt][_int_sequence[j]]);
	  _Alpha_stemend[i][j-i] = temp;
	}else{
	  _Alpha_stemend[i][j-i] = -INF;
	}
      }
    }
  }
  
  //Alpha_Outer
  for(int i = 1;i <= _seq_length;i++){
    double temp = _Alpha_outer[i-1];
    for(int p = max(0,i-_maximal_span-1); p <i;p++){
      if(_Alpha_stem[p][i-p] != -INF){
	int type = BP_pair[_int_sequence[p+1]][_int_sequence[i]];
	double ao = _Alpha_stem[p][i-p]+CalcDangleEnergy(type,p,i);
	temp = logsumexp(temp,ao+_Alpha_outer[p]);
      }
    }
    _Alpha_outer[i] = temp;
  }
}

double Raccess::CalcDangleEnergy(int type, int a, int b){
  double x = 0;
  if (type != 0) {
    if (a>0) x += dangle5[type][_int_sequence[a]];
    if (b<_seq_length) x += dangle3[type][_int_sequence[b+1]];
    if( b == _seq_length && type>2){
      x += TermAU;
    }
  }
  return(x);
}

void Raccess::CalcOutsideVariable(){
  //Beta_outer
  for(int i = _seq_length-1;i >= 0;i--){
    double temp = _Beta_outer[i+1];
    for(int p = i+1; p <= min(i+_maximal_span+1,_seq_length);p++){
      if(_Alpha_stem[i][p-i] != -INF){
	int type = BP_pair[_int_sequence[i+1]][_int_sequence[p]];
	double bo = _Alpha_stem[i][p-i] + CalcDangleEnergy(type,i,p);
	temp = logsumexp(temp,bo+_Beta_outer[p]);
      }
    }
    _Beta_outer[i] = temp;	
  }
  
  for (int q=_seq_length; q>=TURN+1; q--) {
    for (int p=max(0,q-_maximal_span-1); p<= q-TURN; p++) {
      int type = 0;
      int type2 = 0;

      double temp = 0; bool flag = 0;
      if(p != 0 && q != _seq_length){
	//Beta_stemend
	_Beta_stemend[p][q-p] = q-p >= _maximal_span ? -INF : _Beta_stem[p-1][q-p+2];
	
	//Beta_Multi
	flag = 0;
	if(q-p+1 <= _maximal_span+1){
	  if(_Beta_multi[p-1][q-p+1] != -INF){
	    temp = _Beta_multi[p-1][q-p+1] + MLbase;
	    flag = 1;
	  }
	}
	
	type = BP_pair[_int_sequence[p]][_int_sequence[q+1]];
	int tt = rtype[type];
	if(flag == 1){
	  if(_Beta_stemend[p][q-p] != -INF){
	    temp = logsumexp(temp,_Beta_stemend[p][q-p]+MLclosing+MLintern+ dangle3[tt][_int_sequence[p+1]]+dangle5[tt][_int_sequence[q]]);
	  }
	}else{
	  if(_Beta_stemend[p][q-p] != -INF){
	    temp = _Beta_stemend[p][q-p]+MLclosing+MLintern+dangle3[tt][_int_sequence[p+1]]+dangle5[tt][_int_sequence[q]];
	  }else{
	    temp = -INF;
	  }
	}
	_Beta_multi[p][q-p] = temp;
	
	//Beta_Multi1
	temp = 0; flag = 0;
	for(int k = q+1 ; k<= min(_seq_length,p+_maximal_span);k++){
	  if(_Beta_multibif[p][k-p] != -INF && _Alpha_multi2[q][k-q] != -INF){
	    temp = flag == 0 ? _Beta_multibif[p][k-p]+_Alpha_multi2[q][k-q] : logsumexp(temp,_Beta_multibif[p][k-p]+_Alpha_multi2[q][k-q]) ;
	    flag = 1;
	  }
	}
	_Beta_multi1[p][q-p] = flag == 1 ? temp: -INF;
	
	//Beta_Multi2
	temp = 0; flag = 0;
	if(_Beta_multi1[p][q-p] != -INF){
	  temp = _Beta_multi1[p][q-p];
	  flag = 1;
	}
	if(q-p <= _maximal_span){
	  if(_Beta_multi2[p][q-p+1] != -INF){
	    temp = flag == 1 ? logsumexp(temp,_Beta_multi2[p][q-p+1]+MLbase) : _Beta_multi2[p][q-p+1]+MLbase;
	    flag = 1;
	  }
	}
	
	for(int k = max(0,q-_maximal_span); k < p ;k++){
	  if(_Beta_multibif[k][q-k] != -INF && _Alpha_multi1[k][p-k] != -INF){
	    temp = flag == 0 ? _Beta_multibif[k][q-k]+_Alpha_multi1[k][p-k] : logsumexp(temp,_Beta_multibif[k][q-k]+_Alpha_multi1[k][p-k]);
	    flag = 1;
	  }
	}
	_Beta_multi2[p][q-p] = flag == 0 ? -INF : temp;
	
	//Beta_multibif
	if(_Beta_multi1[p][q-p] != -INF && _Beta_multi[p][q-p] != -INF){
	  _Beta_multibif[p][q-p] = logsumexp(_Beta_multi1[p][q-p],_Beta_multi[p][q-p]);
	}else if(_Beta_multi[p][q-p] == -INF){
	  _Beta_multibif[p][q-p] = _Beta_multi1[p][q-p];
	}else if(_Beta_multi1[p][q-p] == -INF){
	  _Beta_multibif[p][q-p] = _Beta_multi[p][q-p];
	}else{
	  _Beta_multibif[p][q-p] = -INF;
	}
	
      }
      
      //Beta_stem
      type2 = BP_pair[_int_sequence[p+1]][_int_sequence[q]];
      if(type2 != 0){
	temp = _Alpha_outer[p]+_Beta_outer[q]+CalcDangleEnergy(type2,p,q);
	
	type2 = rtype[type2];
	for (int i=max(1,p-MAXLOOP); i<=p; i++){
	  for (int j=q; j<=min(q+ MAXLOOP -p+i,_seq_length-1); j++) {
	    type = BP_pair[_int_sequence[i]][_int_sequence[j+1]];
	    if (type != 0 && !(i == p && j == q)) {
	      if(j-i <= _maximal_span+1 && _Beta_stemend[i][j-i] != -INF){
		temp = logsumexp(temp,_Beta_stemend[i][j-i]+LoopEnergy(type,type2,i,j+1,p+1,q));
	      }
	    }
	  }
	}
	
	if(p != 0 && q != _seq_length){
	  type = BP_pair[_int_sequence[p]][_int_sequence[q+1]];
	  if(type != 0){
	    if(q-p+2 <= _maximal_span+1 && _Beta_stem[p-1][q-p+2] != -INF){
	      temp = logsumexp(temp,_Beta_stem[p-1][q-p+2]+LoopEnergy(type,type2,p,q+1,p+1,q));
	    }
	  }
	}
	_Beta_stem[p][q-p] = temp;
	
	if(_Beta_multi2[p][q-p] != -INF){
	  type2 = rtype[type2];
	  temp = _Beta_multi2[p][q-p] + MLintern + CalcDangleEnergy(type2,p,q);
	  _Beta_stem[p][q-p] = logsumexp(temp,_Beta_stem[p][q-p]);
	}
      }else{
	_Beta_stem[p][q-p] = -INF;
      }
    }
  }
}

double Raccess::logsumexp(double x,double y){
  double temp = x > y ? x + log(exp(y-x) + 1.0) : y + log(exp(x-y) + 1.0) ;
  return(temp);
}

void Raccess::CalcAccessibility(){
  double prob = 0.0;
  float accessibility = 0.0;
  vector<float> accessibility_array; accessibility_array.resize(_seq_length, 0.0);

  vector<double> biloop_probability; biloop_probability.resize(_seq_length, 0.0);
  vector<double> conditional_biloop_probability; conditional_biloop_probability.resize(_seq_length, 0.0);
  vector<double> hairpin_probability; hairpin_probability.resize(_seq_length, 0.0);
  vector<double> conditional_hairpin_probability; conditional_hairpin_probability.resize(_seq_length, 0.0);

  double pf = _Alpha_outer[_seq_length];
  if(pf >= -690 && pf <= 690){
    CalcBulgeAndInternalProbability(biloop_probability, conditional_biloop_probability);
  }else{
    CalcLogSumBulgeAndInternalProbability(biloop_probability, conditional_biloop_probability);
  }
  
  CalcHairpinProbability(hairpin_probability, conditional_hairpin_probability);
  
  ofstream of(_db_name.c_str(), ios::out | ios::binary | ios::app);
  int count = _seq_length - _min_accessible_length + 1;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  for(int i = 1; i <=count;i++){
    prob += CalcExteriorProbability(i,_min_accessible_length);
    prob += hairpin_probability[i-1];
    prob += biloop_probability[i-1];
    prob += CalcMultiProbability(i,_min_accessible_length);
    accessibility = (-log(prob)*kT)/1000;
    //cout << accessibility << endl;
    accessibility_array[i-1] = accessibility;
    of.write(reinterpret_cast<const char*>(&accessibility), sizeof(float));
    prob = 0.0;
  }
 
  count = _seq_length;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  float conditional_accessibility = 0.0;
  for(int i = 0; i < _min_accessible_length; i++){
    of.write(reinterpret_cast<const char*>(&conditional_accessibility), sizeof(float));
  }
  
  for(int i = 1; i+_min_accessible_length <=_seq_length;i++){
    prob += CalcExteriorProbability(i,_min_accessible_length+1);
    prob += conditional_hairpin_probability[i-1];
    prob += conditional_biloop_probability[i-1];
    prob += CalcMultiProbability(i,_min_accessible_length+1);
    conditional_accessibility = (-log(prob)*kT)/1000 - accessibility_array[i-1];
    of.write(reinterpret_cast<const char*>(&conditional_accessibility), sizeof(float));
    prob = 0.0;
  }
  of.close();
}

void Raccess::CalcAccessibility(vector<float> &accessibility, vector<float> &conditional_accessibility){
  double prob = 0.0;
  accessibility.resize(_seq_length, 0.0);
  conditional_accessibility.resize(_seq_length, 0.0);

  vector<double> biloop_probability; biloop_probability.resize(_seq_length, 0.0);
  vector<double> conditional_biloop_probability; conditional_biloop_probability.resize(_seq_length, 0.0);

  vector<double> hairpin_probability; hairpin_probability.resize(_seq_length, 0.0);
  vector<double> conditional_hairpin_probability; conditional_hairpin_probability.resize(_seq_length, 0.0);

  double pf = _Alpha_outer[_seq_length];
  if(pf >= -690 && pf <= 690){
    CalcBulgeAndInternalProbability(biloop_probability, conditional_biloop_probability);
  }else{
    CalcLogSumBulgeAndInternalProbability(biloop_probability, conditional_biloop_probability);
  }
  CalcHairpinProbability(hairpin_probability, conditional_hairpin_probability);

  for(int i = 1; i+_min_accessible_length-1 <=_seq_length;i++){
    prob += CalcExteriorProbability(i,_min_accessible_length);
    prob += hairpin_probability[i-1];
    prob += biloop_probability[i-1];
    prob += CalcMultiProbability(i,_min_accessible_length);
    accessibility[i-1] = (float)(-log(prob)*kT)/1000;
    prob = 0.0;
  }
  
  for(int i = 1; i+_min_accessible_length-1 < _seq_length;i++){
    prob += CalcExteriorProbability(i,_min_accessible_length+1);
    prob += conditional_hairpin_probability[i-1];
    prob += conditional_biloop_probability[i-1];
    prob += CalcMultiProbability(i,_min_accessible_length+1);
    conditional_accessibility[i+_min_accessible_length-1] = (float)(-log(prob)*kT)/1000 - accessibility[i-1];
    prob = 0.0;
  }
}

double Raccess::CalcExteriorProbability(int x, int w){
  double probability = exp(_Alpha_outer[x-1]+_Beta_outer[x+w-1]-_Alpha_outer[_seq_length]);
  return(probability);
}

void Raccess::CalcHairpinProbability(vector<double> &hairpin_probability, vector<double> &conditional_hairpin_probability){
  int w = _min_accessible_length;
  for(int x = 1; x+w-1 <=_seq_length;x++){
    double temp = 0.0;
    double c_temp = 0.0;
    int type = 0;
    bool flag = 0;
    bool c_flag = 0;
    double h_energy = 0.0;
    
    for(int i = max(1,x-_maximal_span);i<x ;i++){
      for(int j = x+w; j<=min(i+_maximal_span,_seq_length);j++){
	type = BP_pair[_int_sequence[i]][_int_sequence[j]];
	if(_Beta_stemend[i][j-i-1] != -INF){
	  h_energy = _Beta_stemend[i][j-i-1] + HairpinEnergy(type, i,j);
	  if(j == x+w){
	    temp = flag == 1 ? logsumexp(temp, h_energy) : h_energy;
	    flag = 1;
	  }else{
	    c_temp = c_flag == 1 ? logsumexp(c_temp, h_energy) :  h_energy;
	    c_flag = 1;
	  }
	}
      }
    }

    if(flag == 1 && c_flag == 1){
      temp = logsumexp(temp, c_temp);
    }
    if(flag == 1){
      hairpin_probability[x-1] = exp(temp-_Alpha_outer[_seq_length]);
    }
    if(c_flag == 1){
      conditional_hairpin_probability[x-1] = exp(c_temp-_Alpha_outer[_seq_length]);
    }
  }
}

double Raccess::CalcMultiProbability(int x, int w){
  double probability = 0.0;
  double temp = 0.0;
  bool flag = 0;
  
  for(int i = x+w-1; i<=min(x+_maximal_span,_seq_length);i++){
    if(_Beta_multi[x-1][i-x+1] != -INF && _Alpha_multi[x+w-1][i-x-w+1] != -INF){
      temp = flag == 0 ? _Beta_multi[x-1][i-x+1] + _Alpha_multi[x+w-1][i-x-w+1] : logsumexp(temp,_Beta_multi[x-1][i-x+1] + _Alpha_multi[x+w-1][i-x-w+1]);
      flag = 1;
    }
  }
  
  for(int i = max(0,x+w-1-_maximal_span); i<x;i++){
    if(_Beta_multi2[i][x+w-1-i] != -INF && _Alpha_multi2[i][x-i-1] != -INF){
      temp = flag == 0 ? _Beta_multi2[i][x+w-1-i] + _Alpha_multi2[i][x-i-1] : logsumexp(temp,_Beta_multi2[i][x+w-1-i] + _Alpha_multi2[i][x-i-1]);
      flag = 1;
    }
  }
  if(flag == 1){ probability = exp(temp-_Alpha_outer[_seq_length]); }
  return(probability);
}

void Raccess::CalcBulgeAndInternalProbability(vector<double> &biloop_probability, vector<double> &conditional_biloop_probability){
  double probability = 0;
  double temp = 0;
  int type = 0;
  int type2 = 0;
  vector<bool> b_flag_array; b_flag_array.resize(_seq_length,0);
  vector<bool> c_flag_array; c_flag_array.resize(_seq_length,0);
  int w = _min_accessible_length;
  
  for(int i = 1; i<_seq_length-TURN-2;i++){
    for(int j = i+TURN+3; j<=min(i+_maximal_span,_seq_length);j++){
      type = BP_pair[_int_sequence[i]][_int_sequence[j]];
      if (type!=0) {
	for (int p =i+1; p <= min(i+MAXLOOP+1,j-TURN-2); p++) {
	  int u1 = p-i-1;
	  for (int q=max(p+TURN+1,j-MAXLOOP+u1-1); q<j; q++) {
	    type2 = BP_pair[_int_sequence[p]][_int_sequence[q]];
	    if (type2 != 0 && !(p == i+1 && q == j-1)) {
	      type2 = rtype[type2];
	      if(_Beta_stemend[i][j-i-1] != -INF && _Alpha_stem[p-1][q-p+1] != -INF){
		temp = exp(_Beta_stemend[i][j-i-1] + LoopEnergy(type, type2,i,j,p,q)+_Alpha_stem[p-1][q-p+1]);
		
		for(int k = i+1; k <= p-w;k++){
		  if(k == p-w){
		    biloop_probability[k-1] += temp;		   
		  }else{
		    conditional_biloop_probability[k-1] += temp;		  
		  }
		}
		
		for(int k = q+1; k <= j-w;k++){
		  if(k == j-w){
		    biloop_probability[k-1] += temp;
		  }else{
		    conditional_biloop_probability[k-1] += temp;
		  }
		}
	      } 
	    }
	  }
	}
      }
    }
  }
  
  for(int i=0;i<_seq_length;i++){
    if(biloop_probability[i] != 0){
      biloop_probability[i] = log(biloop_probability[i] + conditional_biloop_probability[i]);
      biloop_probability[i] = exp(biloop_probability[i]-_Alpha_outer[_seq_length]);
    }
    if(conditional_biloop_probability[i] != 0){
      conditional_biloop_probability[i] = log(conditional_biloop_probability[i]);
      conditional_biloop_probability[i] = exp(conditional_biloop_probability[i]-_Alpha_outer[_seq_length]);
    }
  }
}

void Raccess::CalcLogSumBulgeAndInternalProbability(vector<double> &biloop_probability, vector<double> &conditional_biloop_probability){
  double probability = 0;
  double temp = 0;
  int type = 0;
  int type2 = 0;
  vector<bool> b_flag_array; b_flag_array.resize(_seq_length,0);
  vector<bool> c_flag_array; c_flag_array.resize(_seq_length,0);
  int w = _min_accessible_length;
  
  for(int i = 1; i<_seq_length-TURN-2;i++){
    for(int j = i+TURN+3; j<=min(i+_maximal_span,_seq_length);j++){
      type = BP_pair[_int_sequence[i]][_int_sequence[j]];
      if (type!=0) {
	for (int p =i+1; p <= min(i+MAXLOOP+1,j-TURN-2); p++) {
	  int u1 = p-i-1;
	  for (int q=max(p+TURN+1,j-MAXLOOP+u1-1); q<j; q++) {
	    type2 = BP_pair[_int_sequence[p]][_int_sequence[q]];
	    if (type2 != 0 && !(p == i+1 && q == j-1)) {
	      type2 = rtype[type2];
	      if(_Beta_stemend[i][j-i-1] != -INF && _Alpha_stem[p-1][q-p+1] != -INF){
		temp = _Beta_stemend[i][j-i-1] + LoopEnergy(type, type2,i,j,p,q)+_Alpha_stem[p-1][q-p+1];
		
		for(int k = i+1; k <= p-w;k++){
		  if(k == p-w){
		    biloop_probability[k-1] = b_flag_array[k-1] == 1 ? logsumexp(biloop_probability[k-1], temp) : temp;
		    b_flag_array[k-1] = 1;		   
		  }else{
		    conditional_biloop_probability[k-1] = c_flag_array[k-1] == 1 ? logsumexp(conditional_biloop_probability[k-1], temp) : temp;
		    c_flag_array[k-1] = 1;		  
		  }
		}
		
		for(int k = q+1; k <= j-w;k++){
		  if(k == j-w){
		    biloop_probability[k-1] = b_flag_array[k-1] == 1 ? logsumexp(biloop_probability[k-1], temp) : temp;
		    b_flag_array[k-1] = 1;
		  }else{
		    conditional_biloop_probability[k-1] = c_flag_array[k-1] == 1 ? logsumexp(conditional_biloop_probability[k-1], temp) : temp;
		    c_flag_array[k-1] = 1;
		  }
		}
	      } 
	    }
	  }
	}
      }
    }
  }
  
  for(int i=0;i<_seq_length;i++){
    if(b_flag_array[i]==1 && c_flag_array[i]==1){
      biloop_probability[i] = logsumexp(biloop_probability[i], conditional_biloop_probability[i]);
    }
    if(b_flag_array[i]==1){
      biloop_probability[i] = exp(biloop_probability[i]-_Alpha_outer[_seq_length]);
    }
    if(c_flag_array[i]==1){
      conditional_biloop_probability[i] = exp(conditional_biloop_probability[i]-_Alpha_outer[_seq_length]);
    }
  }
}

double Raccess::LoopEnergy(int type, int type2,int i,int j,int p,int q){
  double z=0;
  int u1 = p-i-1;
  int u2 = j-q-1;
  
  if ((u1==0) && (u2==0)){
    z = stack[type][type2];
  }else{
    if ((u1==0)||(u2==0)) {
      int u;
      u = u1 == 0 ? u2 : u1;
      z = u <=30 ? bulge[u] : bulge[30] - lxc37*log( u/30.)*10./kT;
      
      if (u == 1){
	z += stack[type][type2];
      }else {
	if (type>2){ z += TermAU;}
	if (type2>2){ z += TermAU;}
      }
    }else{     
      if (u1+u2==2) {
	z = int11[type][type2][_int_sequence[i+1]][_int_sequence[j-1]];
      }else if ((u1==1) && (u2==2)){
	z = int21[type][type2][_int_sequence[i+1]][_int_sequence[q+1]][_int_sequence[j-1]];
      }else if ((u1==2) && (u2==1)){
	z = int21[type2][type][_int_sequence[q+1]][_int_sequence[i+1]][_int_sequence[p-1]];
      }else if ((u1==2) && (u2==2)){
	z = int22[type][type2][_int_sequence[i+1]][_int_sequence[p-1]][_int_sequence[q+1]][_int_sequence[j-1]];
      }else{
	z = internal[u1+u2]+mismatchI[type][_int_sequence[i+1]][_int_sequence[j-1]]+mismatchI[type2][_int_sequence[q+1]][_int_sequence[p-1]];
	z += ninio[abs(u1-u2)];
      }
    }
  }
  return z;
}

double Raccess::HairpinEnergy(int type, int i, int j) {
  int d = j-i-1;
  double q = 0;

  q = d <= 30 ? hairpin[d] : hairpin[30] - lxc37*log( d/30.) *10./kT;  
  if(d!= 3){
    q += mismatchH[type][_int_sequence[i+1]][_int_sequence[j-1]];
  }else{
    if(type > 2){q += TermAU;}
  }
  return q;
}

void Raccess::Clear(){
  for(int i = 0; i <= _seq_length;i++){
    _Alpha_stem[i].clear();
    _Alpha_stemend[i].clear();
    _Alpha_multi[i].clear();
    _Alpha_multibif[i].clear();
    _Alpha_multi1[i].clear();
    _Alpha_multi2[i].clear();
    _Beta_stem[i].clear();
    _Beta_stemend[i].clear();
    _Beta_multi[i].clear();
    _Beta_multibif[i].clear();
    _Beta_multi1[i].clear();
    _Beta_multi2[i].clear();
  }

  _int_sequence.clear();
  _seq_length = 0;
  _Alpha_outer.clear();
  _Alpha_stem.clear();
  _Alpha_stemend.clear();
  _Alpha_multi.clear();
  _Alpha_multibif.clear();
  _Alpha_multi1.clear();
  _Alpha_multi2.clear();
  
  _Beta_outer.clear();
  _Beta_stem.clear();
  _Beta_stemend.clear();
  _Beta_multi.clear();
  _Beta_multibif.clear();
  _Beta_multi1.clear();
  _Beta_multi2.clear();
}
