// Q20 exact canonical markers, deduplicated per fragment and evidence group.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
int base(char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1;}}
int main(int argc,char**argv){
 if(argc!=4){cerr<<"marker_groups.tsv counts.bin groups.tsv; interleaved FASTQ on stdin\n";return 2;}
 ifstream file(argv[1]);if(!file)return 2;uint64_t key;unsigned mask;vector<unsigned> masks;unordered_map<uint64_t,size_t> idx;
 while(file>>key>>mask){idx.emplace(key,idx.size());masks.push_back(mask);}
 vector<uint32_t> counts(idx.size(),0);uint64_t groups[7]={},aonly=0,bonly=0,otheronly=0,abambig=0,lsambig=0,fragments=0;unordered_set<size_t> hit;string prev,name,s,plus,q;
 auto flush=[&](){if(prev.empty())return;unsigned bits=0;for(auto i:hit){counts[i]++;bits|=masks[i];}for(int j=0;j<7;j++)if(bits&(1<<j))groups[j]++;if((bits&67)==1)aonly++;if((bits&67)==2)bonly++;if((bits&67)==64)otheronly++;if((bits&3)==3)abambig++;if((bits&16)&&(bits&12))lsambig++;fragments++;hit.clear();};
 while(getline(cin,name)){
  if(!getline(cin,s)||!getline(cin,plus)||!getline(cin,q)||s.size()!=q.size()||name.empty()||name[0]!='@')return 3;
  name=name.substr(1,name.find_first_of(" \t")-1);if(name.size()>2&&name[name.size()-2]=='/'&&(name.back()=='1'||name.back()=='2'))name.resize(name.size()-2);
  if(name!=prev){flush();prev=name;}uint64_t f=0,r=0;int valid=0;
  for(size_t j=0;j<s.size();j++){int b=base(s[j]);if(b<0||q[j]-33<20){f=r=0;valid=0;continue;}f=((f<<2)|b)&((1ULL<<62)-1);r=(r>>2)|(uint64_t(3-b)<<60);if(++valid>=31){auto it=idx.find(min(f,r));if(it!=idx.end())hit.insert(it->second);}}
 }flush();ofstream out(argv[2],ios::binary);for(auto c:counts)for(int j=0;j<4;j++)out.put(char((c>>(8*j))&255));if(!out)return 4;
 ofstream g(argv[3]);g<<"metric\tfragments\n";string labels[]={"A_any","B_any","LONG_LEFT","LONG_RIGHT","SHORT","BOUNDARY","NONCANONICAL"};for(int i=0;i<7;i++)g<<labels[i]<<'\t'<<groups[i]<<'\n';g<<"A_only\t"<<aonly<<"\nB_only\t"<<bonly<<"\nOTHER_only\t"<<otheronly<<"\nAB_ambiguous\t"<<abambig<<"\nLS_ambiguous\t"<<lsambig<<"\nall_fragments\t"<<fragments<<'\n';return g?0:4;
}
