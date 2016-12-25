#include <Rcpp.h>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace std;

std::string ssC(int I1,int I2,std::string line,std::string sep){
std::string temp="";
if (I1!=I2){
temp="";
if (line!=""){
if ((I1==0) & (I2==1) & (line[0]==sep[0])){temp="";}
if ((I1!=0) | (I2!=1) | (line[0]!=sep[0])){
int num=0;
bool ok=false;
if (I1==0){ok=true;}
int i;
for(i = 0; i < line.length(); i++)
{
if (line[i]==sep[0]){
num=num+1;
if (I2==num){break;}
if (ok==true){temp=temp+line[i];}
if (I1==num){ok=true;}
}
if (line[i]!=sep[0]){
if (ok==true){temp=temp+line[i];}
}
}
}
}
}
return temp;
}

//' Return subString for each item in list
//'
//' @param I1 integer
//' @param I2 integer
//' @param line string
//' @param sep string
//' @export
// [[Rcpp::export]]
std::vector<std::string> ssLC(int I1,int I2, std::vector<std::string> line,std::string sep) {
int num_strings=line.size();
for( int i=0; i < num_strings; i++ ) {
line[i]=ssC(I1,I2,line[i],sep);
}
return line;
}

