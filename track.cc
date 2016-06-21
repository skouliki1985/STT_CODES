#include "basictracking.h"
#include <iostream>
#include <TApplication.h>
using namespace std;


//int main(int argc, char **argv){
int main(){
//	TApplication theApp("App",&argc,argv);
	basictracking b;
	cout << "START!" << endl;
	b.track();
  //  theApp.Run(kTRUE);
   // theApp.Terminate();

}
