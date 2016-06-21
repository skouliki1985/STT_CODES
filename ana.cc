
//#include <TApplication.h>
#include "myevent.h"
#include <iostream>

using namespace std;

int main(){
	myevent t;
	cout << " *******SOS******* Start of processing the events!!!!!! " << endl;
	t.Loop();
	t.Correction();
	
}
