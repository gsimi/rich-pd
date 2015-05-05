/*
#include <iostream>

using namespace std;

int main()
{
    cout << "Hello world!" << endl;
    return 0;
}*/

#include "Serial.h"


int main(int argc, char* argv[]){
	CSerial serial;
	if (serial.Open(2, 9600))
		AfxMessageBox("Port opened successfully");
	else
		AfxMessageBox("Failed to open port!");
}
/*
int main(int argc, char* argv[]){
  CSerial serial;
  if (serial.Open(2, 9600)){
    AfxMessageBox("Port opened successfully");
    
    static char* szMessage[] = "This is test data";
    int nBytesSent = serial.SendData(szMessage, strlen(szMessage));
    ASSERT(nBytesSent == strlen(szMessage));
    
    AfxMessageBox("Failed to open port!");
    //  Reading Data
    //    CSerial serial;
    
    char* lpBuffer = new char[500];
    int nBytesRead = serial.ReadData(lpBuffer, 500);
    delete []lpBuffer;
  }
  else
    AfxMessageBox("Failed to open port!");
  
}
*/