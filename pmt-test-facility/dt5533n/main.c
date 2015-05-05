#include "dt5533.hh"
#include "CAENHVWrapper.h"
#include "CAENComm.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(){

  unsigned short n=0;
  DT5533 *hvps = new DT5533(n);
  hvps->PrintChannels();
  hvps->Status(0);
  hvps->Status(1);
  hvps->Status(2);
  hvps->Status(3);
  /* hvps->TurnOn(0); */
  /* hvps->Status(0); */
 /* hvps->TurnOff(0);		 */
  /* hvps->Status(0); */
  delete hvps; 
      /*
  unsigned short port=0;
  std::ostringstream arg;
  arg << port << "_0_00000000";
  int m_handle;
  CAENHVRESULT ret = CAENHV_InitSystem(DT55XX, LINKTYPE_USB, (void*)(arg.str().c_str()), "", "", &m_handle);
      */  

}
