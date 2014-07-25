#include "dt5533.hh"
#include "CAENHVWrapper.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

DT5533::DT5533(unsigned short port) : m_nchans(4) {
  m_errors[0] = "No error";
  m_errors[1] = "Operating system error";
  m_errors[2] = "Writing error";
  m_errors[3] = "Reading error";
  m_errors[4] = "Time out error";
  m_errors[5] = "Command Front End application is down";
  m_errors[6] = "Communication with system not yet connected by a Login command";
  m_errors[7] = "Execute Command not yet implemented";
  m_errors[8] = "Get Property not yet implemented";
  m_errors[9] = "Set Property not yet implemented";
  m_errors[10] = "Communication with RS232 not yet implemented";
  m_errors[11] = "User memory not sufficient";
  m_errors[12] = "Value out of range";
  m_errors[13] = "Property not yet implemented";
  m_errors[14] = "Property not found";
  m_errors[15] = "Command not found";
  m_errors[16] = "Not a Property";
  m_errors[17] = "Not a reading Property";
  m_errors[18] = "Not a writing Property";
  m_errors[19] = "Not a Command";
  m_errors[20] = "configuration change";
  m_errors[21] = "Parameter’s Property not found";
  m_errors[22] = "Parameter not found";
  m_errors[23] = "No data present";
  m_errors[24] = "Device already open";
  m_errors[25] = "Too Many devices opened";
  m_errors[26] = "Function Parameter not valid";
  m_errors[27] = "Function not available for the connected device";
  m_errors[28] = "SOCKET ERROR";
  m_errors[29] = "COMMUNICATION ERROR";
  m_errors[30] = "NOT YET IMPLEMENTED";
  m_errors[31] = "OTHER error";

  m_status[0] = "on";
  m_status[1] = "ramping up";
  m_status[2] = "ramping down";
  m_status[3] = "overcurrent";
  m_status[4] = "overvoltage";
  m_status[5] = "undervoltage";
  m_status[6] = "external trip";
  m_status[7] = "max V";
  m_status[8] = "external disable";
  m_status[9] = "internal trip";
  m_status[10] = "calibration error";
  m_status[11] = "unplugged";

  for (unsigned short i=0; i<m_nchans; ++i) {
    m_channels[i]=i;
  }

  std::ostringstream arg;
  arg << port << "_0_00000000";
//   std::cout << arg.str() << "\tport" << std::endl;
  CAENHVRESULT ret = CAENHV_InitSystem(DT55XX, LINKTYPE_USB, (void*)(arg.str().c_str()), "", "", &m_handle);
  if (ret != 0) {
    std::cout << std::hex <<ret << std::dec << "\terror!!!" << std::endl;
    if (ret<32) {
      std::cout << "DT5533 error initializing class: " << m_errors[ret] << endl;
    }
  }
}

DT5533::~DT5533() {
  CAENHVRESULT ret = CAENHV_DeinitSystem(m_handle);
  if (ret != 0) {
    std::cout << "DT5533 error destroing class: " << m_errors[ret] << endl;
  }
}

void DT5533::PrintChannels() const {
  char  (*listNameCh)[MAX_CH_NAME];
  listNameCh = (char(*)[MAX_CH_NAME]) malloc(m_nchans*MAX_CH_NAME);

  CAENHVRESULT ret = CAENHV_GetChName(m_handle, 0, m_nchans, m_channels, listNameCh);
  if( ret != CAENHV_OK ) {
    free(listNameCh);
    printf("CAENHV_GetChName: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  } else {
    printf("CHANNEL NAME\n");
    for(unsigned short n = 0; n < m_nchans; n++ )
      printf("Channel n. %d: %s\n", m_channels[n], listNameCh[n]);
  }
  free(listNameCh);
}

void DT5533::PrintChInfo(unsigned short ch) const {
  int nparms=0;
  char  *parmlist;
//   parmlist = (char(**)) malloc(100*1024);
  CAENHVRESULT ret = CAENHV_GetChParamInfo(m_handle, 0, ch, &parmlist, &nparms);
  if( ret != CAENHV_OK ) {
    printf("CAENHV_GetChParamInfo: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  } else {
//     printf("CHANNEL NAME\n");
//     for(unsigned short n = 0; n < m_nchans; n++ )
//       printf("Channel n. %d: %s\n", m_channels[n], listNameCh[n]);
//     cout << nparms << endl;
    char* p=parmlist;
    for (unsigned int i=0; i<nparms; p+=strlen(p)+1) {
      if( *p != '\0' ) {
        cout << p << " " << endl;//strlen(p) << " " << strlen(parmlist) << endl;
//         cout << "\t----- " << i << endl;
        ++i;
        int iretval;
        float fretval;
        ret = CAENHV_GetChParamProp(m_handle, 0, ch, p, "Minval", (void*)&fretval);
        cout << "\t" << "Minval  " << (float)fretval << endl;
        ret = CAENHV_GetChParamProp(m_handle, 0, ch, p, "Maxval", (void*)&fretval);
        cout << "\t" << "Maxval  " << (float)fretval << endl;
        ret = CAENHV_GetChParamProp(m_handle, 0, ch, p, "Unit", (void*)&iretval);
        cout << "\t" << "Unit  " << (unsigned short)iretval << endl;
        ret = CAENHV_GetChParamProp(m_handle, 0, ch, p, "Exp", (void*)&iretval);
        cout << "\t" << "Exp  " << (short)iretval << endl;
        ret = CAENHV_GetChParamProp(m_handle, 0, ch, p, "Decimal", (void*)&iretval);
        cout << "\t" << "Decimal  " << (unsigned short)iretval << endl;
      }

    }
    free(parmlist);
  }
}

void DT5533::PrintStatus(std::ostream& out, unsigned short ch) const {
  out << *this << endl;

}

void DT5533::TurnOn(unsigned short ch) {
  int state=1;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "Pw", 1, &ch, (void*)&state);
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

void DT5533::TurnOff(unsigned short ch) {
  int state=0;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "Pw", 1, &ch, (void*)&state);
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

void DT5533::SetV(unsigned short ch, float v) {
  unsigned short channel[1];
  channel[0]=ch;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "VSet", 1, channel, (void*)&v);
//   cout << ch << " " << v << endl;
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

void DT5533::SetI(unsigned short ch, float i) {
  unsigned short channel[1];
  channel[0]=ch;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "ISet", 1, &ch, (void*)&i);
//   cout << ch << " " << i << endl;
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

void DT5533::SetRampUp(unsigned short ch, float r) {
  unsigned short channel[1];
  channel[0]=ch;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "RampUp", 1, &ch, (void*)&r);
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

void DT5533::SetRampDown(unsigned short ch, float r) {
  unsigned short channel[1];
  channel[0]=ch;
  CAENHVRESULT ret = CAENHV_SetChParam(m_handle, 0, "RampDwn", 1, &ch, (void*)&r);
  if( ret != CAENHV_OK ) {
    printf("CAENHV_SetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
  }
}

float DT5533::GetFParm(unsigned short ch, const char* parm) const {
  float f;
  CAENHVRESULT ret = CAENHV_GetChParam(m_handle, 0, parm, 1, &ch, (void*)&f);

  if( ret != CAENHV_OK ) {
    printf("CAENHV_GetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
    return -1;
  }
  return f;
}

int DT5533::GetIParm(unsigned short ch, const char* parm) const {
  int i;
  CAENHVRESULT ret = CAENHV_GetChParam(m_handle, 0, parm, 1, &ch, (void*)&i);

  if( ret != CAENHV_OK ) {
    printf("CAENHV_GetChParam: %s (num. 0x%08x)\n\n", CAENHV_GetError(m_handle), ret);
    std::cout << "DT5533 error: " << m_errors[ret] << endl;
    return -1;
  }
  return i;
}

int DT5533::GetStatus(unsigned short ch, bool verbose) const {
  int status = GetIParm(ch, "Status");
  if (status!=0 && status !=1) {
    int s=status;
    int m=1;
//     cout << std::hex << status << std::dec << endl;
    for (unsigned int i=0; i<12; ++i) {
//       cout << i << " " << (s>>i) << " - " << ((s>>i)&m) << " + " << std::hex << status << std::dec << endl;
      int k =(s>>i)&m;
      if (verbose && k!=0) cout << m_status[i] << endl;
    }

  }
  return status;
}


void DT5533::Status(unsigned short ch) const {
  float is = GetFParm(ch, "ISet");
  float vs = GetFParm(ch, "VSet");
  float im = GetFParm(ch, "IMon");
  float vm = GetFParm(ch, "VMon");
  float ru = GetFParm(ch, "RampUp");
  float rd = GetFParm(ch, "RampDwn");
  int status = GetIParm(ch, "Status");
  int power = GetIParm(ch, "Pw");

  cout << "Channel: " << ch << " " << " VSet " << vs  << " ISet " << is << " VMon " << vm
       << " IMon " << im << " ramp UP " << ru << " ramp down " << rd << " Status " << std::hex
      << status << std::dec << " power " << power << endl;

}


std::ostream & operator <<(std::ostream& output, const DT5533 & dt) {

  const char onoff[4][20] = { "\033[22;32mOFF\033[01;37m", "\033[22;31mON\033[01;37m", "Up", "Down" };
  for (unsigned short ch=0; ch<4; ++ch) {
    int status = dt.GetStatus(ch);
    float vs = dt.GetFParm(ch,"VSet");
    float v = dt.GetFParm(ch,"VMon");
    float c = dt.GetFParm(ch,"IMon");
    std::string st=onoff[(status&1)];
    if ((status&6) != 0) {
      if ((status&2) == 0) {
        st=onoff[3];
      } else if ((status&4) == 0) {
        st=onoff[2];
      }
    }
    std::string s("Channel #");
    std::ostringstream oss(s, std::ostringstream::ate);
    oss << ch  << " ---   Vset: " << std::setw(8) << std::fixed << std::setprecision(1) << vs;
    oss << " Vmon: " << std::setw(8) << std::fixed << std::setprecision(1) << v;
    oss << " Imon: " << std::setw(8) << std::fixed << std::setprecision(2) << c;
    oss << "   " << std::setfill(' ') << std::fixed << std::setw(15) << st;
    oss << " status: " << std::setw(15) << std::hex << status;
    output << oss.str() << "\n";
  }
//   output << "\n\n\n1 " << (status&1) << " 2 " << (status&2) << " 4 " << (status&4) << " 6 " << (status&6) << " == " << std::hex << status << endl;
  return output;
}
