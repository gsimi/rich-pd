#ifndef DT5533_hh
#define DT5533_hh

#include <string>
#include <iostream>

class DT5533 {
private:
  int m_handle;
  unsigned short m_nchans;
  unsigned short m_channels[32];
  std::string m_errors[32];
  std::string m_status[12];


public:
  DT5533(unsigned short port);
  ~DT5533();

  void PrintChannels(void) const;
  void PrintChInfo(unsigned short ch) const;
  void PrintStatus(std::ostream& os, unsigned short ch) const;
  void Status(unsigned short ch) const;
  int GetStatus(unsigned short ch, bool verbose=false) const;
  float GetFParm(unsigned short ch, const char* parm) const;
  int GetIParm(unsigned short ch, const char* parm) const;
  void TurnOn(unsigned short ch);
  void TurnOff(unsigned short ch);
  void SetV(unsigned short ch, float v);
  void SetI(unsigned short ch, float i);
  void SetRampUp(unsigned short ch, float r);
  void SetRampDown(unsigned short ch, float r);
  void PrintInfo(const DT5533& hv, int ch);

  friend std::ostream & operator << (std::ostream& output, const DT5533& dt);
};

#endif