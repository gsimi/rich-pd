#!/bin/bash
cp PDconfig.c pmtdaq.c Makefile.pmtdaq ~/CAEN/Wavedump/wavedump-3.5.3/src/
cp pmtdaq.h ~/CAEN/Wavedump/wavedump-3.5.3/inc/
cd ~/CAEN/Wavedump/wavedump-3.5.3/src/
make -f Makefile.pmtdaq pmtdaq 
cp pmtdaq /home/lhcb/rich-pd/pmt-test-facility/daq
cd -
