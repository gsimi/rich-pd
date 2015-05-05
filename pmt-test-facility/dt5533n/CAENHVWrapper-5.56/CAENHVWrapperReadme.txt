/*****************************************************************************/
/*                                                                           */
/*                  --- CAEN SpA - Computing Division ---                    */
/*                                                                           */
/*   CAEN HV Wrapper Library Rel. 5.56 Installation and Use Instructions     */
/*		                                                             */
/*   January  2014		                                             */
/*                                                                           */
/*****************************************************************************/
  
 This archive contains the last release of the CAEN HV Wrapper Library and the
 corresponding Release Notes.

 The complete documentation can be found in the CAEN HV Wrapper Library Manual
 available once installed the Library or on CAEN's Web Site 
 at www.caen.it.


 Content of the archive
 ----------------------

 CAENUVWrapperReadme.txt        :  This file
 CAENHVWrapperReleaseNotes.txt  :  Release Notes of the last software release
 Makefile 	         		        :  makefile for library installation and demo 
                                   program
 
 Lib/
  libcaenhvwrapper.so.x.xx    	:  executable of the library (dynamic)
  hscaenetlib.so.x.x          	:  executable of the HSCAENETLib library
  libcaenhvwrapper.x.xx.a	:  executable of the library (static)
  
 Doc/
  CAENHVWrapper.pdf 	        :  user's manual of the library
 
 include/
  CAENHVWrapper.h             :  include file for use of the library
  caenhvoslib.h		    	      :  accessory include file
 
 HVWrapperDemo/               :  directory with sources/executable of the demo 
                                 program 
 


 System Requirements
 -------------------
 
 - Network Interface Card + TCP/IP protocol (to control SY 1527/ SY 2527 / SY 4527 / SY 5527)
 - USB port (to control V65xx Boards, via V1718 VME Bridge)
 - Optical Link (to control V65xx Boards, via V2718 - VME-PCI Optical Link Bridge)
 - A303A/A1303 H.S. CAENET Controller Card (to control SY 527, SY 127, SY 403,
   N 470, N 570 and N568)
 - A128HS SY127 H.S. CAENET Controller installed on SY 127
 - SY 1527/ SY 2527 firmware version 1.10.0 or later (recommended 1.14.03)
 - SY403 firmware version 1.45 or later (recommended 1.46)
 - SY527 firmware version 4.03 or later
 - Linux gcc 2.9 or greater with gnu C/C++ compiler


 Installation notes
 ------------------

 1. It's necessary to login as root
 
 2. execute: make install
 
 The installation copies and installs the library in /usr/lib, compiles the demo 
 program and installs it in the work directory.
 To use the demo program, change to the demo program directory and launch the
 application typing ./HVWrappdemo.
   

 Note:
 -----
 Control of CAEN Power Supplies via CaeNet link requires the correct
 installation of the A1303 device driver. 
 Control of CAEN Power Supplies via USB/OpycalLink link requires the correct
 installation of the USB/A2818 device driver and CAENComm library. 
