/******************************************************************************
* 
* CAEN SpA - Front End Division
* Via Vetraia, 11 - 55049 - Viareggio ITALY
* +390594388398 - www.caen.it
*
***************************************************************************//**
* \note TERMS OF USE:
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation. This program is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The user relies on the 
* software, documentation and results solely at his own risk.
* -----------------------------------------------------------------------------
* WDconfig contains the functions for reading the configuration file and 
* setting the parameters in the WDcfg structure
******************************************************************************/


#include <CAENDigitizer.h>
#include "pmtdaq.h"


/*! \fn      int ParseConfigFile(FILE *f_ini, WaveDumpConfig_t *WDcfg) 
*   \brief   Read the configuration file and set the WaveDump paremeters
*            
*   \param   f_ini        Pointer to the config file
*   \param   WDcfg:   Pointer to the WaveDumpConfig data structure
*   \return  0 = Success; negative numbers are error codes
*/
int ParseConfigFile(FILE *f_ini, WaveDumpConfig_t *WDcfg) 
{
	char str[1000], str1[1000];
	int i,j, ch=-1, val, Off=0, tr = -1;

	/* Default settings */
	WDcfg->RecordLength = (1024*16);
	WDcfg->PostTrigger = 80;
	WDcfg->NumEvents = 1023;
	WDcfg->EnableMask = 0xFF;
	WDcfg->GWn = 0;
    WDcfg->ExtTriggerMode = CAEN_DGTZ_TRGMODE_ACQ_ONLY;
    WDcfg->InterruptNumEvents = 0;
    WDcfg->TestPattern = 0;
    WDcfg->TriggerEdge = 0;
    WDcfg->DesMode = 0;
	WDcfg->FastTriggerMode = 0; 
    WDcfg->FastTriggerEnabled = 0; 
	WDcfg->FPIOtype = 0;
	strcpy(WDcfg->GnuPlotPath, GNUPLOT_DEFAULT_PATH);
	for(i=0; i<MAX_SET; i++) {
		WDcfg->DCoffset[i] = 0;
		WDcfg->Threshold[i] = 0;
        WDcfg->ChannelTriggerMode[i] = CAEN_DGTZ_TRGMODE_DISABLED;
		WDcfg->GroupTrgEnableMask[i] = 0;
		for(j=0; j<MAX_SET; j++) WDcfg->DCoffsetGrpCh[i][j] = -1;
		WDcfg->FTThreshold[i] = 0;
		WDcfg->FTDCoffset[i] =0;
		WDcfg->GroupWriteMask[i]=0;
    }

    WDcfg->useCorrections = -1;
    WDcfg->EventsToWrite=0;
    WDcfg->PlotMask=0;


	/* read config file and assign parameters */
	while(!feof(f_ini)) {
		int read;
        char *res;
        // read a word from the file
        read = fscanf(f_ini, "%s", str);
        if( !read || (read == EOF) || !strlen(str))
			continue;
        // skip comments
        if(str[0] == '#') {
            res = fgets(str, 1000, f_ini);
			continue;
        }

        if (strcmp(str, "@ON")==0) {
            Off = 0;
            continue;
        }
		if (strcmp(str, "@OFF")==0)
            Off = 1;
        if (Off)
            continue;


        // Section (COMMON or individual channel)
		if (str[0] == '[') {
            if (strstr(str, "COMMON")) {
                ch = -1;
               continue; 
            }
            if (strstr(str, "TR")) {
				sscanf(str+1, "TR%d", &val);
				 if (val < 0 || val >= MAX_SET) {
                    printf("%s: Invalid channel number\n", str);
                } else {
                    tr = val;
                }
            } else {
                sscanf(str+1, "%d", &val);
                if (val < 0 || val >= MAX_SET) {
                    printf("%s: Invalid channel number\n", str);
                } else {
                    ch = val;
                }
            }
            continue;
		}
 
        // OPEN: read the details of physical path to the digitizer
		if (strstr(str, "OPEN")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "USB")==0)
				WDcfg->LinkType = CAEN_DGTZ_USB;
			else if (strcmp(str1, "PCI")==0)
				WDcfg->LinkType = CAEN_DGTZ_PCI_OpticalLink;
            else {
                printf("%s %s: Invalid connection type\n", str, str1);
				return -1; 
            }
			read = fscanf(f_ini, "%d", &WDcfg->LinkNum);
            if (WDcfg->LinkType == CAEN_DGTZ_USB)
                WDcfg->ConetNode = 0;
            else
			    read = fscanf(f_ini, "%d", &WDcfg->ConetNode);
			read = fscanf(f_ini, "%x", &WDcfg->BaseAddress);
			continue;
		}

		// Generic VME Write (address offset + data, both exadecimal)
		if ((strstr(str, "WRITE_REGISTER")!=NULL) && (WDcfg->GWn < MAX_GW)) {
			read = fscanf(f_ini, "%x", (int *)&WDcfg->GWaddr[WDcfg->GWn]);
			read = fscanf(f_ini, "%x", (int *)&WDcfg->GWdata[WDcfg->GWn]);
			WDcfg->GWn++;
			continue;
		}

        // Acquisition Record Length (number of samples)
		if (strstr(str, "RECORD_LENGTH")!=NULL) {
			read = fscanf(f_ini, "%d", &WDcfg->RecordLength);
			continue;
		}

        // Correction Level (mask)
		if (strstr(str, "CORRECTION_LEVEL")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
            if( strcmp(str1, "AUTO") == 0 )
                WDcfg->useCorrections = -1;
            else
                WDcfg->useCorrections = atoi(str1);
			continue;
		}

        // Test Pattern
		if (strstr(str, "TEST_PATTERN")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "YES")==0)
				WDcfg->TestPattern = 1;
			else if (strcmp(str1, "NO")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}

        // Trigger Edge
		if (strstr(str, "TRIGGER_EDGE")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "FALLING")==0)
				WDcfg->TriggerEdge = 1;
			else if (strcmp(str1, "RISING")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}

        // External Trigger (DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT)
		if (strstr(str, "EXTERNAL_TRIGGER")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "DISABLED")==0)
                WDcfg->ExtTriggerMode = CAEN_DGTZ_TRGMODE_DISABLED;
			else if (strcmp(str1, "ACQUISITION_ONLY")==0)
                WDcfg->ExtTriggerMode = CAEN_DGTZ_TRGMODE_ACQ_ONLY;
			else if (strcmp(str1, "ACQUISITION_AND_TRGOUT")==0)
                WDcfg->ExtTriggerMode = CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT;
            else
                printf("%s: Invalid Parameter\n", str);
            continue;
		}

        // Max. number of events for a block transfer (0 to 1023)
		if (strstr(str, "MAX_NUM_EVENTS_BLT")!=NULL) {
			read = fscanf(f_ini, "%d", &WDcfg->NumEvents);
			continue;
		}

		// GNUplot path
		if (strstr(str, "GNUPLOT_PATH")!=NULL) {
			read = fscanf(f_ini, "%s", WDcfg->GnuPlotPath);
			continue;
		}

		// Post Trigger (percent of the acquisition window)
		if (strstr(str, "POST_TRIGGER")!=NULL) {
			read = fscanf(f_ini, "%d", &WDcfg->PostTrigger);
			continue;
		}

        // DesMode (Double sampling frequency for the Mod 731 and 751)
		if (strstr(str, "ENABLE_DES_MODE")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "YES")==0)
				WDcfg->DesMode = 1;
			else if (strcmp(str1, "NO")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}

		// Output file format (BINARY or ASCII)
		if (strstr(str, "OUTPUT_FILE_FORMAT")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "BINARY")==0)
				WDcfg->OutFileFlags|= OFF_BINARY;
			else if (strcmp(str1, "ASCII")!=0)
				printf("%s: invalid output file format\n", str1);
			continue;
		}

		// Header into output file (YES or NO)
		if (strstr(str, "OUTPUT_FILE_HEADER")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "YES")==0)
				WDcfg->OutFileFlags|= OFF_HEADER;
			else if (strcmp(str1, "NO")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}

        // Interrupt settings (request interrupt when there are at least N events to read; 0=disable interrupts (polling mode))
		if (strstr(str, "USE_INTERRUPT")!=NULL) {
			read = fscanf(f_ini, "%d", &WDcfg->InterruptNumEvents);
			continue;
		}
		
		if (!strcmp(str, "FAST_TRIGGER")) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "DISABLED")==0)
                WDcfg->FastTriggerMode = CAEN_DGTZ_TRGMODE_DISABLED;
			else if (strcmp(str1, "ACQUISITION_ONLY")==0)
                WDcfg->FastTriggerMode = CAEN_DGTZ_TRGMODE_ACQ_ONLY;
            else
                printf("%s: Invalid Parameter\n", str);
            continue;
		}
		
		if (strstr(str, "ENABLED_FAST_TRIGGER_DIGITIZING")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "YES")==0)
				WDcfg->FastTriggerEnabled= 1;
			else if (strcmp(str1, "NO")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}
		
		// DC offset (percent of the dynamic range, -50 to 50)
		if (!strcmp(str, "DC_OFFSET")) {
            float dc;
			read = fscanf(f_ini, "%f", &dc);
			if (tr != -1) {
// 				WDcfg->FTDCoffset[tr] = dc;
 				WDcfg->FTDCoffset[tr*2] = (uint32_t)dc;
 				WDcfg->FTDCoffset[tr*2+1] = (uint32_t)dc;
				continue;
			}
            val = (int)((dc+50) * 65535 / 100);
            if (ch == -1)
                for(i=0; i<MAX_SET; i++)
                    WDcfg->DCoffset[i] = val;
            else
                WDcfg->DCoffset[ch] = val;
			continue;
		}
		
		if (strstr(str, "GRP_CH_DC_OFFSET")!=NULL) {
            float dc[8];
			read = fscanf(f_ini, "%f,%f,%f,%f,%f,%f,%f,%f", &dc[0], &dc[1], &dc[2], &dc[3], &dc[4], &dc[5], &dc[6], &dc[7]);
            for(i=0; i<MAX_SET; i++) {
				val = (int)((dc[i]+50) * 65535 / 100); 
				WDcfg->DCoffsetGrpCh[ch][i] = val;
            }
			continue;
		}

		// Threshold
		if (strstr(str, "TRIGGER_THRESHOLD")!=NULL) {
			read = fscanf(f_ini, "%d", &val);
			if (tr != -1) {
//				WDcfg->FTThreshold[tr] = val;
 				WDcfg->FTThreshold[tr*2] = val;
 				WDcfg->FTThreshold[tr*2+1] = val;

				continue;
			}
            if (ch == -1)
                for(i=0; i<MAX_SET; i++)
                    WDcfg->Threshold[i] = val;
            else
                WDcfg->Threshold[ch] = val;
			continue;
		}

		// Group Trigger Enable Mask (hex 8 bit)
		if (strstr(str, "GROUP_TRG_ENABLE_MASK")!=NULL) {
			read = fscanf(f_ini, "%x", &val);
            if (ch == -1)
                for(i=0; i<MAX_SET; i++)
                    WDcfg->GroupTrgEnableMask[i] = val & 0xFF;
            else
                 WDcfg->GroupTrgEnableMask[ch] = val & 0xFF;
			continue;
		}

        // Channel Auto trigger (DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT)
		if (strstr(str, "CHANNEL_TRIGGER")!=NULL) {
            CAEN_DGTZ_TriggerMode_t tm;
			read = fscanf(f_ini, "%s", str1);
            if (strcmp(str1, "DISABLED")==0)
                tm = CAEN_DGTZ_TRGMODE_DISABLED;
            else if (strcmp(str1, "ACQUISITION_ONLY")==0)
                tm = CAEN_DGTZ_TRGMODE_ACQ_ONLY;
			else if (strcmp(str1, "ACQUISITION_AND_TRGOUT")==0)
                tm = CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT;
            else {
                printf("%s: Invalid Parameter\n", str);
                continue;
            }
            if (ch == -1)
                for(i=0; i<MAX_SET; i++)
                    WDcfg->ChannelTriggerMode[i] = tm;
            else
                WDcfg->ChannelTriggerMode[ch] = tm;
		    continue;
		}

        // Front Panel LEMO I/O level (NIM, TTL)
		if (strstr(str, "FPIO_LEVEL")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
			if (strcmp(str1, "TTL")==0)
				WDcfg->FPIOtype = 1;
			else if (strcmp(str1, "NIM")!=0)
				printf("%s: invalid option\n", str);
			continue;
		}

        // Channel Enable (or Group enable for the V1740) (YES/NO)
        if (strstr(str, "ENABLE_INPUT")!=NULL) {
			read = fscanf(f_ini, "%s", str1);
            if (strcmp(str1, "YES")==0) {
                if (ch == -1)
                    WDcfg->EnableMask = 0xFF;
                else
                    WDcfg->EnableMask |= (1 << ch);
			    continue;
            } else if (strcmp(str1, "NO")==0) {
                if (ch == -1)
                    WDcfg->EnableMask = 0x00;
                else
                    WDcfg->EnableMask &= ~(1 << ch);
			    continue;
            } else {
                printf("%s: invalid option\n", str);
            }
			continue;
		}

	
        // Write Channel Enable (for pmtdaq)
        if (strstr(str, "GROUP_WRITE_MASK")!=NULL) {
	  read = fscanf(f_ini, "%x", &val);
	  if (ch == -1)
	    for(i=0; i<MAX_SET; i++)
	      WDcfg->GroupWriteMask[i] = val & 0xFFFF;
	  else
	    WDcfg->GroupWriteMask[ch] = val & 0xFFFF;
	  continue;
	}

	if (strstr(str, "EVENTS_TO_WRITE")!=NULL) {
	  read = fscanf(f_ini, "%d", &WDcfg->EventsToWrite);
	  continue;
	}

	if (strstr(str, "CHANNELS_TO_PLOT")!=NULL) {
	  read = fscanf(f_ini, "%x", &WDcfg->PlotMask);
	  continue;
	}


	if (strstr(str, "SAMPLING_FREQ")!=NULL) {
	  read = fscanf(f_ini, "%f", &WDcfg->Ts); // 1.0=1GHz, 0.4=2.5GHz, 0.2=5GHz
	  continue;
	}
	
        printf("%s: invalid setting\n", str);
	}
	return 0;
}

