#ifndef _COMMON_H_
#define _COMMON_H_

/* **************************************************************************************************** */
/* IRA Istituto di Radioastronomia                                                                      */
/* $Id: Common.h,v 1.1 2011-03-14 14:15:07 a.orlati Exp $																								*/
/*                                                                                                      */
/* This code is under GNU General Public Licence (GPL).                                                 */
/*                                                                                                      */
/* Who                                when           What                                              */
/* Andrea Orlati(aorlati@ira.inaf.it)     15/10/2008      Creation				*/                                 


// if define the data are not sent to the bulk data transfer but printed to the screen. For production version it must be 
// commented out
//#define BKD_DEBUG 

// some hard coded, hardware dependant constants
#define STARTFREQUENCY 0.0 // this is the start frequency of all IF inputs
#define MIN_FREQUENCY 0.0 // min frequency of all IF inputs
#define MAX_FREQUENCY 150.0 // max frequency of all IF inputs
#define BINSNUMBER 1000
#define SAMPLETYPE float
#define SAMPLESIZE sizeof(SAMPLETYPE)  //

#define SENDBUFFERSIZE 100
#define RECBUFFERSIZE 512
#define DATABUFFERSIZE 64000
#define MAX_SECTION_NUMBER 14
#define MAX_BOARDS_NUMBER 14

#define DEFAULT_SAMPLE_RATE 3000.0   //one sample every 40 milliseconds
//#define DEFAULT_SAMPLE_RATE 0.001   //one sample every 1 milliseconds
#define DEFAULT_INTEGRATION 40  // milliseconds of integration
#define MIN_INTEGRATION 5  // milliseconds of integration
#define MAX_INTEGRATION 1000  // milliseconds of integration
#define DEFAULT_DIODE_SWITCH_PERIOD 0 // the cal diode is not switched

#define MAX_BAND_WIDTH 2300.0
#define MIN_BAND_WIDTH 0.0

#define MAX_ATTENUATION 15.0
#define MIN_ATTENUATION 0.0

#define MAX_BINS 16384
#define MIN_BINS 1024

#define MAX_SAMPLE_RATE 4600.0 // MHz

#define FLOW_NUMBER 1
#endif /*COMMON_H_*/
