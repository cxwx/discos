
#include "Protocol.h"
#include "Common.h"

WORD CProtocol::setConfiguration_broadcast(char *sBuff,const TInputs& input,const double &att,const double& bw)
{
	strcpy(sBuff,PROT_SET_CONF_BROAD);
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeInput(input));
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeAttenuationLevel(att));
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeBandWidth(bw));
	strcat(sBuff,PROT_TERMINATOR);
	return strlen (sBuff);
}

WORD CProtocol::setConfiguration(char *sBuff,long id,const TInputs& input,const double &att,const double& bw,long *boards)
{
	char tempBuff[128];
	strcpy(sBuff,PROT_SET_CONF);
	strcat(sBuff,PROT_SEPARATOR);
	if (boards==NULL) {
		sprintf(tempBuff,"%ld",id+1); //the inner enumeration goes from 0 to SECTIONNUMBER-1,
														// in the backend goes from 1 to BOARDNUMBER
	}
	else {
		sprintf(tempBuff,"%ld",boards[id]+1); //the inner enumeration goes from 0 to SECTIONNUMBER-1,
	}
	strcat(sBuff,tempBuff);
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeInput(input));
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeAttenuationLevel(att));
	strcat(sBuff,PROT_SEPARATOR);
	strcat(sBuff,(const char *)encodeBandWidth(bw));
	strcat(sBuff,PROT_TERMINATOR);
	return strlen (sBuff);	
}

WORD CProtocol::askBackendConfiguration(char *sBuff)
{
	strcpy(sBuff,PROT_GET_CONF);
	strcat(sBuff,PROT_TERMINATOR);
	return strlen(sBuff);
}

WORD CProtocol::checkBackendTime(char *sBuff)
{
	TIMEVALUE now,tm;
	TIMEDIFFERENCE diff;
	long long seconds;
	char tempBuff[128];
	strcpy(sBuff,PROT_CHECK_TIME);
	strcat(sBuff,PROT_SEPARATOR);
	getReferenceTime(tm);
	IRA::CIRATools::getTime(now);
	diff.value(now.difference(tm.value()));
	seconds=diff.day()*86400+diff.hour()*3600+diff.minute()*60+diff.second();
	sprintf(tempBuff,"%lld%c%ld",seconds,PROT_SEPARATOR_CH,(long)diff.microSecond());
	strcat(sBuff,tempBuff);
	strcat(sBuff,PROT_TERMINATOR);
	return strlen(sBuff);	
}

WORD CProtocol::setBackendTime(char *sBuff)
{
	TIMEVALUE now,tm;
	TIMEDIFFERENCE diff;
	long micro;
	long long seconds;
	char tempBuff[128];
	strcpy(sBuff,PROT_SET_TIME);
	strcat(sBuff,PROT_SEPARATOR);
	getReferenceTime(tm);
	IRA::CIRATools::getTime(now);
	micro=now.microSecond();
	if (micro<300000) {
		micro=(400000-micro); // take 1 more tens of second. 
		IRA::CIRATools::Wait(0,micro);
		IRA::CIRATools::getTime(now);
	}
	else if (micro>700000) {
		micro=1400000-micro; // take 1 more tens of second.
		IRA::CIRATools::Wait(0,micro);
		IRA::CIRATools::getTime(now);
	}
	diff.value(now.difference(tm.value()));
	seconds=diff.day()*86400+diff.hour()*3600+diff.minute()*60+diff.second();
	sprintf(tempBuff,"%lld%c%ld",seconds,PROT_SEPARATOR_CH,(long)diff.microSecond());
	strcat(sBuff,tempBuff);
	strcat(sBuff,PROT_TERMINATOR);
	return strlen(sBuff);
}

 WORD CProtocol::startAcquisition(char *sBuff,const double& sampleRate,const long& cal,const WORD& port,const IRA::CString& address)
 {
	 double sR;
	 char tempBuff[200];
	 strcpy(sBuff, PROT_START_ACQ);
	 strcat(sBuff,PROT_SEPARATOR);
	 sR=sampleRate*1000.0; //sample rate in KHz
	 sR=1/sR; // sample rate given in milliseconds (sample period)
	 sprintf(tempBuff,"%ld%c%ld%c%hd%c%s%c%hd",(long)round(sR),PROT_SEPARATOR_CH,cal,PROT_SEPARATOR_CH,
			 0,PROT_SEPARATOR_CH,(const char*)address,PROT_SEPARATOR_CH,port); 
	 strcat(sBuff,tempBuff);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);
 }
 
 WORD CProtocol::setIntegrationTime(char *sBuff,const long& integrationTime)
 {
	 char tempBuff[200];
	 strcpy(sBuff, PROT_SET_SAMPLE_RATE);
	 strcat(sBuff,PROT_SEPARATOR);
	 sprintf(tempBuff,"%ld%c",integrationTime,PROT_SEPARATOR_CH); 
	 strcat(sBuff,tempBuff);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);	 
 }
 
WORD CProtocol::AutoGainControl(char *sBuff,const WORD& channels,const WORD& level)
{
	 char tempBuff[200];
	 strcpy(sBuff,PROT_AUTO_GAIN);
	 strcat(sBuff,PROT_SEPARATOR);
	 sprintf(tempBuff,"%d%c%d%c",channels,PROT_SEPARATOR_CH,level,PROT_SEPARATOR_CH);
	 strcat(sBuff,tempBuff);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);
}
 
WORD CProtocol::getSample(char *sBuff)
{
	 strcpy(sBuff,PROT_GET_SAMPLE);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);	 
}
 
 WORD CProtocol::resumeAcquisition(char *sBuff)
 {
	 strcpy(sBuff,PROT_RESUME_ACQ);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);	 
 }
 
 WORD CProtocol::stopAcquisition(char *sBuff)
{
	 strcpy(sBuff,PROT_STOP_ACQ);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);
}
 
WORD CProtocol::setZeroInput(char *sBuff,bool activate)
{
	 strcpy(sBuff,PROT_SET_ZERO);
	 strcat(sBuff,PROT_SEPARATOR);
	 if (activate) strcat(sBuff,"1");
	 else strcat(sBuff,"0");
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);
}
 
 WORD CProtocol::suspendAcquisition(char *sBuff)
 {
	 strcpy(sBuff,PROT_SUSPEND_ACQ);
	 strcat(sBuff,PROT_TERMINATOR);
	 return strlen(sBuff);
 }

void CProtocol::swap(DWORD* word)
 {
	 DWORD p;
	 p=((*word)<<16)|((*word)>>16);
	 *word=p;
 }

bool CProtocol::status(char *rBuff,double &timestamp, char *statusCode, bool &acquiring)
{

}

bool CProtocol::setConfiguration(const char *rBuff)
{
	IRA::CString str(rBuff);
	IRA::CString ret;
	int start=0;

	// get !set-configuration string
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_COMMA,ret)) {
		return false;
	}
	// get !set-configuration string
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_COMMA,ret)) {
		return false;
	}

}

bool CProtocol::decodeBackendConfiguration(const char *rBuff,const long& sectionNumber,const DWORD& boardsNumber,double *att,double *bw,TInputs *in,TIMEVALUE& tm,long& currentSR,long * boards)
{
	IRA::CString str(rBuff);
	IRA::CString ret;
	int start=0;
	double boardBw[MAX_BOARDS_NUMBER],boardAtt[MAX_BOARDS_NUMBER];
	TInputs boardIn[MAX_BOARDS_NUMBER];
	// get fox cpu time (seconds)...not used yet
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	// get fox cpu time (microseconds)...not used yet
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	// get FPGA time ...just seconds are reported
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	else {
		if (!decodeFPGATime(ret,tm)) {
			return false;
		}
	}
	// get the sample rate in millisecs.......
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	else {	
		currentSR=ret.ToLong();
	}
	// get the cal diode cycle period (multiple of the sample rate).
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	else {
	}
	// get the zero cycle period (multiple of the caldiode cycle)
	if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
		return false;
	}
	else {
	}
	for (DWORD i=0;i<boardsNumber;i++) {
		// read the input source.....
		if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
			return false;
		}
		else {
			if (!decodeInput(ret,boardIn[i])) {
				return false;
			}
		}
		// read the attenuation level.....
		if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
			return false;
		}
		else {
			if (!decodeAttenuationLevel(ret,boardAtt[i])) {
				return false;
			}
		}
		// read the bandwidth
		if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_CH,ret)) {
			return false;
		}
		else {
			if (!decodeBandWidth(ret,boardBw[i])) {
				return false;
			}
		}
	}
	// now I will build the output arrays
	for (long j=0;j<sectionNumber;j++) {
		if (boards==NULL) {
			att[j]=boardAtt[j];
			in[j]=boardIn[j];
			bw[j]=boardBw[j];
		}
		else {
			att[j]=boardAtt[boards[j]];
			in[j]=boardIn[boards[j]];
			bw[j]=boardBw[boards[j]];
		}
	}
	return true;
}

bool CProtocol::decodeBackendTime(const char *rBuff,bool& res)
{
	long hostMicro,foxMicro;
	long long hostSec,fpgaSec,foxSec;
	if (sscanf(rBuff,"%lld, %ld, %lld, %ld, %lld",&hostSec,&hostMicro,&foxSec,&foxMicro,&fpgaSec)!=5) {
		return false;
	}
	else {
		res=(fpgaSec==hostSec);
		return true;
	}
}

bool CProtocol::checkBackendTime(const char * rBuff,const long threshold,bool& res)
{
	long hostMicro,foxMicro;
	long long hostSec,fpgaSec,foxSec;
	long long micro;
	res=true;
	//printf("%s\n",rBuff);
	if (sscanf(rBuff,"%lld, %ld, %lld, %ld, %lld",&hostSec,&hostMicro,&foxSec,&foxMicro,&fpgaSec)!=5) { //both fgpa and fox timer are checkd
		return false; //parsing error
	}
	else { 
		if ((hostMicro>threshold) && (hostMicro<(1000000-threshold))) { //if the time check request has been issued around the middle of a second...
			if ((hostSec!=foxSec) && (hostSec!=fpgaSec)) {
				res=false;
			}
		}
		else {
			if (hostSec>foxSec) {
				micro=(hostSec-foxSec)*1000000;
				micro+=(hostMicro-foxMicro);
			}
			else if (foxSec>hostSec) {
				micro=(foxSec-hostSec)*1000000;
				micro+=(foxMicro-hostMicro);
			}
			else {
				micro=ABS(hostMicro,foxMicro);
			}
			if (micro>threshold) {
				res=false;
			}
		}
	}
	return true;
}

bool CProtocol::isAck(const char *rBuff)
{
	return (strcmp(rBuff,PROT_ACK)==0);
}

bool CProtocol::isTpi(const WORD& flag)
{
	return (((flag & 0x30) >> 4)==0);
}

bool CProtocol::isCal(const WORD& flag)
{
	return (((flag & 0x30) >> 4)==1);
}

bool CProtocol::isNewStream(const WORD& previousStatus,const WORD& previousCounter,const WORD& currentCounter)
{
	if (previousStatus==0) {
		return true;
	}
	else if (ABS(previousCounter,currentCounter)>1) {
		if ((previousCounter==65535) && (currentCounter==0)) {
			return false;
		}
		else return true;
	}
	else return false;
}

bool CProtocol::decodeData(const char *rBuff,DWORD *data,const DWORD& boardsNumber,const long& sectionNumber,long *boards)
{
	IRA::CString str(rBuff);
	IRA::CString ret;
	int start=0;
	DWORD boardOutput[MAX_BOARDS_NUMBER];

    printf("outtpi = %s\n",rBuff);

	// get fox cpu time (seconds)...not used yet
	if (!CIRATools::getNextToken(str,start,PROT_GET_TPI_STR,ret)) {
        printf("uno\n");
		return false;
	}
	// get sample counter....it should be always zero
	//if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_COMMA,ret)) {
    //    printf("due\n");
	//	return false;
	//}
	// get sample flag.....not used for this purpose
	//if (!CIRATools::getNextToken(str,start,PROT_OK,ret)) {
    //    printf("tre\n");
	//	return false;
	//}
	//if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_COMMA,ret)) {
	//	return false;
	//}
	//for (DWORD  j=0;j<boardsNumber;j++) {
	for (DWORD  j=0;j<2;j++) {
		if (!CIRATools::getNextToken(str,start,PROT_SEPARATOR_COMMA,ret)) {
            printf("quattro\n");
			return false;
		}
		boardOutput[j]=ret.ToLong();
        printf("tpi = %d\n",boardOutput[j]);
	}
	//for (long  k=0;k<sectionNumber;k++) {
	for (long  k=0;k<2;k++) {
		//if (boards==NULL) {
			data[k]=boardOutput[k];
		//}
		//else {
		//	data[k]=boardOutput[boards[k]];
		//}
	}
	return true;
}

void CProtocol::decodeData(BYTE *buff,const double& sampleRate,const WORD& prevStatus,const WORD& prevCounter,
		TIMEVALUE& tm,WORD& counter,WORD& flag,DWORD *data,const long& sectionNumber,long *boards)
{
	static DWORD previousClock=0;
	static WORD sampleCounter=0;
	DWORD clock;
	DWORD *boardData;
	//CProtocol::swap((DWORD *)buff);
	counter=*((WORD*)(buff+4));
	flag=*((WORD*)(buff+6));
	clock=*((DWORD*)buff);
	if (clock!=previousClock) { // stores the sample counter when the clock changes
		if (isNewStream(prevStatus,prevCounter,counter) ) {
			sampleCounter=1;   // this is a workaround the samples belonging to the first second of data are 1 less than expected and the missing one is the first
		}
		else {
			sampleCounter=0;
		}
	}
	else {
		sampleCounter++;
	}
	previousClock=clock;
	decodeFPGATime(clock,sampleRate,sampleCounter,tm);
	boardData=(DWORD *)(buff+8);
	for(long k=0;k<sectionNumber;k++) {
		if (boards==NULL) {
			data[k]=boardData[k];
		}
		else {
			data[k]=boardData[boards[k]];
		}
	}
}	

bool CProtocol::decodeAttenuationLevel(const IRA::CString& str,double& val)
{
	val=str.ToDouble();
	return true;
}

IRA::CString CProtocol::encodeAttenuationLevel(const double& att)
{
	IRA::CString temp;
	temp.Format("%d",(int)att);
	return temp;
}

bool CProtocol::decodeInput(const IRA::CString& str,TInputs& val)
{
	if (str=="NA") {
		val=OHM50;
	}
	else if (str=="PRIM") {
		val=PRIMARY;
	}
	else if (str=="BWG") {
		val=BWG;
	}
	else if (str=="GREG") {
		val=GREGORIAN;
	}
	else if (str=="50_OHM") {
		val=OHM50;
	}
	else {
		return false;
	}
	return true;
}

IRA::CString CProtocol::encodeInput(const TInputs& in)
{
	if (in==PRIMARY) {
		return IRA::CString("P");
	}
	else if (in==BWG) {
		return IRA::CString("B");
	}
	else if (in==GREGORIAN) {
		return IRA::CString("G");
	}
	else  { //OHM50
		return IRA::CString("Z");
	}
}

bool CProtocol::decodeBandWidth(const IRA::CString& str,double& val)
{
	long bw;
	if (str=="NA") {
		val=0.0;
	}
	else {
		bw=str.ToLong();
		val=(double)bw;
	}
	return true;
}

IRA::CString CProtocol::encodeBandWidth(const double& bw)
{
	IRA::CString temp;
	int num;
	if (bw<=300.0) {
		num=4;
	}
	else if (bw<=730.0) {
		num=3;
	}
	else if (bw<=1250.0) {
		num=2;
	}
	else {  //2350.0
		num=1;
	}
	temp.Format("%d",num);
	return temp;
}

bool CProtocol::decodeFPGATime(const DWORD& clock,const double& sampleRate,const WORD& counter,TIMEVALUE& tm)
{
	DWORD clockT;
	double period;
	long second,minute,hour,day,micro;
	ACS::TimeInterval temp;
	TIMEDIFFERENCE bkndTime;
	bkndTime.reset(); 
	bkndTime.normalize(true);
	clockT=clock;
	day=clockT/86400;
	clockT%=86400;
	hour=clockT/3600;
	clockT%=3600;
	minute=clockT/60;
	second=clockT % 60;
	period=1.0/sampleRate; // sample rate is MHz so period is in microseconds
	//micro=(long)(period*((double)counter-0.5)); // take the mean time
	micro=(long)(period*((double)counter)); // I'd like to assign the time to the start of the sample!
	bkndTime.day(day);
	bkndTime.hour(hour);
	bkndTime.minute(minute);
	bkndTime.second(second);
	bkndTime.microSecond(micro);
	temp=bkndTime.value().value;
	temp-=(long long)(period*10);   // we have to shift back the time of one period...because the sample is referred to the period before.
	bkndTime.value(temp);
	getReferenceTime(tm);
	tm+=bkndTime.value();
	return true;
}

bool CProtocol::decodeFPGATime(const IRA::CString& str,TIMEVALUE& tm)
{
	long long seconds;
	long second,minute,hour,day;
	TIMEDIFFERENCE bkndTime;
	bkndTime.reset(); 
	bkndTime.normalize(true);
	seconds=str.ToLongLong();
	day=seconds/86400;
	seconds%=86400;
	hour=seconds/3600;
	seconds%=3600;
	minute=seconds/60;
	second=seconds % 60;
	bkndTime.day(day);
	bkndTime.hour(hour);
	bkndTime.minute(minute);
	bkndTime.second(second);
	getReferenceTime(tm);
	tm+=bkndTime.value();
	return true;
}

void CProtocol::getReferenceTime(TIMEVALUE& tm)
{
	tm.reset();
	tm.normalize(true);
	tm.year(1970);
	tm.month(1);
	tm.day(1);
}
