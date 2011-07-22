// $Id: Configuration.cpp,v 1.2 2011-05-12 14:14:31 a.orlati Exp $

#include "Configuration.h"

#define _GET_DWORD_ATTRIBUTE(ATTRIB,DESCR,FIELD) { \
	DWORD tmpw; \
	if (!CIRATools::getDBValue(Services,ATTRIB,tmpw)) { \
		_EXCPT(ComponentErrors::CDBAccessExImpl,dummy,"CConfiguration::Init()"); \
		dummy.setFieldName(ATTRIB); \
		throw dummy; \
	} \
	else { \
		FIELD=tmpw; \
		ACS_DEBUG_PARAM("CConfiguration::Init()",DESCR" %lu",tmpw); \
	} \
}

#define _GET_STRING_ATTRIBUTE(ATTRIB,DESCR,FIELD) { \
	CString tmps; \
	if (!CIRATools::getDBValue(Services,ATTRIB,tmps)) { \
		_EXCPT(ComponentErrors::CDBAccessExImpl,dummy,"::CConfiguration::Init()"); \
		dummy.setFieldName(ATTRIB); \
		throw dummy; \
	} \
	else { \
		FIELD=tmps; \
		ACS_DEBUG_PARAM("CConfiguration::Init()",DESCR" %s",(const char*)tmps); \
	} \
}

CConfiguration::CConfiguration()
{
	m_configurationTable=NULL;
}

CConfiguration::~CConfiguration()
{
	if (m_configurationTable!=NULL) {
		m_configurationTable->closeTable();
		delete m_configurationTable;
	}
}

bool CConfiguration::getSetupFromID(const IRA::CString setupID,TBackendSetup& setup)
{
	bool done=false;
	m_configurationTable->First();
	for (int i=0;i<m_configurationTable->recordCount();i++) {
		CString id((*m_configurationTable)["configurationID"]->asString());
		if (setupID==id) {
			done=true;
			break;
		}
		m_configurationTable->Next();
	}
	if (done) {
		setup.setupID=(*m_configurationTable)["configurationID"]->asString();
		setup.sections=(*m_configurationTable)["sections"]->asLongLong();
		if (setup.sections>MAX_INPUT_NUMBER) {
			return false;
		}
		if ((*m_configurationTable)["calSwitchingEnabled"]->asString()=="TRUE") {
			setup.calSwitchEnabled=true;
		}
		else {
			setup.calSwitchEnabled=false;
		}
		if ((*m_configurationTable)["inputPort"]->asString()=="PRIMARY") {
			setup.inputPort=CProtocol::PRIMARY;
		}
		else if ((*m_configurationTable)["inputPort"]->asString()=="BWG") {
			setup.inputPort=CProtocol::BWG;
		}
		else { // GREGORIAN
			setup.inputPort=CProtocol::GREGORIAN;
		}
		setup.beams=(*m_configurationTable)["beams"]->asLongLong();
		setup.bandWidth=(*m_configurationTable)["bandWidth"]->asDouble();
		setup.attenuation=(*m_configurationTable)["attenuation"]->asDouble();
		int start_si=0,start_pol=0,start_feed=0,start_ifs=0;
		IRA::CString sections_inputs((*m_configurationTable)["sections_inputs"]->asString());
		IRA::CString polarizations((*m_configurationTable)["polarizations"]->asString());
		IRA::CString feed((*m_configurationTable)["feed"]->asString());
		IRA::CString ifs((*m_configurationTable)["IF"]->asString());
		IRA::CString ret;
		for(int i=0;i<setup.sections;i++) {
			if (!IRA::CIRATools::getNextToken(sections_inputs,start_si,' ',ret)) return false;
			setup.sections_inputs[i]=ret.ToLong();
			if (!IRA::CIRATools::getNextToken(polarizations,start_pol,' ',ret)) return false;
			if (ret=="L") {
				setup.polarizations[i]=Backends::BKND_LCP;
			}
			else {
				setup.polarizations[i]=Backends::BKND_RCP;
			}
			if (!IRA::CIRATools::getNextToken(feed,start_feed,' ',ret)) return false;
			setup.feed[i]=ret.ToLong();
			if (!IRA::CIRATools::getNextToken(ifs,start_ifs,' ',ret)) return false;
			setup.ifs[i]=ret.ToLong();
			ACS_LOG(LM_FULL_INFO,"CConfiguration::getSetupFromID()",(LM_DEBUG,"Sections: %d - Polarization: %d - Feed: %d - Ifs: %d",setup.sections_inputs[i],setup.polarizations[i],setup.feed[i],setup.ifs[i]));	
		}
		return true;
	}
	else {
		return false;
	}
}

void CConfiguration::init(maci::ContainerServices *Services) throw (ComponentErrors::CDBAccessExImpl,ComponentErrors::MemoryAllocationExImpl,ComponentErrors::IRALibraryResourceExImpl)
{
	_GET_STRING_ATTRIBUTE("IPAddress","TCP/IP address is: ",m_sAddress);
	_GET_DWORD_ATTRIBUTE("Port","Port is: ",m_wPort);
	_GET_DWORD_ATTRIBUTE("CommandLineTimeout","Time out of the command line is (uSec): ",m_dwCommandLineTimeout);
	_GET_DWORD_ATTRIBUTE("ConnectTimeout","Time out of the command line connection  is (uSec): ",m_dwConnectTimeout);
	_GET_DWORD_ATTRIBUTE("PropertyRefreshTime","The property refresh time is (microseconds): ",m_dwPropertyRefreshTime);
	_GET_STRING_ATTRIBUTE("Configuration","Backend configuration is: ",m_sConfig);
	_GET_DWORD_ATTRIBUTE("TimeTollerance","The time tollerance is  (microseconds): ",m_dwTimeTollerance);	
	_GET_DWORD_ATTRIBUTE("RepetitionCacheTime","Log repetition filter cache time (uSec)",m_dwRepetitionCacheTime);
	_GET_DWORD_ATTRIBUTE("RepetitionExpireTime","Log repetition filter expire time  (uSec)",m_dwRepetitionExpireTime);
	_GET_DWORD_ATTRIBUTE("DataPort","Port of data line is: ",m_wDataPort);
	_GET_STRING_ATTRIBUTE("DataIPAddress","Data line TCP/IP address is: ",	m_sDataAddress);
	_GET_DWORD_ATTRIBUTE("DataLatency","The latency of the backend data line is (uSec)",m_dwDataLatency);
	_GET_DWORD_ATTRIBUTE("SenderSleepTime","The sender thread sleep time is (uSec)",m_dwSenderSleepTime);
	_GET_DWORD_ATTRIBUTE("SenderResponseTime","The sender thread response time is (uSec)",m_dwSenderResponseTime);
	_GET_DWORD_ATTRIBUTE("ControlSleepTime","The control thread sleep time is (uSec)",m_dwControlSleepTime);
	_GET_DWORD_ATTRIBUTE("ControlResponseTime","The control thread response time is (uSec)",m_dwControlResponseTime);
	_GET_DWORD_ATTRIBUTE("DeviceNumber","The number of installed boards is ",m_dwDeviceNumber);
	_GET_DWORD_ATTRIBUTE("DataBufferSize","Size of data packet is (bytes)",m_dwDataBufferSize);
	
	// read the configurations
	IRA::CError error;
	try {
		m_configurationTable=(CDBTable *)new CDBTable(Services,"Configuration","DataBlock/TotalPower");
	}
	catch (std::bad_alloc& ex) {
		_EXCPT(ComponentErrors::MemoryAllocationExImpl,dummy,"TotalPowerImpl::initialize()");
		throw dummy;
	}
	error.Reset();
	if (!m_configurationTable->addField(error,"configurationID", CDataField::STRING)) {
		error.setExtra("Error adding field configurationID", 0);
	}
	else if (!m_configurationTable->addField(error,"sections",CDataField::LONGLONG)) {
		error.setExtra("Error adding field sections", 0);
	}
	else if (!m_configurationTable->addField(error,"calSwitchingEnabled",CDataField::STRING)) {
		error.setExtra("Error adding field calSwitchingEnabled", 0);
	}
	else if (!m_configurationTable->addField(error,"inputPort", CDataField::STRING)) {
		error.setExtra("Error adding field inputPort", 0);
	}
	else if (!m_configurationTable->addField(error,"beams", CDataField::LONGLONG))	{
		error.setExtra("Error adding field beams", 0);
	}
	else if (!m_configurationTable->addField(error,"bandWidth",CDataField::DOUBLE)) {
		error.setExtra("Error adding field bandWidth", 0);
	}
	else if (!m_configurationTable->addField(error,"attenuation",CDataField::DOUBLE)) {
		error.setExtra("Error adding field attenuation", 0);
	}
	else if (!m_configurationTable->addField(error,"sections_inputs", CDataField::STRING)) {
		error.setExtra("Error adding field sections_inputs", 0);
	}
	else if (!m_configurationTable->addField(error,"polarizations", CDataField::STRING))	{
		error.setExtra("Error adding field polarizations", 0);
	}
	else if (!m_configurationTable->addField(error,"feed",CDataField::STRING)) {
		error.setExtra("Error adding field feed", 0);
	}
	else if (!m_configurationTable->addField(error,"IF",CDataField::STRING)) {
		error.setExtra("Error adding field IF", 0);
	}
	if (!error.isNoError()) {
		_EXCPT_FROM_ERROR(ComponentErrors::IRALibraryResourceExImpl,dummy,error);
		dummy.setCode(error.getErrorCode());
		dummy.setDescription((const char*)error.getDescription());
		throw dummy;
	}
	if (!m_configurationTable->openTable(error))	{
		_EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
		throw dummy;
	}
	
}

