#include "NotoWeatherStationImpl.h"

using namespace baci;

using namespace SimpleParser;



NotoWeatherStationImpl::NotoWeatherStationImpl(
				const ACE_CString &name,
			     maci::ContainerServices * containerServices) :
    CharacteristicComponentImpl(name, containerServices),m_temperature(this),
		       m_winddir(this),
		       m_windspeed(this),
		       m_humidity(this),
		       m_pressure(this),
             m_windspeedPeak(this)
{	
	AUTO_TRACE("NotoWeatherStationImpl::NotoWeatherStationImpl");
}

NotoWeatherStationImpl::~NotoWeatherStationImpl()
{

        /*if(m_controlThread_p!=0)
        getContainerServices()->getThreadManager()->destroy(m_controlThread_p);*/

        AUTO_TRACE("NotoWeatherStationImpl::~NotoWeatherStationImpl");
//	deleteAll();
}


void NotoWeatherStationImpl::cleanUp() throw (ACSErr::ACSbaseExImpl)
{
	if (m_controlThread_p!=NULL) {
		m_controlThread_p->suspend();
		getContainerServices()->getThreadManager()->destroy(m_controlThread_p);
	}
	CharacteristicComponentImpl::cleanUp();	
   AUTO_TRACE("NotoWeatherStationImpl::cleanUp()");
}

void NotoWeatherStationImpl::aboutToAbort() throw (ACSErr::ACSbaseExImpl)
{
	if (m_controlThread_p!=NULL) {
		m_controlThread_p->suspend();
		getContainerServices()->getThreadManager()->destroy(m_controlThread_p);
	}	
	AUTO_TRACE("NotoWeatherStationImpl::aboutToAbort()");
}

void NotoWeatherStationImpl::execute() throw (ACSErr::ACSbaseExImpl)
{
	ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::execute()",(LM_INFO,"COMPSTATE_OPERATIONAL"));
}

void NotoWeatherStationImpl::initialize() throw (ACSErr::ACSbaseExImpl)
{

	MeteoSocket *sock;
	try {
		  if 			(CIRATools::getDBValue(getContainerServices(),"IPAddress",ADDRESS) && CIRATools::getDBValue(getContainerServices(),"port",PORT))
		{
			ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::initialize()",(LM_INFO,"IP address %s, Port %d ",(const char *) ADDRESS,PORT));


		} else
		{
			 ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::initialize()",(LM_ERROR,"Error getting IP address from CDB" ));
		_EXCPT(ComponentErrors::MemoryAllocationExImpl,dummy,"NotoWeatherStationImpl::initialize()");
		throw dummy;
		}    


		
		sock=new MeteoSocket(ADDRESS,PORT);
 		m_socket =new CSecureArea<MeteoSocket>(sock);
		m_temperature=new RWdouble(getContainerServices()->getName()+":temperature", getComponent(), new DevIOTemperature(m_socket),true);
		
		m_winddir=new RWdouble(getContainerServices()->getName()+":winddir", getComponent(), new DevIOWinddir(m_socket),true);
		m_windspeed=new RWdouble(getContainerServices()->getName()+":windspeed", getComponent(), new DevIOWindspeed(m_socket),true);
		m_humidity=new RWdouble(getContainerServices()->getName()+":humidity", getComponent(), new DevIOHumidity(m_socket),true);
		m_pressure=new RWdouble(getContainerServices()->getName()+":pressure", getComponent(), new DevIOPressure(m_socket),true);
                m_windspeedPeak=new RWdouble(getContainerServices()->getName()+":windspeedpeak", getComponent(), new DevIOWindspeedPeak(m_socket),true);
                
                m_controlThread_p = getContainerServices()->getThreadManager()->create<CMeteoParamUpdaterThread, MeteoSocket*>("MeteoParam Updater",sock );
                m_controlThread_p->setSleepTime  (5*10000000);
//              m_controlThread_p->setResponseTime(60*1000000);
                m_controlThread_p->resume();
                
		m_parser=new CParser<MeteoSocket>(sock,10);
		m_parser->add("getWindSpeed",new function0<MeteoSocket,non_constant,double_type >(sock,&MeteoSocket::getWindSpeed),0 );
		m_parser->add("getTemperature",new function0<MeteoSocket,non_constant,double_type >(sock,&MeteoSocket::getTemperature),0 );
		m_parser->add("getHumidity",new function0<MeteoSocket,non_constant,double_type >(sock,&MeteoSocket::getHumidity),0 );
		m_parser->add("getPressure",new function0<MeteoSocket,non_constant,double_type >(sock,&MeteoSocket::getPressure),0 );
	}
	catch (std::bad_alloc& ex) {
		_EXCPT(ComponentErrors::MemoryAllocationExImpl,dummy,"NotoWeatherStationImpl::NotoWeatherStationImpl()");
		throw dummy;
	}catch (ComponentErrors::ComponentErrorsExImpl& E) {
		E.log(LM_DEBUG);
		throw E;
	}

#if 0


	try {
		CSecAreaResourceWrapper<MeteoSocket> sock=m_socket->Get();
		if (!sock->isConnected())
		{
		sock->connection();
// 		cout << "Connected  to Meteo Station @"<<ADDRESS<<":"<<PORT <<endl;
		ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::Disconnect()",(LM_DEBUG,"Connected  to Meteo Station @%s:%d  ",(const char *) ADDRESS,PORT));

		}

	} catch (ComponentErrors::SocketErrorExImpl &x)
	
        {
                 ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::initialize()",
                       (LM_ERROR,"Can not connect  to WeatherStation @%s:%d  ",
                       (const char *) ADDRESS,PORT));

		_THROW_EXCPT(ComponentErrors::SocketErrorExImpl,"NotoWeatherStationImpl::initialize()");
 
	}
#endif 


        AUTO_TRACE("NotoWeatherStationImpl::initialize");
}



CORBA::Boolean NotoWeatherStationImpl::command(const char *cmd,CORBA::String_out answer) throw (CORBA::SystemException)
{
	IRA::CString out;
	bool res;
	CSecAreaResourceWrapper<MeteoSocket> line=m_socket->Get();
	try {
		m_parser->run(cmd,out);
		res=true;
	}
	catch (ParserErrors::ParserErrorsExImpl &ex) {
		res=false;
	}
	catch (ACSErr::ACSbaseExImpl& ex) {
		ex.log(LM_ERROR); // the errors resulting from the execution are logged here as stated in the documentation of CommandInterpreter interface, while the parser errors are never logged.
		res=false;
	}
	answer=CORBA::string_dup((const char *)out);
	return res;
}


void NotoWeatherStationImpl::deleteAll()
{
        AUTO_TRACE("NotoWeatherStationImpl::deleteAll");
   	CError err;
               
 	try{
 	CSecAreaResourceWrapper<MeteoSocket> sock=m_socket->Get();
                
                if (sock->isConnected())
 		{
			sock->disconnection();
 			delete m_socket;
 		}
 	} catch (...)
 	{
 		cout << "unknown exception in closing component " << endl;


 	}
 
	ACS_LOG(LM_FULL_INFO,"NotoWeatherStationImpl::deleteAll()",(LM_DEBUG,"Disconnecting from socket @%s  ",(const char *)err.getFullDescription()));
	//	delete m_socket;
 }


CORBA::Double NotoWeatherStationImpl::getTemperature() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getTemperature");

	double temperature;
	ACSErr::Completion_var completion;
	temperature = m_temperature->get_sync(completion.out());

	return temperature;
}

CORBA::Double NotoWeatherStationImpl::getWindspeed() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getTemperature");

	double windspeed;
	ACSErr::Completion_var completion;
	windspeed = m_windspeed->get_sync(completion.out());

	return windspeed;

}
          
CORBA::Double NotoWeatherStationImpl::getWindDir() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getTemperature");

        double winddir;
        ACSErr::Completion_var completion;
        winddir = m_winddir->get_sync(completion.out());
        return winddir;
        
}
             
          
CORBA::Double NotoWeatherStationImpl::getWindSpeedAverage() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getTemperature");

        double windspeed;
        ACSErr::Completion_var completion;
        windspeed = m_windspeed->get_sync(completion.out());

        return windspeed;

}
          
          
          
CORBA::Double NotoWeatherStationImpl::getWindspeedPeak() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getWindspeedPeak");
        double windspeed;
        ACSErr::Completion_var completion;
        windspeed = m_windspeed->get_sync(completion.out());

        return windspeed;

}
          
          
          
          
CORBA::Double NotoWeatherStationImpl::getPressure() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getpressure()");

	double pressure;
	ACSErr::Completion_var completion;
	pressure = m_pressure->get_sync(completion.out());
	return pressure;



}
CORBA::Double NotoWeatherStationImpl::getHumidity() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
        AUTO_TRACE("NotoWeatherStationImpl::getHumidity()");

	double humidity;
	ACSErr::Completion_var completion;
	humidity = m_humidity->get_sync(completion.out());
	return humidity;



}

Weather::parameters NotoWeatherStationImpl::getData() throw (ACSErr::ACSbaseEx, CORBA::SystemException)
{
	Weather::parameters mp;
        AUTO_TRACE("NotoWeatherStationImpl::getData");
	CError err;
	CString rdata;
	CSecAreaResourceWrapper<MeteoSocket> sock=m_socket->Get();
	double temperature;
	double winddir;
	double windspeed;
	double pressure;
	double humidity;
	
//	sock->sendCMD(err,CString("wx\n"));
//	sock->receiveData(err,rdata);
        

	
	ACSErr::Completion_var completion;

	temperature = m_temperature->get_sync(completion.out());
	winddir     = m_winddir->get_sync(completion.out());
	windspeed   = m_windspeed->get_sync(completion.out());
	pressure    = m_pressure->get_sync(completion.out());
	humidity    = m_humidity->get_sync(completion.out());
	

//  	cout <<"received"<< len << endl;
//	string ss;
//	string srecv;
//	srecv=string((const char*)rdata);
//	vector<string> vrecv;
//
//	istringstream  ist(srecv); // string stream
//	while (ist >> ss) vrecv.push_back(ss) ;// split the string
//	int ndata=vrecv.size();
//	if (ndata > 3)
//
//	{
// 		temperature = atof(vrecv[ndata-3].c_str());
//		pressure = atof(vrecv[ndata-2].c_str());
//		humidity  = atof(vrecv[ndata-1].c_str());
//
////		windspeed  = atof(vrecv[ndata-1].c_str());
//
//	} else
//	{
//		 ACS_LOG(LM_FULL_INFO,"MeteoSocket::update()",(LM_ERROR,"Not enough data from meteo server"));
//	}
	

	mp.temperature=temperature;
	mp.humidity   =humidity;
 	mp.wind       =windspeed;
	mp.pressure   =pressure;
	

	return mp;
}
 

//pdate();



ACS::RWdouble_ptr
NotoWeatherStationImpl::temperature ()
    throw (CORBA::SystemException)
{
    if (m_temperature == 0)
	{
	return ACS::RWdouble::_nil();
	}
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_temperature->getCORBAReference());
    return prop._retn();
}

ACS::RWdouble_ptr
NotoWeatherStationImpl::winddir ()
    throw (CORBA::SystemException)
{
    if (m_winddir == 0)
	{
	return ACS::RWdouble::_nil();
	}
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_winddir->getCORBAReference());
    return prop._retn();
}

ACS::RWdouble_ptr
NotoWeatherStationImpl::windspeed ()
    throw (CORBA::SystemException)
{
    if (m_windspeed == 0)
	{
	return ACS::RWdouble::_nil();
	}
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_windspeed->getCORBAReference());
    return prop._retn();
}



ACS::RWdouble_ptr
NotoWeatherStationImpl::humidity ()
    throw (CORBA::SystemException)
{
    if (m_humidity == 0)
	{
	return ACS::RWdouble::_nil();
	}
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_humidity->getCORBAReference());
    return prop._retn();
}
ACS::RWdouble_ptr
NotoWeatherStationImpl::pressure ()
    throw (CORBA::SystemException)
{
    if (m_pressure == 0)
	{
	return ACS::RWdouble::_nil();
	}
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_pressure->getCORBAReference());
    return prop._retn();
}

ACS::RWdouble_ptr
NotoWeatherStationImpl::windspeedpeak ()
    throw (CORBA::SystemException)
{
    if (m_windspeedPeak == 0)
        {
        return ACS::RWdouble::_nil();
        }
    
    ACS::RWdouble_var prop = ACS::RWdouble::_narrow(m_windspeedPeak->getCORBAReference());
    return prop._retn();
}



/* --------------- [ MACI DLL support functions ] -----------------*/
#include <maciACSComponentDefines.h>
MACI_DLL_SUPPORT_FUNCTIONS(NotoWeatherStationImpl)
/* ----------------------------------------------------------------*/


