#ifndef NotoWeatherSTATIONIMPLIMPL_H
#define NotoWeatherSTATIONIMPLIMPL_H

#include <baci.h>
#include <baciCharacteristicComponentImpl.h>
// #include "NotoWeatherStationImpl.h"
 


///CORBA generated servant stub
#include <GenericWeatherStationS.h>



///Includes for each BACI property used in this example
#include <baciROdouble.h>
#include <baciRWdouble.h>

//*******DEVIOs*********************************************
#include "DevIOtemperature.h"
#include "DevIOwinddir.h"
#include "DevIOpressure.h"
#include "DevIOwindspeed.h"
#include "DevIOhumidity.h"
#include "DevIOwindspeedPeak.h"

#include <SP_parser.h>
///Include the smart pointer for properties
#include <baciSmartPropertyPointer.h>
#include <ComponentErrors.h>
#include "acstime.h"
 #include <IRA>
#include <cstdlib>
#include "MeteoSocket.h"
#include "MeteoData.h"
#include <ManagementErrors.h>

#include <acsThread.h>


/*
const CString ADDRESS="192.167.8.13"; //DEBUG
const int THREADSLEEPSECONDS=10000; // 1 s in unit of 100ns, from thread sleep time*/


using namespace maci;
using namespace baci; 
using namespace std;

using namespace IRA;

 
/** 
 * @mainpage Noto Meteo Station
 * @date 04/12/2014
 * @version 1.1.0
 * @author <a href=mailto:spoppi@oa-cagliari.inaf.it>Sergio Poppi</a>
 * @remarks Last compiled under ACS 8.2
 * @remarks compiler version is 4.1.2
*/


 /**
The class implements the Noto Weather Station Station.
Not all the paramters from the station have been implemented.


@todo add all the info returned by the HW
*/

using ACS::ThreadBase;

class NotoWeatherStationImpl;

class CMeteoParamUpdaterThread : public ACS::Thread
{
public:
    CMeteoParamUpdaterThread (const ACE_CString& name,
                              NotoWeatherStationImpl* weatherStation,
                              const ACS::TimeInterval& responseTime=ThreadBase::defaultResponseTime,
                              const ACS::TimeInterval& sleepTime=ThreadBase::defaultSleepTime);

    ~CMeteoParamUpdaterThread() { ACS_TRACE("CMeteoParamUpdaterThread::~CMeteoParamUpdaterThread"); }
    virtual void onStop();
    virtual void onStart();
    virtual void runLoop();

private:
    int loopCounter_m;
    int m_threshold; // wind threshold in km/h
    NotoWeatherStationImpl* m_weatherstation_p;
};


class  NotoWeatherStationImpl:     public virtual CharacteristicComponentImpl,     //Standard component superclass  
			 public virtual POA_Weather::GenericWeatherStation
{
public:
	/**
	Constructor 
	
	* @param CompName component's name. This is also the name that will be used to find the configuration data for the component in the Configuration Database.
	* @param containerServices pointer to the class that exposes all services offered by container
	
	
	*/
	NotoWeatherStationImpl(const ACE_CString &CompName,
			maci::ContainerServices *containerServices);


	/**
	 * Destructor.
	*/
	virtual ~NotoWeatherStationImpl(); 


	virtual Weather::parameters getData();
	virtual void updateData() throw (ACSErr::ACSbaseEx, CORBA::SystemException);

	virtual CORBA::Double getTemperature() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
	virtual	CORBA::Double getPressure() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
	virtual	CORBA::Double getHumidity() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
	//virtual	CORBA::Double getWindspeed() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
        virtual CORBA::Double getWindspeedPeak() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
        virtual CORBA::Double getWindSpeedAverage() throw (ACSErr::ACSbaseEx, CORBA::SystemException);
        virtual CORBA::Double getWindDir() throw (ACSErr::ACSbaseEx, CORBA::SystemException);

	virtual CORBA::Boolean command(const char *cmd,CORBA::String_out answer) throw (CORBA::SystemException);


	virtual void initialize() throw (ACSErr::ACSbaseExImpl);
	virtual void	 cleanUp()throw (ACSErr::ACSbaseExImpl);
	virtual void aboutToAbort() throw (ACSErr::ACSbaseExImpl);
        virtual void execute() 	throw (ACSErr::ACSbaseExImpl);

	/** 
	 * Returns a reference to the air temperature  property Implementation of IDL interface.
	 * @return pointer to read-write double temperature
	*/
	 virtual ACS::RWdouble_ptr  temperature ()
	throw (CORBA::SystemException);

	 /** 
	 * Returns a reference to the wind direction  property Implementation of IDL interface.
	 * @return pointer to read-write double wind direction
	*/
 	 virtual ACS::RWdouble_ptr  winddir ()
	throw (CORBA::SystemException);
	 /** 
	 * Returns a reference to the wind speed  property Implementation of IDL interface.
	 * @return pointer to read-write double wind speed
	*/

 	 virtual ACS::RWdouble_ptr  windspeed ()
	throw (CORBA::SystemException);
         virtual ACS::RWdouble_ptr  windspeedpeak ()
        throw (CORBA::SystemException);

	 /** 
	 * Returns a reference to the relative humidity  property Implementation of IDL interface.
	 * @return pointer to read-write double wind speed
	*/

  	 virtual ACS::RWdouble_ptr  humidity ()
	throw (CORBA::SystemException);
	 /** 
	 * Returns a reference to the wind speed  property Implementation of IDL interface.
	 * @return pointer to read-write double relative humidity
	*/
  	 virtual ACS::RWdouble_ptr  pressure ()
	throw (CORBA::SystemException);

	 /** 
	 * Returns a reference to the autoparkThreshold  property Implementation of IDL interface.
	 * @return pointer to read-only double relative humidity
	*/
  	 virtual ACS::ROdouble_ptr  autoparkThreshold()
	throw (CORBA::SystemException);

	/**
	*return  the all the data parameters formatted into a 
	*meteoparameters structure
	* @return  meteoparameters
	*/
 
private:
    void deleteAll();
    CSecureArea<MeteoSocket> *m_socket;
    CError err;
    CString ADDRESS;
    unsigned int PORT;
    SmartPropertyPointer<RWdouble> m_temperature;
    SmartPropertyPointer<RWdouble> m_winddir;
    SmartPropertyPointer<RWdouble> m_windspeed;
    SmartPropertyPointer<RWdouble> m_humidity;
    SmartPropertyPointer<RWdouble> m_pressure;
    SmartPropertyPointer<RWdouble> m_windspeedPeak;
    SmartPropertyPointer<ROdouble> m_autoparkThreshold;
    CMeteoParamUpdaterThread *m_controlThread_p;
    SimpleParser::CParser<MeteoSocket> * m_parser;
    void operator=(const NotoWeatherStationImpl&);
    Weather::parameters m_parameters;
    BACIMutex m_meteoParametersMutex;
};


#endif
