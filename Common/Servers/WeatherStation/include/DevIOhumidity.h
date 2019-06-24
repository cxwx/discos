#ifndef DevIOHumidity_H_
#define DevIOHumidity_H_

/* **************************************************************************************************** */
/*INAF - OACA                                                                      */
/*  					           */
/*                                                                                                      */
/* This code is under GNU General Public Licence (GPL).                                                 */
/*                                                                                                      */
/* Who                                when            What                                              */
/* Sergio  Poppi(spoppi@oa-cagliari.inaf.it)  	02/01/2009     Creation                                         */

#include <baci.h>
#include <baciDevIO.h>
#include <IRA>
using namespace IRA;

/**
 * This class is derived from template DevIO and it is used by the humidity property of the MeteoStation component
 * @author <a href=mailto:spoppi@oa-cagliari.inaf.it>Sergio Poppi</a>,
 * Istituto di Radioastronomia, Italia<br> 
*/
class DevIOHumidity : public DevIO<CORBA::Double>
{
public:
 	/** 
	 * Constructor
	*/
	DevIOHumidity()
	{		
		AUTO_TRACE("DevIOHumidity::DevIOHumidity()");		
	}

	/**
	 * Destructor
	*/ 
	~DevIOHumidity()
	{
		ACS_TRACE("DevIOHumidity::~DevIOHumidity()");		
	}

	/** 
	 * @return true to initialize the property with default value from CDB.
	*/
	bool initializeValue()
	{		
		AUTO_TRACE("DevIOHumidity::initializeValue()");		
		return false;
	}
	
	/**
	 * Used to read the property value.
	 * @param timestamp epoch when the operation completes
	*/ 
	CORBA::Double read(ACS::Time& timestamp)
	{
		timestamp = getTimeStamp();
        baci::ThreadSyncGuard guard(&m_mutex);
        return m_val;
	}

	/**
	 * It writes values into controller.
     * @param value the new value to set
	 * @param timestamp epoch when the operation completes
	*/ 
	void write(const CORBA::Double& value, ACS::Time& timestamp)
	{
		timestamp = getTimeStamp();
        baci::ThreadSyncGuard guard(&m_mutex);
        m_val = value;
	}
	
private:
	CORBA::Double m_val;
    BACIMutex m_mutex;
};



#endif /*DevIOHumidity_H_*/
