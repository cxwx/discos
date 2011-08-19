#ifndef _DEVIOVACUUM_H_
#define _DEVIOVACUUM_H_

/** **************************************************************************************************** */
/* IRA Istituto di Radioastronomia                                                                      */
/*                                                                                                      */
/* This code is under GNU General Public Licence (GPL).                                                 */
/*                                                                                                      */
/* Who                                when            What                                              */
/* Andrea Orlati(aorlati@ira.inaf.it) 18/08/2011     Creation                                         */


#include <baciDevIO.h>
#include <IRA>

/**
 * This class is derived from template DevIO and it is used by the vacuum  property of the  component.
 * @author <a href=mailto:a.orlati@ira.inaf.it>Andrea Orlati</a>,
 * Istituto di Radioastronomia, Italia<br>
*/
class DevIOVacuum : public DevIO<CORBA::Double>
{
public:

	/**
	 * Constructor
	 * @param core pointer to the component core
	*/
	DevIOVacuum(CComponentCore* core) :  m_pCore(core)
	{
		AUTO_TRACE("DevIOVacuum::DevIOVacuum()");
	}

	/**
	 * Destructor
	*/
	~DevIOVacuum()
	{
		ACS_TRACE("DevIOVacuum::~DevIOVacuum()");
	}

	/**
	 * @return true to initialize the property with default value from CDB.
	*/
	bool initializeValue()
	{
		AUTO_TRACE("DevIOVacuum::DevIOVacuum()");
		return true; // initialize with the default in order to avoid the alarm system when the component start and the value has not been read at least once
	}

	/**
	 * Used to read the property value.
	 * @param timestamp epoch when the operation completes
	*/
	CORBA::Double read(ACS::Time& timestamp) throw (ACSErr::ACSbaseExImpl)
	{
		m_val=m_pCore->getVacuum();
		timestamp=getTimeStamp();  //Completion time
		return m_val;
	}
	/**
	 * It writes values into controller. Unused because the properties are read-only.
	*/
	void write(const CORBA::Double& value, ACS::Time& timestamp) throw (ACSErr::ACSbaseExImpl)
	{
		timestamp=getTimeStamp();
		return;
	}

private:
	CComponentCore* m_pCore;
	double  m_val;
};

#endif
