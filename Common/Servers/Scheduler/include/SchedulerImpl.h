#ifndef SCHEDULERIMPL_H_
#define SCHEDULERIMPL_H_

/* ************************************************************************************************************* */
/* IRA Istituto di Radioastronomia                                                                               */
/* $Id: SchedulerImpl.h,v 1.20 2011-06-21 16:39:52 a.orlati Exp $										         */
/*                                                                                                               */
/* This code is under GNU General Public Licence (GPL).                                                          */
/*                                                                                                               */
/* Who                                                when                     What                                                       */
/* Andrea Orlati(aorlati@ira.inaf.it) 18/12/2008       Creation                                                  */


#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include <baciCharacteristicComponentImpl.h>
#include <baciSmartPropertyPointer.h>
#include <enumpropROImpl.h>
#include <baciROdouble.h>
#include <baciROstring.h>
#include <baciROlong.h>
#include <SchedulerS.h>
#include <ComponentErrors.h>
#include <ManagementErrors.h>
#include <Definitions.h>
#include "Configuration.h"
#include "Core.h"

/** 
 * @mainpage Scheduler component Implementation 
 * @date 13/10/2010
 * @version 1.12.0
 * @author <a href=mailto:a.orlati@ira.inaf.it>Andrea Orlati</a>
 * @remarks Last compiled under ACS 7.0.2
 * @remarks compiler version is 4.1.2
*/

using namespace baci;

/**
 * This class is the implementation of the Scheduler component.  
 * @author <a href=mailto:a.orlati@ira.inaf.it>Orlati Andrea</a>
 * Istituto di Radioastronomia, Italia
 * <br> 
 */
class SchedulerImpl: public CharacteristicComponentImpl,
				       public virtual POA_Management::Scheduler
{

public:
	
	/** 
	* Constructor.
	* @param CompName component's name. This is also the name that will be used to find the configuration data for the component 
	*                in the Configuration Database.
	* @param containerServices pointer to the class that exposes all services offered by container
	*/
	SchedulerImpl(const ACE_CString &CompName,maci::ContainerServices *containerServices);

	/**
	 * Destructor.
	*/
	virtual ~SchedulerImpl(); 

	/** 
	 * Called to give the component time to initialize itself. The component reads in configuration files/parameters or 
	 * builds up connection to devices or other components. 
	 * Called before execute. It is implemented as a synchronous (blocking) call.
	 * @throw ACSErr::ACSbaseExImpl
	*/
	virtual void initialize() throw (ACSErr::ACSbaseExImpl);

	/**
 	 * Called after <i>initialize()</i> to tell the component that it has to be ready to accept incoming 
 	 * functional calls any time. 
	 * Must be implemented as a synchronous (blocking) call. In this class the default implementation only 
	 * logs the COMPSTATE_OPERATIONAL
	 * @throw ACSErr::ACSbaseExImpl
	*/
	virtual void execute() throw (ACSErr::ACSbaseExImpl);
	
	/** 
	 * Called by the container before destroying the server in a normal situation. This function takes charge of 
	 * releasing all resources.
	*/
	virtual void cleanUp();
	
	/** 
	 * Called by the container in case of error or emergency situation. This function tries to free all resources 
	 * even though there is no warranty that the function is completely executed before the component is destroyed.
	*/	
	virtual void aboutToAbort();
	
	/**
     * Returns a reference to the status property implementation of IDL interface.
	 * @return pointer to read-only ROTSystemStatus property status
	*/
	virtual Management::ROTSystemStatus_ptr status() throw (CORBA::SystemException);
	
	/**
     * Returns a reference to the scheduleName property implementation of IDL interface.
	 * @return pointer to read-only ROstring property status
	*/
	virtual ACS::ROstring_ptr scheduleName() throw (CORBA::SystemException);
	
	/**
     * Returns a reference to the scheduleLine property implementation of IDL interface.
	 * @return pointer to read-only ROlong property status
	*/
	virtual ACS::ROlong_ptr scanNumber() throw (CORBA::SystemException);
	
	/**
     * Returns a reference to the tracking property implementation of IDL interface.
	 * @return pointer to read-only ROTBoolean property tracking
	*/
	virtual Management::ROTBoolean_ptr tracking() throw (CORBA::SystemException);
	
	/**
     * Returns a reference to the currentDevice property implementation of IDL interface.
	 * @return pointer to read-only ROlong property
	*/	
	virtual ACS::ROlong_ptr currentDevice() throw (CORBA::SystemException);
	
	/**
	 * This method implements the command line interpreter. The interpreter allows to ask for services or to issue commands
	 * to the control system by human readable command lines.
	 * @throw CORBA::SystemException
	 * @throw ManagementErrors::CommandLineErrorEx
	 * @param cmd string that contains the command line
	 * @return the string that contains the answer to the command.
	 */
	virtual char * command(const char *cmd) throw (CORBA::SystemException,ManagementErrors::CommandLineErrorEx);
	
	/**
	 * This method can load a new schedule file. After a succesfull parse of the schedule the itself is started
	 * from the specified line
	 * @throw CORBA::SystemException
	 * @throw ComponentErrors::ComponentErrorsEx
	 * @throw ManagementErrors::ManagementErrorsEx
	 * @param fileName name of the file of the schedule
	 * @param startLine line number to start from
	 */
	virtual void startSchedule(const char * fileName,CORBA::Long startLine) throw (CORBA::SystemException,
			ComponentErrors::ComponentErrorsEx,ManagementErrors::ManagementErrorsEx);
	
	/**
	 * This method performs the system temperature measurment. In order to do that it uses the currently configured backend and frontend. The measure is performed
	 * for aech of the section of the backend. The sequence of operation is: 
	 * 	@arg \c a) call the backend in order to collect the informations for each section (start fequency, bandwidth, polarization and feed).
	 *     @arg \c b) call the receiverboss in order to know the value of the calibration diode for each of the sections
	 *     @arg \c c) call the backend ion order to get the TPI and TPIzero measurments
	 *     @arg \c d) call the receiverBoss in order to turn the calibration diode on
	 *     @arg \c  e) call the backend in order to get the TPIcal  measurment
	 *     @arg \c  f) call the freceiver boss in order to turn the calibration diode off
	 *     @arg \c g) compute the system temperature
	 * @throw CORBA::SystemException
	 * @throw ComponentErrors::ComponentErrorsEx
	 * @throw ManagementErrors::ManagementErrorsEx
	 * @return the system temperature in Kelvin related with each of the backend sections.
	 */
	virtual ACS::doubleSeq *systemTemperature() throw (CORBA::SystemException,ComponentErrors::ComponentErrorsEx,ManagementErrors::ManagementErrorsEx);
	
	/**
	 * This method implements the cross scan operation. The cross scan is done around a previously commanded target (prerequisite) using the main drives of the telescope and 
	 * the Antenna subsystem. The operation consists of a tsys measurment and the the lon/lat On-The-Fly scans.
 	 * @throw CORBA::SystemException
	 * @throw ComponentErrors::ComponentErrorsEx
	 * @throw ManagementErrors::ManagementErrorsEx
	 * @param coordFrame allow to select the frame along which the single scan is done
	 * @param span length of the scan
	 * @param duration how long the single subscan has to take
	 */ 
	virtual void crossScan(Management::TCoordinateFrame coordFrame,CORBA::Double span,ACS::TimeInterval duration) throw (CORBA::SystemException,
			ComponentErrors::ComponentErrorsEx,ManagementErrors::ManagementErrorsEx);
	
	/**
	 * It allows to change the name of the instance of the current default backend
	 * @throw CORBA::SystemException 
	 * @param bckInstance name of the instance 
	*/ 
	virtual void chooseDefaultBackend(const char *bckInstance) throw (CORBA::SystemException);
	
	/**
	 * It allows to change the name of the instance of the current data recorder component
	 * @throw CORBA::SystemException 
	 * @param rvcInstance name of the instance 
	 */
	virtual void chooseDefaultDataRecorder(const char *rcvInstance) throw (CORBA::SystemException);
	
	/**
	 * This method stops the current runnig schedule, if any. In order to start a new schedule the current one must be stopped.
	 * @throw CORBA::SystemException
	 **/
	virtual void stopSchedule() throw (CORBA::SystemException);
	
	/**
	 * This method halts the current runnig schedule, if any. The schedule is stopped after that the currently running scan is completed. 
	 * In order to start a new schedule the current one must be stopped.
	 * @throw CORBA::SystemException
	 **/
	virtual void haltSchedule() throw (CORBA::SystemException);
	
	/**
	 * This method will resetthe property <i>status</i> to MNG_OK.
	 * @throw CORBA::SystemException 
	 */ 
	virtual void clearStatus() throw (CORBA::SystemException);
	
	/**
	 * This method will set the current device. The current device is a section of the currently used backend. When a new device is set some configuration are done.
	 * first of all the configuration of the section is read, if everything works the new beamsize is computed.The current device is the backend section also used for calibration purposes.
	 * By default the device (if not esplicitly set) is the number zero, but if this method is not called there is no guarantee that the backend configuration, the beamsize, the receiver configuration,
	 * and the current device are coherent. the last set value can be read in the <i>currentDevice</i> attribute.
	 * @param deviceID identifier of the section (of the current backend). The value is checked to be inside valid sections identifiers.  
	 */
	virtual void setDevice(CORBA::Long deviceID) throw (CORBA::SystemException,ComponentErrors::ComponentErrorsEx,ManagementErrors::ManagementErrorsEx);
	
	
private:
	SmartPropertyPointer<ROstring> m_pscheduleName;
	SmartPropertyPointer < ROEnumImpl<ACS_ENUM_T(Management::TSystemStatus),POA_Management::ROTSystemStatus> > m_pstatus;
	SmartPropertyPointer<ROlong> m_pscanNumber;
	SmartPropertyPointer< ROEnumImpl<ACS_ENUM_T(Management::TBoolean),POA_Management::ROTBoolean> > m_ptracking;
	SmartPropertyPointer<ROlong> m_pcurrentDevice;
	CConfiguration m_config;
	CCore *m_core;
};


#endif /*SCHEDULERIMPL_H_*/
