/*******************************************************************************\
 *  Author Infos
 *  ============
 *  Name:         Marco Buttu
 *  E-mail:       mbuttu@oa-cagliari.inaf.it
 *  Personal Web: http://www.pypeople.com/
\*******************************************************************************/

#ifndef __MINORSERVOIMPL_H__
#define __MINORSERVOIMPL_H__

#ifndef __cplusplus
#error This is a C++ include file and cannot be used from plain C
#endif

#include <baciCharacteristicComponentImpl.h>
#include <baciSmartPropertyPointer.h>
#include <baciRWdouble.h>
#include <baciROdouble.h>
#include <baciRWdoubleSeq.h>
#include <baciROdoubleSeq.h>
#include <baciROlong.h>
#include <baciROpattern.h>
#include <baciROuLongLong.h>
#include <acsncSimpleSupplier.h>
#include <enumpropROImpl.h>
#include <ComponentErrors.h>
#include <MinorServoErrors.h>
#include <MinorServoS.h>
#include "RequestScheduler.h"
#include "SocketListener.h"
#include "WPStatusUpdater.h"
#include "MSParameters.h"
#include <map>

using namespace baci;

// Use the following declarations to avoid a forward declaration
struct ThreadParameters;
class WPServoTalker;
class WPServoSocket;

class WPServoImpl: public CharacteristicComponentImpl,  public virtual POA_MinorServo::WPServo {

public:
    
    WPServoImpl(const ACE_CString &CompName, maci::ContainerServices *containerServices);

    virtual ~WPServoImpl(); 

    /**
     * Get the parameter from CDB and create a WPServoSocket and a CSecureArea. 
     * Initialize the socket calling its Init method.
     *
     * @throw ComponentErrors::CDBAccessExImpl
     * @throw ComponentErrors::MemoryAllocationExImpl
     */
    virtual void initialize() throw (
            ComponentErrors::CDBAccessExImpl, 
            ComponentErrors::MemoryAllocationExImpl,
            ComponentErrors::ThreadErrorExImpl
    );


    /**
     * @throw ComponentErrors::MemoryAllocationExImpl
     */
    virtual void execute() throw (
            ComponentErrors::MemoryAllocationExImpl,
            ComponentErrors::SocketErrorExImpl
    );
    

    /** 
     * Called by the container before destroying the server in a normal situation. 
     * This function takes charge of releasing all resources.
     */
     virtual void cleanUp();
    

    /** 
     * Called by the container in case of error or emergency situation. 
     * This function tries to free all resources even though there is no warranty that the function 
     * is completely executed before the component is destroyed.
     */ 
    virtual void aboutToAbort();

    
    /** 
     * Return a reference to actPos property (actual position, ROdoubleSeq) 
     *
     * The actPos property is a ROdoubleSeq of coordinates. The reference
     * system of coordinates is a "virtual system", that is a system binds
     * with minor servo, and in general it is not bind with minor servo 
     * positioner. 
     *
     * @return pointer to ROdoubleSeq actPos property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr actPos() throw (CORBA::SystemException);


    /** 
     * Set the cmdPos property (commanded position, RWdoubleSeq) 
     *
     * The cmdPos property is a RWdoubleSeq of coordinates. The reference
     * system of coordinates is a "virtual system", that is a system binds
     * with minor servo, and in general it is not bind with minor servo 
     * positioner. This property is set by a call to setpos method.
     *
     * @return pointer to RWdoubleSeq cmdPos property
     * @throw CORBA::SystemException
     */
     virtual ACS::RWdoubleSeq_ptr cmdPos() throw (CORBA::SystemException);


    /** 
     * Return a reference to posDiff property (position difference, ROdoubleSeq) 
     *
     * The posDiff property is a ROdoubleSeq of coordinates. The reference
     * system of coordinates is a "virtual system", that is a system binds
     * with minor servo, and in general it is not bind with minor servo 
     * positioner. Every item is a difference (in the virtual reference system) 
     * between the related items of commanded and actual position 
     *
     * @return pointer to ROdoubleSeq posDiff property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr posDiff() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to engTemperature property (ROdoubleSeq) 
     *
     * The engTemperature property is a ROdoubleSeq of temperatures of every
     * minor servo slave. 
     *
     * @return pointer to ROdoubleSeq engTemperature property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr engTemperature() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to counturingErr property (ROdoubleSeq) 
     *
     * The counturingErr property is a ROdoubleSeq of counturing errors of every
     * minor servo slave. 
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr counturingErr() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to torquePerc property (ROdoubleSeq) 
     *
     * The torquePerc property is a ROdoubleSeq of torque percentages of every
     * minor servo slave. 
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr torquePerc() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to engCurrent property (ROdoubleSeq) 
     *
     * The engCurrent property is a ROdoubleSeq of engine currents of every
     * minor servo slave. 
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr engCurrent() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to engVoltage property (ROdoubleSeq) 
     *
     * The engVoltage property is a ROdoubleSeq of engine voltages of every
     * minor servo slave. 
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr engVoltage() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to driTemperature property (ROdoubleSeq) 
     *
     * The driTemperature property is a ROdoubleSeq of temperatures of every
     * slave driver. 
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr driTemperature() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to utilizationPerc property (ROdoubleSeq) 
     *
     * The utilizationPerc property is a ROdoubleSeq of utilization percentages of every
     * minor servo slave. 
     * Each item is the utilization percentage of minor servo slave.
     * That percentage is the ratio between actual power and maximum supplying power.
     *
     * @return pointer to ROdoubleSeq property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr utilizationPerc() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to dcTemperature property (ROdoubleSeq) 
     *
     * @return pointer to ROdoubleSeq dcTemperature property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr dcTemperature() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to driverStatus property (ROdoubleSeq) 
     *
     * @return pointer to ROdoubleSeq driverStatus property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr driverStatus() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to errorCode property (ROdoubleSeq) 
     *
     * @return pointer to ROdoubleSeq errorCode property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROdoubleSeq_ptr errorCode() throw (CORBA::SystemException);
    

    /** 
     * Return a reference to status property (ROpattern) 
     *  When the component is not initialized (by a setup) the pattern will be 00000X;
     *  The returned pattern is a 4 bits field with the following meaning:
     *
     *  @arg \c 0 Ready: 1 means the MinorServo is ready and can be positioned
     *  @arg \c 1 Warning: 1 means there is an abnormal condition that should not affect the observation
     *  @arg \c 2 Failure: 1 means there is a failure condition, the observation will not be succesful
     *  @arg \c 3 Tracking: 1 means the MinorServo is tracking the commanded position
     *
     * @return pointer to ROpattern status property
     * @throw CORBA::SystemException
     */
     virtual ACS::ROpattern_ptr status() throw (CORBA::SystemException);


    /**
     * Set the position of MinorServo. The time of execution command is an input parameter.
     * 
     * @arg doubleSeq position sequence of virtual axis
     * @arg long  exe_time execution time
     * @throw MinorServoErrors::PositioningErrorEx
     * @throw MinorServoErrors::CommunicationErrorEx
     */ 
     virtual void setPosition(const ACS::doubleSeq &position, const ACS::Time exe_time) 
         throw (MinorServoErrors::PositioningErrorEx, MinorServoErrors::CommunicationErrorEx);


    /**
     * Accomplish a setup of minor servo.
     * 
     * @arg long  exe_time execution time
     * @throw MinorServoErrors::SetupErrorEx
     * @throw MinorServoErrors::CommunicationErrorEx
     */ 
     virtual void setup(const ACS::Time exe_time) 
         throw (MinorServoErrors::SetupErrorEx, MinorServoErrors::CommunicationErrorEx);
     
     
     virtual void setStatusUpdating(bool flag);
     
     virtual bool isStatusThreadEn();


    /**
     * Accomplish a stow of minor servo.
     * 
     * @arg long  exe_time execution time
     * @throw MinorServoErrors::StowErrorEx
     * @throw MinorServoErrors::CommunicationErrorEx
     */ 
     virtual void stow(const ACS::Time exe_time) 
         throw (MinorServoErrors::StowErrorEx, MinorServoErrors::CommunicationErrorEx);


    /**
     * Accomplish a calibration of minor servo.
     * 
     * @arg long  exe_time execution time
     * @throw MinorServoErrors::CalibrationErrorEx
     * @throw MinorServoErrors::CommunicationErrorEx
     */ 
     virtual void calibrate(const ACS::Time exe_time) 
         throw (MinorServoErrors::CalibrationErrorEx, MinorServoErrors::CommunicationErrorEx);

     
     /** Return true when the minor servo is parked */
     virtual bool isParked(); 

     
     /** Return true when the minor servo is tracking */
     virtual bool isTracking(); 

     
     /** Return true when the minor servo is ready */
     virtual bool isReady(); 

     
     /** Return true when the minor servo is starting */
     virtual bool isStarting(); 

     /**
      * Return the data required as a sequence, each item is related to the corresponding axis
      * 
      * @param string data_name is the name of the data we would retrieve. Possible values are:
      * <ul>
      *     <li>POS_LIMIT</li>
      *     <li>NEG_LIMIT</li>
      *     <li>ACCELERATION</li>
      *     <li>MAX_SPEED</li>
      * </ul>
      * return doubleSeq one data for each axis.
      * @throw MinorServoErrors::CommunicationErrorEx
      */
     virtual ACS::doubleSeq * getData(const char *data_name) throw (MinorServoErrors::CommunicationErrorEx);


private:

    /** 
     * This attribute is a map of flags used by the StatusUpdater to enable the servo status updating.
     * The key is the servo address and the value is the corresponding flag.
     */
    static map<int, bool> m_status_thread_en;

    /** Structure containing the CDB parameters */
    CDBParameters *m_cdb_ptr;

    /* A list of vectors of commanded positions. The item are indexed by servo address */
    static CSecureArea< map<int, vector<PositionItem> > > *m_cmdPos_list;

    /** Structure containing the ExpireTime values */
    static ExpireTime m_expire;

    /** Secure area for the socket reference */
    WPServoSocket *m_wpServoLink_ptr;

    /** Actual position property */
    SmartPropertyPointer<ROdoubleSeq> m_actPos;

    /** Commanded position property */
    SmartPropertyPointer<RWdoubleSeq> m_cmdPos;

    /** Position difference property */
    SmartPropertyPointer<ROdoubleSeq> m_posDiff;

    /** Engine Temperature property */
    SmartPropertyPointer<ROdoubleSeq> m_engTemperature;

    /** Counturing Error property */
    SmartPropertyPointer<ROdoubleSeq> m_counturingErr;

    /** Torque Percentage property */
    SmartPropertyPointer<ROdoubleSeq> m_torquePerc;

    /** Engine Current property */
    SmartPropertyPointer<ROdoubleSeq> m_engCurrent;

    /** Engine Voltage property */
    SmartPropertyPointer<ROdoubleSeq> m_engVoltage;

    /** Driver Temperature property */
    SmartPropertyPointer<ROdoubleSeq> m_driTemperature;

    /** Utilization Percentage property */
    SmartPropertyPointer<ROdoubleSeq> m_utilizationPerc;

    /** Drive Cabinet property */
    SmartPropertyPointer<ROdoubleSeq> m_dcTemperature;

    /** Driver Status property */
    SmartPropertyPointer<ROdoubleSeq> m_driverStatus;

    /** Error Code property */
    SmartPropertyPointer<ROdoubleSeq> m_errorCode;

    /** Status property */
    SmartPropertyPointer<ROpattern> m_status;

    /** @var static pointer to scheduler thread */
    static RequestScheduler *m_scheduler_ptr;

    /** @var static pointer to listener thread */
    static SocketListener *m_listener_ptr;

    /** @var static pointer to status thread */
    static WPStatusUpdater *m_status_ptr;

    /** @var pointer to the WPServoTalker of minor servo */
    WPServoTalker *m_wpServoTalker_ptr;

    /** @var map of WPServoTalkers; the int key is the minor servo address. 
     *  The minor servos share this map. 
     */
    static map<int, WPServoTalker *> m_talkers;

    /** @var thread parameters */
    static ThreadParameters m_thread_params;

    /** @var Instance counter  */
    static CSecureArea<unsigned short> *m_instance_counter;

    void operator=(const WPServoImpl &);

};

#endif
