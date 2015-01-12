/*******************************************************************************\
 *  Author Infos
 *  ============
 *  Name:         Marco Bartolini
 *  E-mail:       bartolini@ira.inaf.it
\*******************************************************************************/

#include "MinorServoBossImpl.h"
#include <acsThread.h>
#include <acsThreadManager.h>

using namespace SimpleParser;

MinorServoBossImpl::MinorServoBossImpl(
        const ACE_CString &CompName, maci::ContainerServices *containerServices
        ) : 
    CharacteristicComponentImpl(CompName, containerServices),
    m_services(containerServices),
    m_actual_config(NULL),
    m_actual_conf("unconfigured"),
    m_status(this),
    m_ready(this),
    m_actualSetup(this),
    m_starting(this),
    m_asConfiguration(this),
    m_elevationTrack(this),
    m_scanActive(this),
    m_scanning(this),
    // m_verbose_status(this),
    m_nchannel(NULL),
    m_publisher_thread_ptr(NULL),
    m_tracking_thread_ptr(NULL),
    m_position_monitoring_thread_ptr(NULL)
{   
    AUTO_TRACE("MinorServoBossImpl::MinorServoBossImpl()");
    m_commanded_conf = "none";
    //m_scan_thread_ptr = NULL;
    m_servo_scanned = "none";
    m_parser = new SimpleParser::CParser<MinorServoBossImpl>(this, 10);
}


MinorServoBossImpl::~MinorServoBossImpl() { 
    AUTO_TRACE("MinorServoBossImpl::~MinorServoBossImpl()"); 
    if(m_parser != NULL)
        delete m_parser;
}


void MinorServoBossImpl::initialize() throw (ComponentErrors::CouldntGetComponentExImpl, ManagementErrors::ConfigurationErrorExImpl)
{
    AUTO_TRACE("MinorServoBossImpl::initialize()");

    /** 
     * INITIALIZE NOTIFICATION CHANNEL
     */
    try {
        m_nchannel = new nc::SimpleSupplier(MinorServo::MINORSERVO_DATA_CHANNEL, this);
        CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::initialize",
                   (LM_DEBUG, "Initialized NC supplier"));
    }
    catch (...) {
        _THROW_EXCPT(ComponentErrors::UnexpectedExImpl,"MinorServoBoss::initialize()");
    }

    /**
     * READ POSSIBLE SETUPS FROM CDB
     */
    m_config = get_configuration_from_CDB(m_services);
    CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::initialize",
              (LM_DEBUG, "Read setups from CDB"));
    /**
     * READ ATTRIBUTES FROM CDB
     */
    if(!IRA::CIRATools::getDBValue(m_services, "server_ip", m_server_ip))
        THROW_EX(ComponentErrors,CDBAccessEx, 
                 "cannot read server_ip from CDB",
                 false);
    CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::initialize",
              (LM_DEBUG, "server ip: %s", (const char*)m_server_ip));
    /**
     * INITIALIZE STATUS
     */
    m_servo_status.starting = false;
    m_servo_status.ready = false;
    m_servo_status.scan_active = false;
    m_servo_status.scanning = false;
    m_servo_status.parking = false;
    m_servo_status.parked = false;

    /**
     * INITIALIZE SERVO CONTROL
     */
    //TODO: add exception management here
    m_control = get_servo_control(m_server_ip);
    CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::initialize",
              (LM_DEBUG, "Instantiated new minor servo control"));

    /**
     * INITIALIZE PARSER WITH COMMANDS
     */
    m_parser->add(
            "servoPark", 
            new function0<MinorServoBossImpl, non_constant, void_type>(this, &MinorServoBossImpl::parkImpl), 
            0
    );
    m_parser->add(
            "servoSetup", 
            new function1<MinorServoBossImpl, non_constant, void_type, I<string_type> >(this, &MinorServoBossImpl::setupImpl), 
            1
    );
    m_parser->add(
            "setServoElevationTracking", 
            new function1<MinorServoBossImpl, non_constant, void_type, I<string_type> >(this, &MinorServoBossImpl::setElevationTrackingImpl), 
            1
    );
    m_parser->add(
            "setServoASConfiguration", 
            new function1<MinorServoBossImpl, non_constant, void_type, I<string_type> >(this, &MinorServoBossImpl::setASConfigurationImpl), 
            1
    );
    m_parser->add(
            "clearServoOffsets", 
            new function0<MinorServoBossImpl, non_constant, void_type>(this, &MinorServoBossImpl::clearOffsetsFromOI), 
            0
    );
 
    m_parser->add(
            "setServoOffset", 
             new function2<MinorServoBossImpl, non_constant, void_type, I<string_type>, I<double_type> >(this, &MinorServoBossImpl::setUserOffsetFromOI), 
             2
    );
}


void 
MinorServoBossImpl::execute() throw (ComponentErrors::MemoryAllocationExImpl)
{
    AUTO_TRACE("MinorServoBossImpl::execute()");

    try {       
        m_status = new ROEnumImpl<ACS_ENUM_T(Management::TSystemStatus), 
                                  POA_Management::ROTSystemStatus>(
                                    getContainerServices()->getName() + ":status", 
                                    getComponent(), 
                                    new DevIOStatus(&m_servo_status), 
                                    true
                                  );  
        m_ready = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":ready", getComponent(), \
                   new DevIOReady(&m_servo_status), true);
		m_actualSetup = new baci::ROstring(getContainerServices()->getName() + ":actualSetup",
                getComponent(), new DevIOActualSetup(&m_actual_config), true);
        m_starting = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":starting",getComponent(), \
                   new DevIOStarting(&m_servo_status), true);
        m_asConfiguration = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":asConfiguration",getComponent(), \
                   new DevIOASConfiguration(&m_servo_status), true);
        m_elevationTrack = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":elevationTrack",getComponent(), \
                   new DevIOElevationTrack(&m_servo_status), true);
        m_scanActive = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":scanActive",getComponent(), \
                   new DevIOScanActive(&m_servo_status), true);
        m_scanning = new ROEnumImpl<ACS_ENUM_T(Management::TBoolean), POA_Management::ROTBoolean>\
                  (getContainerServices()->getName()+":scanning",getComponent(), \
                   new DevIOScanning(&m_servo_status), true);
    }
    catch (std::bad_alloc& ex) {
        THROW_EX(ComponentErrors, MemoryAllocationEx, "MinorServoBoss::execute(): 'new' failure", false);
    }
    catch (...) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::execute: unexpected exception allocating status properties"));
    }

    if(m_publisher_thread_ptr != NULL) {
        m_publisher_thread_ptr->suspend();
        m_publisher_thread_ptr->terminate();
        m_publisher_thread_ptr = NULL;
    }
    try {
        PublisherThreadParameters thread_params(&m_servo_status, 
                                                m_nchannel,
                                                m_control);
        m_publisher_thread_ptr = getContainerServices()->
                                 getThreadManager()->
                                 create<MSBossPublisher, PublisherThreadParameters>
                                 (PUBLISHER_THREAD_NAME, thread_params);
        m_publisher_thread_ptr->resume();
        CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::execute",
                   (LM_DEBUG, "Started publisher thread"));
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "Error: MSBossPublisher already exists", true);
    }
    m_position_monitoring_thread_ptr = getContainerServices()->
                                       getThreadManager()->
                                       create<PositionMonitoringThread, MedMinorServoControl_sp>
                                       (POSITION_MONITORING_THREAD_NAME, m_control);
    m_position_monitoring_thread_ptr->resume();
    CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::execute",
              (LM_DEBUG, "Position Monitoring Thread started"));
    startPropertiesMonitoring();
}

void 
MinorServoBossImpl::cleanUp() 
{
    AUTO_TRACE("MinorServoBossImpl::cleanUp()");
    stopPropertiesMonitoring();

    try {
        turnTrackingOff(); // Raises ConfigurationError
    }
    catch(...){
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::cleanUp(): some problems turning the elevation tracking off."));
    }
    if (m_publisher_thread_ptr != NULL) {
        m_publisher_thread_ptr->suspend();
        m_publisher_thread_ptr->terminate();
    }
    if (m_position_monitoring_thread_ptr != NULL) {
        m_position_monitoring_thread_ptr->suspend();
        m_position_monitoring_thread_ptr->terminate();
    }
    if(m_nchannel !=NULL ) {
        m_nchannel->disconnect();
        m_nchannel = NULL;
    }
    disconnect();
    CharacteristicComponentImpl::cleanUp(); 
}


void 
MinorServoBossImpl::aboutToAbort()
{
    AUTO_TRACE("MinorServoBossImpl::aboutToAbort()");
    stopPropertiesMonitoring();
    try {
        turnTrackingOff();
    }
    catch(...){
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::aboutToAbort(): some problems turning the elevation tracking off."));
    }
    if (m_publisher_thread_ptr != NULL) {
        m_publisher_thread_ptr->suspend();
        m_publisher_thread_ptr->terminate();
    }
    if (m_position_monitoring_thread_ptr != NULL) {
        m_position_monitoring_thread_ptr->suspend();
        m_position_monitoring_thread_ptr->terminate();
    }
    if(m_nchannel != NULL ) {
        m_nchannel->disconnect();
        m_nchannel = NULL;
    }
    disconnect();
    CharacteristicComponentImpl::aboutToAbort(); 
}

void 
MinorServoBossImpl::setup(const char *config) 
throw (CORBA::SystemException, ManagementErrors::ConfigurationErrorEx)
{
    AUTO_TRACE("MinorServoBossImpl::setup()");
    try {
        setupImpl(config);
        CUSTOM_LOG(LM_FULL_INFO, "MinorServo::MinorServoBossImpl::setup",
                   (LM_NOTICE, "Minor Servo Setup: %s", config));
    }
    catch (ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
}

void 
MinorServoBossImpl::setupImpl(const char *config) 
throw (ManagementErrors::ConfigurationErrorExImpl)
{
    try {
        setElevationTrackingImpl(IRA::CString("OFF"));
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "cannot turn the tracking off", false);
    }
    //TODO: define a new Setup Error
    if(m_servo_status.starting)
        THROW_EX(ManagementErrors, ConfigurationErrorEx, 
                 "The system is executing another setup", false);
    if(m_servo_status.parking)
        THROW_EX(ManagementErrors, ConfigurationErrorEx, 
                 "The system is executing a park", false);
    if(m_servo_status.scanning)
        THROW_EX(ManagementErrors, ConfigurationErrorEx, 
                 "The system is performing a scan", false);
    if(m_config.count(std::string(config)) == 0)
    {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, 
                 "Cannot find requested configuration", false);
    }else{
        m_servo_status.reset();
        m_servo_status.starting = true;
        m_actual_config = &(m_config[std::string(config)]);
        /**
         * Get the setup position at 45 deg
         */
        MedMinorServoPosition setup_position = m_actual_config->get_position();
        m_offset.initialize(m_actual_config->is_primary_focus());
        //clearUserOffset();
        //clearSystemOffset();
        m_control->set_position(setup_position);
        m_servo_status.ready = true;
        m_servo_status.starting = false;
        m_actual_conf = string(config);
    }
}

void 
MinorServoBossImpl::park() 
throw (CORBA::SystemException, ManagementErrors::ParkingErrorEx)
{
    try {
        parkImpl();
    }
    catch (ManagementErrors::ParkingErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getParkingErrorEx();       
    }
}

void 
MinorServoBossImpl::parkImpl() 
throw (ManagementErrors::ParkingErrorExImpl)
{
    AUTO_TRACE("MinorServoBossImpl::parkImpl()");

    if(m_servo_status.starting) {
        THROW_EX(ManagementErrors, ParkingErrorEx, "The system is executing a setup.", false);
    }

    if(m_servo_status.parking)
        THROW_EX(ManagementErrors, ParkingErrorEx, "The system is executing another park.", false);

    try { 
        setElevationTrackingImpl(IRA::CString("OFF")); // Raises ConfigurationError
    }
    catch(...) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::park(): some problems turning the elevation tracking off."));
    }

    m_servo_status.parking = true;

    /*try {
        if(m_park_thread_ptr != NULL)
            m_park_thread_ptr->restart();
        else {
            m_park_thread_ptr = getContainerServices()->getThreadManager()->create<ParkThread,
               MSBossConfiguration *>("ParkThread", m_configuration);
            m_park_thread_ptr->resume();
        }
    }
    catch(...) {
        THROW_EX(ManagementErrors, ParkingErrorEx, "The MinorServoBoss is attempting to execute a previous park", false);
    }*/

    m_servo_status.parking = false;
}


void 
MinorServoBossImpl::getAxesInfo(ACS::stringSeq_out axes, ACS::stringSeq_out units)
throw (CORBA::SystemException, ManagementErrors::ConfigurationErrorEx) 
{

    //TODO: substitue ConfigurationError with something state-related
    if(!isReady())
        THROW_EX(ManagementErrors, 
                 ConfigurationErrorEx, 
                 "getAxesInfo(): the system is not ready", 
                 true);

    ACS::stringSeq_var axes_res = new ACS::stringSeq;
    ACS::stringSeq_var units_res = new ACS::stringSeq;

    std::vector<std::string> a = m_actual_config->getAxes();
    std::vector<std::string> u = m_actual_config->getUnits();
    axes_res->length(a.size());
    units_res->length(u.size());

    //TODO: rename this exception
    if(a.size() != u.size())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getAxesInfo(): mismatch between axes and units length", true);

    for(size_t i=0; i<a.size(); i++) {
        axes_res[i] = (a[i]).c_str();
        units_res[i] = (u[i]).c_str();
    }

    axes = axes_res._retn();
    units = units_res._retn();
}


ACS::doubleSeq * 
MinorServoBossImpl::getAxesPosition(ACS::Time time) 
    throw (CORBA::SystemException, ManagementErrors::ConfigurationErrorEx, ComponentErrors::UnexpectedEx)
{
    //TODO: substitue ConfigurationError with something state-related
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getAxesPosition(): the system is not ready", true);
    ACS::doubleSeq_var positions_res = new ACS::doubleSeq;
    vector<double> vpositions;
    MedMinorServoPosition position;
    if(time == 0)
        position = m_control->get_position();
    else
        position = m_control->get_position_at(time);
    vpositions = position.get_axes_positions();
    positions_res->length(vpositions.size());
    for(size_t i=0; i<vpositions.size(); i++)
        positions_res[i] = vpositions[i];
    return positions_res._retn();
}


bool MinorServoBossImpl::isStarting() { return m_servo_status.starting; }
bool MinorServoBossImpl::isASConfiguration() { return m_servo_status.as_configuration; }
bool MinorServoBossImpl::isParking() { return m_servo_status.parking; }
bool MinorServoBossImpl::isReady() { return m_servo_status.ready; }
bool MinorServoBossImpl::isScanning() { return m_servo_status.scanning; }
bool MinorServoBossImpl::isScanActive() { return m_servo_status.scan_active; }


CORBA::Boolean 
MinorServoBossImpl::command(const char *cmd, CORBA::String_out answer) 
throw (CORBA::SystemException)
{
    AUTO_TRACE("MinorServoBossImpl::command()");

    IRA::CString out;
    bool res;
    try {
        m_parser->run(cmd, out);
        res = true;
    }
    catch(ParserErrors::ParserErrorsExImpl& ex) {
        res = false;
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_ERROR); 
        res = false;
    }
    catch(ACSErr::ACSbaseExImpl& ex) {
        ex.log(LM_ERROR); 
        res = false;
    }
    catch(...) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::command(): unknown exception."));
        res = false;
    }
    answer = CORBA::string_dup((const char *)out);
    return res;
}

bool 
MinorServoBossImpl::isParked() 
throw (ManagementErrors::ConfigurationErrorEx) 
{
    if(isStarting() || isParking() || isReady())
        return false;
    return true;
}

void MinorServoBossImpl::stopScan() throw (ManagementErrors::SubscanErrorEx)
{
/*    if(m_configuration->m_isScanActive) {
        string comp_name((m_configuration->m_scan).comp_name);
        if(m_component_refs.count(comp_name)) {
            MinorServo::WPServo_var component_ref = m_component_refs[comp_name];
            try {
                MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();
                component_ref = m_services->getComponent<MinorServo::WPServo>(("MINORSERVO/" + comp_name).c_str());
                if(!CORBA::is_nil(component_ref)) {
                    if(m_scan_thread_ptr != NULL) {
                        m_scan_thread_ptr->terminate();
                        m_scan_thread_ptr->suspend();
                    }
                    component_ref->cleanPositionsQueue(NOW);
                    m_configuration->m_isScanning = false;
                    if((m_configuration->m_scan).wasElevationTrackingEn)
                        turnTrackingOn();
                    else {
                        if(component_ref->isReady()) {
                            component_ref->setPosition((m_configuration->m_scan).plainCentralPos, NOW);
                        }
                        else {
                            THROW_EX(ManagementErrors, SubscanErrorEx, "stopScan(): component not ready", true);
                        }
                    }

                    m_configuration->m_isScanActive = false;
                }
                else {
                    THROW_EX(ManagementErrors, SubscanErrorEx, "stopScan(): nil component reference", true);
                }
            }
            catch(...) {
                THROW_EX(ManagementErrors, SubscanErrorEx, "stopScan(): unexpected exception", true);
            }
        }
        else {
            THROW_EX(ManagementErrors, SubscanErrorEx, "stopScan(): cannot get the component reference", true);
        }
    }
    else {
        THROW_EX(ManagementErrors, SubscanErrorEx, "stopScan(): no scan active", true);
    }*/
}


bool MinorServoBossImpl::checkFocusScan(const ACS::Time starting_time, const double range, const ACS::Time total_time) 
        throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx) 
{    
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "checkFocusScan: the system is not ready", true);
    if(isScanning())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "checkFocusScan: the system is executing another scan", true);
    
    /*string servo_name =  m_configuration->getActivePFocusServo(); 
    string axis_code = servo_name + string("_TZ");
    return checkScanImpl(starting_time, range, total_time, axis_code);*/
    return true;
}


bool 
MinorServoBossImpl::checkScan(
        const ACS::Time starting_time, 
        const double range, 
        const ACS::Time total_time, 
        const char *axis_code
        )
    throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx)
{
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "StartScan: the system is not ready", true);

    if(isScanning())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "StartScan: the system is executing another scan", true);

    return checkScanImpl(starting_time, range, total_time, string(axis_code));
}


ACS::Time 
MinorServoBossImpl::getMinScanStartingTime(
        double range, 
        const string axis_code, 
        double & acceleration, 
        double & max_speed
    )
    throw (ManagementErrors::ConfigurationErrorExImpl, ManagementErrors::SubscanErrorExImpl)
{
    ACS::Time min_starting_time = 0;
    /*
    InfoAxisCode info;
    info = m_configuration->getInfoFromAxisCode(axis_code); // Throws ConfigurationErrorExImpl
    size_t axis = info.axis_id;
    string comp_name = info.comp_name;

    TIMEVALUE dtime(SCAN_DELTA_TIME);

    MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();

    if((m_configuration->m_component_refs).count(comp_name)) {
        component_ref = (m_configuration->m_component_refs)[comp_name];
        if(!CORBA::is_nil(component_ref))
            if(!component_ref->isReady()) {
                THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): the component is not ready", false);
            }
            else {
                CDB::DAL_ptr dal_p = getContainerServices()->getCDB();
                CDB::DAO_ptr dao_p = dal_p->get_DAO_Servant(("alma/MINORSERVO/" + comp_name).c_str());
                size_t number_of_axis = dao_p->get_long("number_of_axis");
                bool virtual_rs = dao_p->get_long("virtual_rs");
                long servo_address = dao_p->get_long("servo_address");
                double zero = dao_p->get_double("zero");
                if(axis > number_of_axis - 1) {
                    THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): axis out of range", false);
                }

                ACS::doubleSeq *acc = component_ref->getData("ACCELERATION");
                ACS::doubleSeq *mspeed = component_ref->getData("MAX_SPEED");
                size_t idx = (*acc).length() - number_of_axis + axis;
                if((*acc).length() - 1 < idx || (*mspeed).length() - 1 < idx) {
                    THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): index out of range", false);
                }
                acceleration = (*acc)[idx]; 
                max_speed= (*mspeed)[idx]; 
                if(acceleration != 0 && max_speed != 0) {
                    ACSErr::Completion_var completion;          
                    ACS::doubleSeq act_pos = *((component_ref->actPos())->get_sync(completion.out()));

                    if(act_pos.length() <= axis) {
                        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): wrong position indexing", false);
                    }
                    // Get the central position                  
                    ACS::doubleSeq centralPos = m_configuration->isScanActive() ? (m_configuration->m_scan).centralPos : act_pos;

                    ACS::doubleSeq positions_left, positions_right;
                    positions_left = positions_right = centralPos;
                    positions_left[axis] = positions_left[axis] - range / 2;
                    positions_right[axis] = positions_right[axis] + range / 2;
                    // Check if the positions are allowed
                    ACS::doubleSeq_var max_values = component_ref->getMaxPositions();
                    ACS::doubleSeq_var min_values = component_ref->getMinPositions();
                    ACS::doubleSeq_var system_offset = component_ref->getSystemOffset();
                    if(system_offset->length() < axis) {
                        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): wrong system offset indexing", false);
                    }
                    if(positions_left[axis] + system_offset[axis] <= min_values[axis]) {
                        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): min position out of range", false);
                    }
                    if(positions_right[axis] + system_offset[axis] >= max_values[axis]) {
                        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getMinScanStartingTime(): max position out of range", false);
                    }
                    
                    double left_distance = abs_(act_pos[axis] - positions_left[axis]);
                    double right_distance = abs_(act_pos[axis] - positions_right[axis]);
                    double positioning_distance = left_distance <= right_distance ? left_distance : right_distance;

                    // If the component has a virtual reference system
                    if(virtual_rs) {
                        // Conversion from Virtual reference system to the Real one
                        virtual2real(positions_left, servo_address, zero);
                        virtual2real(positions_right, servo_address, zero);
                        ACS::doubleSeq diff;
                        diff.length(act_pos.length());

                        for(size_t i = 0; i < act_pos.length(); i++)
                            diff[i] = abs_(positions_left[i] - positions_right[i]);

                        double max_diff = 0;
                        for(size_t i = 0; i < act_pos.length(); i++)
                            if(diff[i] > max_diff)
                                max_diff = diff[i];

                        range = max_diff;

                        // Conversion from Virtual reference system to the Real one
                        virtual2real(act_pos, servo_address, zero);
                        ACS::doubleSeq diff_left, diff_right;
                        diff_left.length(act_pos.length());
                        diff_right.length(act_pos.length());

                        for(size_t i = 0; i < act_pos.length(); i++) {
                            diff_left[i] = abs_(act_pos[i] - positions_left[i]);
                            diff_right[i] = abs_(act_pos[i] - positions_right[i]);
                        }

                        double max_diff_left = 0;
                        double max_diff_right = 0;
                        for(size_t i = 0; i < act_pos.length(); i++) {
                            if(diff_left[i] > max_diff_left)
                                max_diff_left = diff_left[i];
                            if(diff_right[i] > max_diff_right)
                                max_diff_right = diff_right[i];
                        }
                        positioning_distance = max_diff_left <= max_diff_right ? max_diff_left : max_diff_right;
                    }

                    ACS::Time positioning_time = get_min_time(positioning_distance, acceleration, max_speed);
                    
                    
                    min_starting_time = static_cast<ACS::Time>(
                            getTimeStamp() + SCAN_SHIFT_TIME + positioning_time * (1 + SCAN_GUARD_COEFF)
                    );
                }
                else {
                    THROW_EX(
                            ManagementErrors, 
                            ConfigurationErrorEx, 
                            "getMinScanStartingTime(): Acceleration or maximum speed wrong", 
                            false
                    );
                }
            }
    }
    else {
        THROW_EX(
                ManagementErrors, 
                ConfigurationErrorEx, 
                ("some problems getting the component " + comp_name).c_str(),
                false
        );
    }
    */
    return min_starting_time;
}


bool MinorServoBossImpl::checkScanImpl(
        const ACS::Time starting_time, 
        double range, 
        const ACS::Time total_time, 
        const string axis_code
    ) throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx)
{
    /*                    
    double acceleration = 0.0;
    double max_speed = 0.0;
    TIMEVALUE stime(starting_time);
    TIMEVALUE ttime(total_time);
    // (ManagementErrors::ConfigurationErrorExImpl, ManagementErrors::SubscanErrorExImpl)
    ACS::Time min_starting_time = 0;
    try {
        min_starting_time = getMinScanStartingTime(range, axis_code, acceleration, max_speed);
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        return false;
    }
    catch(ManagementErrors::SubscanErrorExImpl& ex) {
        return false;
    }
    catch(...) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::checkScan(): unexpected exception getting the min starting time"));
        return false;
    }
    
    TIMEVALUE min_stime(min_starting_time); // Minimum starting time
    if(min_starting_time > starting_time) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::checkScan(): starting time too much close to actual time"));
        return false;
    }

    ACS::Time min_scan_time = get_min_time(range, acceleration, max_speed);
    ACS::Time guard_min_scan_time = static_cast<ACS::Time>(min_scan_time * (1 + SCAN_GUARD_COEFF));
    TIMEVALUE gmst(guard_min_scan_time);
    if(CIRATools::timeSubtract(ttime, gmst) <= 0) {
        ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::checkScan(): total time too short for performing the scan."));
        return false;
    }
    else
        return true;

    ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::checkScan - last line!"));
    return false;
    */
    return true;
}


void MinorServoBossImpl::startScan(
        ACS::Time & starting_time, 
        const double range, 
        const ACS::Time total_time, 
        const char *axis_code
    ) throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx)
{
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "StartScan: the system is not ready", true);

    if(isScanning())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "StartScan: the system is executing another scan", true);

    startScanImpl(starting_time, range, total_time, string(axis_code));
}


void MinorServoBossImpl::startFocusScan(ACS::Time & starting_time, const double range, const ACS::Time total_time) 
        throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx) 
{    
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "startFocusScan: the system is not ready", true);

    if(isScanning())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "startFocusScan: the system is executing another scan", true);
    
    /**
     * we call both primary and secondary focus axes "Z"
     */
    startScanImpl(starting_time, range, total_time, "Z");
}


void MinorServoBossImpl::startScanImpl(
        ACS::Time & starting_time, 
        const double range, 
        const ACS::Time total_time, 
        string axis_code
    ) throw (ManagementErrors::ConfigurationErrorEx, ManagementErrors::SubscanErrorEx)
{
    /*
    m_configuration->m_isScanning = true;
    InfoAxisCode info;
    try {
        info = m_configuration->getInfoFromAxisCode(axis_code);
    }
    catch (ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "Unexpected exception getting the axis information.", true);
    }

    try {
        if(starting_time == 0) {
            double acceleration = 0.0;
            double max_speed = 0.0;
            starting_time = getMinScanStartingTime(range, axis_code, acceleration, max_speed);
        }
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(ManagementErrors::SubscanErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getSubscanErrorEx();     
    }
    catch(...) {
        THROW_EX(ManagementErrors, SubscanErrorEx, "Unexpected exception getting the min starting time", true);
    }

    size_t axis = info.axis_id;
    string comp_name = info.comp_name;

    try {
        TIMEVALUE now(0.0L);
        IRA::CIRATools::getTime(now);

        // If the total time is less than the delta time (the time between two consecutive positioning)
        TIMEVALUE ttime(total_time);
        TIMEVALUE dtime(SCAN_DELTA_TIME);
        if(CIRATools::timeSubtract(ttime, dtime) < 0)
            THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: total time is too short", true);

        TIMEVALUE stime(starting_time);
        if(CIRATools::timeSubtract(stime, now) <= 0)
            THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: starting time is not valid", true);

        m_servo_scanned = comp_name;
        MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();
        if((m_configuration->m_component_refs).count(comp_name)) {
            component_ref = (m_configuration->m_component_refs)[comp_name];

            if(CORBA::is_nil(component_ref)) 
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScanImpl: cannot get the reference of the component.", true);

            if(!component_ref->isReady()) 
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScanImpl: the component is not ready.", true);

            CDB::DAL_ptr dal_p = getContainerServices()->getCDB();
            CDB::DAO_ptr dao_p = dal_p->get_DAO_Servant(("alma/MINORSERVO/" + comp_name).c_str());
            size_t number_of_axis = dao_p->get_long("number_of_axis");
                                                        
            if(axis > number_of_axis - 1)               
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: axis index error", true);
                                                        
            bool wasTrackingEn = isElevationTrackingEn();
            turnTrackingOff();

            // Get the actual positions
            ACSErr::Completion_var completion;          
            ACS::doubleSeq actPos = *((component_ref->actPos())->get_sync(completion.out()));
            ACS::doubleSeq plainActPos = *((component_ref->plainActPos())->get_sync(completion.out()));

            // Compute the central positions
            ACS::doubleSeq plainCentralPos = m_configuration->isScanActive() ? (m_configuration->m_scan).plainCentralPos : plainActPos;
            ACS::doubleSeq centralPos;
            centralPos.length(component_ref->numberOfAxes());
            if(m_configuration->isScanActive()) {
                ACS::doubleSeq_var user_offset = component_ref->getUserOffset();
                if(user_offset->length() != centralPos.length()) {
                    THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: mismatch between offset and central position length.", true);
                }
                if(plainCentralPos.length() != centralPos.length()) {
                    THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: mismatch between central positions length.", true);
                }
                for(size_t i=0; i<centralPos.length(); i++) 
                    centralPos[i] = plainCentralPos[i] + user_offset[i];
            }
            else {
                if(centralPos.length() != actPos.length()) {
                    THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: mismatch between actual and central position length.", true);
                } 
                for(size_t i=0; i<actPos.length(); i++)
                    centralPos[i] = actPos[i];
            }

            // Get the virual elongations 
            ACS::doubleSeq virtualCentralElongation = m_configuration->isScanActive() ? 
                                                      (m_configuration->m_scan).virtualCentralElongation : 
                                                      plainActPos;
                                                    
            if(actPos.length() <= axis) {          
                ACS_SHORT_LOG((LM_WARNING, "MinorServoBoss::startScan: wrong actual position length"));
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: wrong actual position length.", true);
            }                                       

            ACS::doubleSeq_var max_values = component_ref->getMaxPositions();
            ACS::doubleSeq_var min_values = component_ref->getMinPositions();

            // Check if the positions are allowed
            if(virtualCentralElongation[axis] - range/2 <= min_values[axis]) {
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: min axis position out of range.", true);
            }
            if(virtualCentralElongation[axis] + range/2 >= max_values[axis]) {
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: max axis position out of range.", true);
            }

            m_configuration->setScan(
                    starting_time, 
                    total_time, 
                    SCAN_DELTA_TIME, 
                    range, 
                    comp_name, 
                    axis, 
                    axis_code,
                    actPos, 
                    centralPos,
                    plainCentralPos,
                    virtualCentralElongation,
                    wasTrackingEn
            );
        }
        else {
            THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: cannot get the component reference.", true);
        }
        
        try {
            if(m_scan_thread_ptr != NULL)
                m_scan_thread_ptr->restart();
            else {
                m_scan_thread_ptr = getContainerServices()->getThreadManager()->create<ScanThread, MSBossConfiguration *> 
                    ("ScanThread", m_configuration);
                m_scan_thread_ptr->resume();
            }
        }
        catch(...) {
            THROW_EX(ManagementErrors, ConfigurationErrorEx, "The MinorServoBoss is attempting to execute a previous scan", false);
        }
        m_configuration->m_isScanActive = true;
    }
    catch(...) {
        m_configuration->m_isScanning = false;
        m_configuration->m_isScanActive = false;
        throw;
    }
    */
}


CORBA::Double MinorServoBossImpl::getCentralScanPosition() throw (ManagementErrors::SubscanErrorEx)
{
    /*
    ACS::doubleSeq_var position_res = new ACS::doubleSeq;
    string comp_name = (m_configuration->m_scan).comp_name;

    if(isScanActive()) {
        MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();
        if((m_configuration->m_component_refs).count(comp_name)) {
            component_ref = (m_configuration->m_component_refs)[comp_name];
            if(CORBA::is_nil(component_ref)) 
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScanImpl: cannot get the reference of the component.", true);

            ACS::doubleSeq plainCentralPos = (m_configuration->m_scan).plainCentralPos;

            ACS::doubleSeq_var user_offset = component_ref->getUserOffset();
            position_res->length(user_offset->length());
            if(user_offset->length() != plainCentralPos.length()) {
                THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: mismatch between offset and central position length.", true);
            }
            for(size_t i=0; i<plainCentralPos.length(); i++) 
                position_res[i] = plainCentralPos[i] + user_offset[i];
        }
        else {
            THROW_EX(ManagementErrors, SubscanErrorEx, "startScan: cannot get the name of" + comp_name, true);
        }
    }
    else {
        THROW_EX(ManagementErrors, SubscanErrorEx, "getCentralScanPosition(): scan not active", true);
    }

    return position_res[(m_configuration->m_scan).axis_index];
    */
    return 0.0;
}


char * MinorServoBossImpl::getScanAxis() {
    if(isScanActive()) {
        //return CORBA::string_dup(((m_configuration->m_scan).axis_code).c_str());
        return CORBA::string_dup("Z");
    }
    else {
        return CORBA::string_dup("");
    }
}

void MinorServoBossImpl::turnTrackingOn() throw (ManagementErrors::ConfigurationErrorEx) 
{
    if(isStarting())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "turnTrackingOn: the system is starting.", true);
    if(isParking())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "turnTrackingOn: the system is parking.", true);
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "turnTrackingOn: the system is not ready.", true);

    if(m_tracking_thread_ptr != NULL) {
        m_tracking_thread_ptr->suspend();
        m_tracking_thread_ptr->terminate();
        m_tracking_thread_ptr = NULL;
    }

    m_servo_status.elevation_tracking = true;
    try {
        TrackerThreadParameters params(&m_servo_status,
                                        m_control,
                                        &m_actual_config,
                                        &m_offset,
                                        m_services);

        m_tracking_thread_ptr = getContainerServices()->
                                 getThreadManager()->
                                 create<MSBossTracker, TrackerThreadParameters>
                                 (TRACKING_THREAD_NAME, params);
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "Error in TrackingThread", true);
    }
    m_tracking_thread_ptr->resume();
}

void MinorServoBossImpl::turnTrackingOff() throw (ManagementErrors::ConfigurationErrorEx) 
{
    if(m_tracking_thread_ptr != NULL) {
        m_tracking_thread_ptr->suspend();
        m_tracking_thread_ptr->terminate();
        m_tracking_thread_ptr = NULL;
    }
    m_servo_status.elevation_tracking = false;
}


void MinorServoBossImpl::clearUserOffset(const char *servo) throw (MinorServoErrors::OperationNotPermittedEx){
    //TODO: add control and exception
    m_offset.clearUserOffset();
}


// Clear offset from Operator Input
void MinorServoBossImpl::clearOffsetsFromOI() throw (MinorServoErrors::OperationNotPermittedExImpl){
    try {
        clearUserOffset();
    }
    catch(...) {
        THROW_EX(
                MinorServoErrors, 
                OperationNotPermittedEx, 
                string("Cannot clear the offsets"), 
                false 
        )
    }
}


void MinorServoBossImpl::clearSystemOffset(const char *servo) throw (MinorServoErrors::OperationNotPermittedEx){
    //TODO: add control and exception
    m_offset.clearSystemOffset();
}

//void MinorServoBossImpl::clearOffset(const char *servo, string offset_type) throw (MinorServoErrors::OperationNotPermittedEx)
//{
    /* 
    string comp_name = get_component_name(string(servo));
    MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();
    if(comp_name == "ALL") {
        vector<string> slaves = m_configuration->getServosToMove();
        for(vector<string>::iterator iter = slaves.begin(); iter != slaves.end(); iter++) {
            if(m_component_refs.count(*iter)) {
                component_ref = MinorServo::WPServo::_nil();
                component_ref = m_component_refs[*iter];
                if(!CORBA::is_nil(component_ref))
                    if(offset_type == "user") {
                        if(component_ref->isReady()) {
                            try {
                                component_ref->clearUserOffset(true);
                            }
                            catch(...) {
                                THROW_EX(
                                        MinorServoErrors, 
                                        OperationNotPermittedEx, 
                                        string("Cannot clear the WPServo user offset."), 
                                        true
                                );
                            }
                        }
                        else if(!component_ref->isParked()) {
                            string msg("MinorServoBossImpl::clearOffset(): cannot clear the user offset, " + *iter + " not ready.");
                            ACS_SHORT_LOG((LM_WARNING, msg.c_str()));
                        }
                    }
                    else
                        if(offset_type == "system") {
                            if(component_ref->isReady()) {
                                try {
                                    component_ref->clearSystemOffset(true);
                                }
                                catch(...) {
                                    THROW_EX(
                                            MinorServoErrors, 
                                            OperationNotPermittedEx, 
                                            string("Cannot clear the WPServo system offset."), 
                                            true
                                    );
                                }
                            }
                            else if(!component_ref->isParked()) {
                                string msg("MinorServoBossImpl::clearOffset(): cannot clear the system offset, " + *iter + " not ready.");
                                ACS_SHORT_LOG((LM_WARNING, msg.c_str()));
                            }
                        }
                        else {
                            THROW_EX(
                                    MinorServoErrors, 
                                    OperationNotPermittedEx, 
                                    string("The offset ") + offset_type + string(" doesn't exist"), 
                                    true
                            );
                        }
            }
        }
    }
    else {
        if(!slave_exists(comp_name)) {
            THROW_EX(
                    MinorServoErrors, 
                    OperationNotPermittedEx, 
                    string("The component ") + comp_name + string(" doesn't exist"), 
                    true
            );
        }
        if(m_component_refs.count(comp_name)) {
            component_ref = m_component_refs[comp_name];
            if(!CORBA::is_nil(component_ref)) {
                    if(offset_type == "user") {
                        if(component_ref->isReady()) {
                            try {
                                component_ref->clearUserOffset(true);
                            }
                            catch(...) {
                                THROW_EX(
                                        MinorServoErrors, 
                                        OperationNotPermittedEx, 
                                        string("Cannot clear the WPServo user offset."), 
                                        true
                                );
                            }
                        }
                        else if(!component_ref->isParked()) {
                            ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::clearOffset(): cannot clear the user offset."));
                            ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::clearOffset(): SRP not ready"));
                        }
                    }
                    else
                        if(offset_type == "system") {
                            if(component_ref->isReady()) {
                                try {
                                    component_ref->clearSystemOffset(true);
                                }
                                catch(...) {
                                    THROW_EX(
                                            MinorServoErrors, 
                                            OperationNotPermittedEx, 
                                            string("Cannot clear the WPServo system offset."), 
                                            true
                                    );
                                }
                            }
                            else if(!component_ref->isParked()) {
                                ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::clearOffset(): cannot clear the system offset."));
                                ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::clearOffset(): SRP not ready"));
                            }
                        }
                        else {
                            THROW_EX(
                                    MinorServoErrors, 
                                    OperationNotPermittedEx, 
                                    string("The offset ") + offset_type + string(" doesn't exist"), 
                                    true
                            );
                        }
            }
            else {
                THROW_EX(
                        MinorServoErrors, 
                        OperationNotPermittedEx, 
                        string("The reference to component ") + comp_name + string(" is NULL"), 
                        true
                );
            }
        }
        else {
            THROW_EX(
                    MinorServoErrors, 
                    OperationNotPermittedEx, 
                    string("The component ") + comp_name + string(" is not active"), 
                    true
            );
        }
    }
    */
//}


void MinorServoBossImpl::setUserOffset(const char *axis_code, const double offset) 
    throw (MinorServoErrors::OperationNotPermittedEx, ManagementErrors::ConfigurationErrorEx) 
{
    try {
        setOffsetImpl(string(axis_code), offset, string("user"));
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(MinorServoErrors::OperationNotPermittedExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getOperationNotPermittedEx();     
    }
}

void MinorServoBossImpl::setUserOffsetFromOI(const char * axis_code, const double & offset)
         throw (MinorServoErrors::OperationNotPermittedExImpl)
{
    try {
        setUserOffset(axis_code, offset);
    }
    catch(...) {
        THROW_EX(
                MinorServoErrors, 
                OperationNotPermittedEx, 
                string("Cannot set the offset"), 
                false 
        )
    }
}

void MinorServoBossImpl::setSystemOffset(const char *axis_code, const double offset) 
    throw (MinorServoErrors::OperationNotPermittedEx, ManagementErrors::ConfigurationErrorEx) 
{
    try {
        setOffsetImpl(string(axis_code), offset, string("system"));
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(MinorServoErrors::OperationNotPermittedExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getOperationNotPermittedEx();     
    }
}


void MinorServoBossImpl::setOffsetImpl(string axis_code, const double offset_value, string offset_type) 
    throw (MinorServoErrors::OperationNotPermittedExImpl, ManagementErrors::ConfigurationErrorExImpl)
{
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "setOffsetImpl(): the system is not ready", false);

    int axis_mapping = m_actual_config->getAxisMapping(axis_code);
    if(axis_mapping < 0)
        THROW_EX(
                MinorServoErrors, 
                OperationNotPermittedEx, 
                string("Wrogn offset axis"),
                false
        );
    if(offset_type == "user")
        m_offset.setUserOffset(axis_mapping, offset_value);
    else
        m_offset.setSystemOffset(axis_mapping, offset_value);
}


ACS::doubleSeq * MinorServoBossImpl::getUserOffset() 
     throw (MinorServoErrors::OperationNotPermittedEx, ManagementErrors::ConfigurationErrorEx)
{
    vector<double> offset;
    try {
        offset = getOffsetImpl("user");
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(MinorServoErrors::OperationNotPermittedExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getOperationNotPermittedEx();     
    }

    ACS::doubleSeq_var offset_res = new ACS::doubleSeq;
    offset_res->length(offset.size());
    for(size_t i=0; i<offset_res->length(); i++)
        offset_res[i] = offset[i];
    return offset_res._retn();
}


ACS::doubleSeq * MinorServoBossImpl::getSystemOffset() 
     throw (MinorServoErrors::OperationNotPermittedEx, ManagementErrors::ConfigurationErrorEx)
{
    vector<double> offset;
    try {
        offset = getOffsetImpl("system");
    }
    catch(ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
    catch(MinorServoErrors::OperationNotPermittedExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getOperationNotPermittedEx();     
    }

    ACS::doubleSeq_var offset_res = new ACS::doubleSeq;
    offset_res->length(offset.size());
    for(size_t i=0; i<offset_res->length(); i++)
        offset_res[i] = offset[i];
    return offset_res._retn();
}


vector<double> 
MinorServoBossImpl::getOffsetImpl(string offset_type)
     throw (MinorServoErrors::OperationNotPermittedExImpl, 
            ManagementErrors::ConfigurationErrorExImpl)
{
    if(!isReady())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "getOffsetImpl(): the system is not ready", false);

    vector<double> axes_values;
    if(offset_type == "user")
        return m_offset.getUserOffset();
    else
        return m_offset.getSystemOffset();
}

void 
MinorServoBossImpl::setElevationTracking(const char * value) 
throw (ManagementErrors::ConfigurationErrorEx) {
    try {
        setElevationTrackingImpl(value);
    }
    catch (ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
}


void 
MinorServoBossImpl::setElevationTrackingImpl(const char * value) 
throw (ManagementErrors::ConfigurationErrorExImpl) {
    string flag(value);
    try{
    if(flag == "ON")
        turnTrackingOn();
    else
        turnTrackingOff();
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, string("setElevationTracking(): cannot turn the tracking") + string(flag), false);
    }
}


void MinorServoBossImpl::setASConfiguration(const char * value) throw (ManagementErrors::ConfigurationErrorEx) {
    try {
        setASConfigurationImpl(value);
    }
    catch (ManagementErrors::ConfigurationErrorExImpl& ex) {
        ex.log(LM_DEBUG);
        throw ex.getConfigurationErrorEx();     
    }
}


void MinorServoBossImpl::setASConfigurationImpl(const char * value) throw (ManagementErrors::ConfigurationErrorExImpl) {
    /*

    if(m_configuration->isStarting())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "Cannot set the ASConfiguration because the system is starting", false);

    if(m_configuration->isParking())
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "Cannot set the ASConfiguration because the system is parking", false);

    IRA::CString flag(value);
    flag.MakeUpper();
    m_configuration->setASConfiguration(flag); 

    bool wasTrackingEn = isElevationTrackingEn();
    try {
        if(isReady()) {
            bool turnedOff = false;
    
            try {
                clearUserOffset("ALL");
                clearSystemOffset("ALL");
            }
            catch(...) {
                THROW_EX(
                        ManagementErrors, 
                        ConfigurationErrorEx, 
                        string("Cannot clear the offsets"), 
                        false 
                );
            }

            try {
                turnTrackingOff(); // It raises ConfigurationErrorEx
                turnedOff = true;
            }
            catch(...) {
                THROW_EX(ManagementErrors, ConfigurationErrorEx, "cannot turn the tracking off", false);
            }
            m_configuration->init(m_configuration->m_baseSetup, true); // Keep the actual setup
            if(wasTrackingEn) {
                try {
                    turnTrackingOn(); // It raises ConfigurationErrorEx
                    turnedOff = false;
                }
                catch(...) {
                    THROW_EX(ManagementErrors, ConfigurationErrorEx, "cannot turn the tracking on", false);
                }
            }
            if(turnedOff) {
                string comp_name("SRP");
                MinorServo::WPServo_var component_ref = MinorServo::WPServo::_nil();
                if(!m_component_refs.count(comp_name)) {
                    try {
                        component_ref = m_services->getComponent<MinorServo::WPServo>(("MINORSERVO/" + comp_name).c_str());
                    }
                    catch(...) {
                        ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::setASConfiguration(): cannot get the SRP component"));
                    }
                }
                else {
                    component_ref = m_component_refs[comp_name];
                }

                ACS::doubleSeq positions = m_configuration->getPosition("SRP", getTimeStamp());
                // set the position
                if(positions.length()) {
                    if(!CORBA::is_nil(component_ref)) {
                        if(component_ref->isReady()) {
                            component_ref->setPosition(positions, NOW);
                        }
                        else {
                            ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::setASConfiguration(): cannot set the position"));
                            ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::setASConfiguration(): SRP not ready"));
                        }
                    }
                    else {
                        ACS_SHORT_LOG((LM_WARNING, "MinorServoBossImpl::setASConfiguration(): cannot set the position"));
                    }
                }
                else {
                    ACS_SHORT_LOG((LM_WARNING, "setASConfiguration(): cannot get the position to set."));
                }
            }
        }
    }
    catch(...) {
        THROW_EX(ManagementErrors, ConfigurationErrorEx, "setASConfiguration(): cannot change the actual configuration.", false);
    }
    */
}


bool MinorServoBossImpl::isElevationTrackingEn() {
    return m_actual_config->can_track_elevation();
}


bool MinorServoBossImpl::isElevationTracking() {
    return m_servo_status.elevation_tracking;
}


bool MinorServoBossImpl::isTracking() {
    return m_control->is_tracking();
}


char * MinorServoBossImpl::getActualSetup() {
    //return CORBA::string_dup((m_actual_config->get_name()).c_str());
    return CORBA::string_dup(m_actual_conf.c_str());
}


char * MinorServoBossImpl::getCommandedSetup() {
    //return CORBA::string_dup((m_configuration->m_commandedSetup).c_str());
    return CORBA::string_dup((m_actual_config->get_name()).c_str());
}

void 
MinorServoBossImpl::connect()
throw (MinorServoErrors::CommunicationErrorExImpl)
{
    try{
        m_control->connect();
    }catch(ServoTimeoutError& ste){
        THROW_EX(MinorServoErrors,CommunicationErrorEx, ste.what(), false);
    }catch(const ServoConnectionError& sce){
        THROW_EX(MinorServoErrors, CommunicationErrorEx, sce.what(), false);
    }
}

void 
MinorServoBossImpl::disconnect()
throw (MinorServoErrors::CommunicationErrorExImpl)
{
    try{
        m_control->disconnect();
    }catch(ServoTimeoutError& ste){
        THROW_EX(MinorServoErrors,CommunicationErrorEx, ste.what(), false);
    }catch(const ServoConnectionError& sce){
        THROW_EX(MinorServoErrors, CommunicationErrorEx, sce.what(), false);
    }
}

void 
MinorServoBossImpl::reset()
throw (MinorServoErrors::CommunicationErrorExImpl)
{
    try{
        m_control->reset();
    }catch(const ServoTimeoutError& ste){
        THROW_EX(MinorServoErrors, CommunicationErrorEx, ste.what(), false);
    }catch(const ServoConnectionError& sce){
        THROW_EX(MinorServoErrors, CommunicationErrorEx, sce.what(), false);
    }
}

ACS::Time get_min_time(double range, double acceleration, double max_speed) {
    if(max_speed == 0 || acceleration == 0)
        return static_cast<ACS::Time>(pow((double)10, (double)15));
    /*
     * s = 1/2 * a * T**2
     * T = Vmax / a
     * s = 1/2 Vmax**2 / a
     * Acceleration and deceleration -> 2 * s = Vmax**2 / a
     * Time in acce and decc: 2*Vmax/acce
    */
    double threshold = pow(max_speed, 2) / acceleration; // 
    double min_time = range <= threshold ? sqrt(range / acceleration) : \
           2 * max_speed / acceleration + (range - threshold) / max_speed;
    return static_cast<ACS::Time>(min_time * pow((double)10, (double)7));
}

GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTSystemStatus, m_status, status);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_ready, ready);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, ACS::ROstring, m_actualSetup, actualSetup);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_starting, starting);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_asConfiguration, asConfiguration);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_elevationTrack, elevationTrack);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_scanActive, scanActive);
GET_PROPERTY_REFERENCE(MinorServoBossImpl, Management::ROTBoolean, m_scanning, scanning);

/* --------------- [ MACI DLL support functions ] -----------------*/
#include <maciACSComponentDefines.h>
MACI_DLL_SUPPORT_FUNCTIONS(MinorServoBossImpl)

