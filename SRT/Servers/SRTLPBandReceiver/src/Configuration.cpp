// $Id: Configuration.cpp,v 1.5 2011-06-01 18:24:44 a.orlati Exp $

#include "Configuration.h"
#include <string>

using namespace IRA;

#define _GET_DOUBLE_ATTRIBUTE(ATTRIB,DESCR,FIELD,NAME) { \
    double tmpd; \
    if (!CIRATools::getDBValue(Services,ATTRIB,tmpd,"alma/",NAME)) { \
        _EXCPT(ComponentErrors::CDBAccessExImpl,dummy,"CConfiguration::Init()"); \
        dummy.setFieldName(ATTRIB); \
        throw dummy; \
    } \
    else { \
        FIELD=tmpd; \
        ACS_DEBUG_PARAM("CConfiguration::Init()",DESCR" %lf",tmpd); \
    } \
}

#define _GET_DWORD_ATTRIBUTE(ATTRIB,DESCR,FIELD,NAME) { \
    DWORD tmpw; \
    if (!CIRATools::getDBValue(Services,ATTRIB,tmpw,"alma/",NAME)) { \
        _EXCPT(ComponentErrors::CDBAccessExImpl,dummy,"CConfiguration::Init()"); \
        dummy.setFieldName(ATTRIB); \
        throw dummy; \
    } \
    else { \
        FIELD=tmpw; \
        ACS_DEBUG_PARAM("CConfiguration::Init()",DESCR" %lu",tmpw); \
    } \
}

#define _GET_STRING_ATTRIBUTE(ATTRIB,DESCR,FIELD,NAME) { \
    CString tmps; \
    if (!CIRATools::getDBValue(Services,ATTRIB,tmps,"alma/",NAME)) { \
        _EXCPT(ComponentErrors::CDBAccessExImpl,dummy,"::CConfiguration::Init()"); \
        dummy.setFieldName(ATTRIB); \
        throw dummy; \
    } \
    else { \
        FIELD=tmps; \
        ACS_DEBUG_PARAM("CConfiguration::Init()",DESCR" %s",(const char*)tmps); \
    } \
}

#define CONFIG_PATH "DataBlock/SRTLPBandReceiver"
#define MARKTABLE_PATH CONFIG_PATH"/MarkCoefficients"
#define FEEDTABLE_PATH CONFIG_PATH"/Feeds"
#define TAPERTABLE_PATH CONFIG_PATH"/Taper"
#define DEFAULTMODE "L1L1"
#define DEFAULTMODE_PATH CONFIG_PATH"/Modes/"DEFAULTMODE


CConfiguration::CConfiguration()
{
    m_markVector = NULL;
    m_markVectorLen = 0;
    m_LBandPolarizations = NULL;
    m_PBandPolarizations = NULL;
    m_feedsTable = NULL;
    m_feedVector = NULL;
    m_taperTable = NULL;
    m_taperVector = NULL;
    m_taperVectorLen = 0;
    m_LBandRFMin = m_PBandRFMin = NULL;
    m_LBandRFMax = m_PBandRFMax = NULL;
    m_LBandIFMin = m_PBandIFMin = NULL;
    m_LBandLO = m_PBandLO = NULL;
    m_LBandIFBandwidth = m_PBandIFBandwidth = NULL;
}

CConfiguration::~CConfiguration()
{
    if (m_markTable) {
        delete m_markTable;
    }
    if (m_feedsTable) {
        delete m_feedsTable;
    }
    if (m_taperTable) {
        delete m_taperTable;
    }
    if (m_markVector) {
        delete [] m_markVector;
    }
    if (m_taperVector) {
        delete [] m_taperVector;
    }
    if (m_LBandPolarizations) {
        delete [] m_LBandPolarizations;
    }
    if (m_PBandPolarizations) {
        delete [] m_PBandPolarizations;
    }
    if (m_feedVector) {
        delete [] m_feedVector;
    }
    if (m_LBandRFMin) {
        delete [] m_LBandRFMin;
    }
    if (m_PBandRFMin) {
        delete [] m_PBandRFMin;
    }
    if (m_LBandRFMax) {
        delete [] m_LBandRFMax;
    }
    if (m_PBandRFMax) {
        delete [] m_PBandRFMax;
    }
    if (m_LBandIFMin) {
        delete [] m_LBandIFMin;
    }
    if (m_PBandIFMin) {
        delete [] m_PBandIFMin;
    }
    if (m_LBandIFBandwidth) {
        delete [] m_LBandIFBandwidth;
    }
    if (m_PBandIFBandwidth) {
        delete [] m_PBandIFBandwidth;
    }
}

void CConfiguration::init(maci::ContainerServices *Services) throw (
        ComponentErrors::CDBAccessExImpl,
        ComponentErrors::MemoryAllocationExImpl, 
        ReceiversErrors::ModeErrorExImpl
        )
{
    m_services = Services;
    IRA::CError error;
    IRA::CString field;
    WORD len;
    // read component configuration
    _GET_STRING_ATTRIBUTE("DewarIPAddress","Dewar IP address:",m_dewarIPAddress,"");
    _GET_STRING_ATTRIBUTE("LNAIPAddress","LNA IP address:",m_LNAIPAddress,"");
    _GET_STRING_ATTRIBUTE("SwitchIPAddress","Switch board IP address:",m_switchIPAddress,"");
    _GET_DWORD_ATTRIBUTE("DewarPort","Dewar port:",m_dewarPort,"");
    _GET_DWORD_ATTRIBUTE("LNAPort","LNA port:",m_LNAPort,"");
    _GET_DWORD_ATTRIBUTE("SwitchPort","Switch board port:",m_switchPort,"");
    _GET_DWORD_ATTRIBUTE("WatchDogResponseTime","Response time of watch dog thread (uSec):",m_watchDogResponseTime,"");
    _GET_DWORD_ATTRIBUTE("WatchDogSleepTime","Sleep time of the watch dog thread (uSec):",m_watchDogSleepTime,"");
    _GET_DWORD_ATTRIBUTE("LNASamplingTime","Time needed to collect LNA information from control boards (uSec):",m_LNASamplingTime,"");
    _GET_DWORD_ATTRIBUTE("RepetitionCacheTime","Log repetition filter, caching time (uSec):",m_repetitionCacheTime,"");
    _GET_DWORD_ATTRIBUTE("RepetitionExpireTime","Log repetition filter, expire time (uSec):",m_repetitionExpireTime,"");
    // now read the receiver configuration
    _GET_STRING_ATTRIBUTE("Mode","mode name:", m_mode, DEFAULTMODE_PATH);
    _GET_DWORD_ATTRIBUTE("Feeds","Number of feeds:", m_feeds, DEFAULTMODE_PATH);
    _GET_DWORD_ATTRIBUTE("IFs","Number of IFs for each feed:", m_IFs, DEFAULTMODE_PATH);

    try {
        m_LBandPolarizations = new Receivers::TPolarization[m_IFs];
        m_PBandPolarizations = new Receivers::TPolarization[m_IFs];

        m_LBandRFMin = new double[m_IFs];
        m_PBandRFMin = new double[m_IFs];

        m_LBandRFMax = new double[m_IFs];
        m_PBandRFMax = new double[m_IFs];

        m_LBandIFMin = new double[m_IFs];
        m_PBandIFMin = new double[m_IFs];

        m_LBandLO = new double[m_IFs];
        m_PBandLO = new double[m_IFs];

        m_LBandIFBandwidth = new double[m_IFs];
        m_PBandIFBandwidth = new double[m_IFs];
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl, dummy, "CConfiguration::init()");
        throw dummy;
    }

    // Set the default operating mode
    setMode(DEFAULTMODE); 

    // The noise mark
    try {
        m_markTable = new IRA::CDBTable(Services, "MarkEntry", MARKTABLE_PATH);
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl,dummy, "CConfiguration::init()");
        throw dummy;
    }
    error.Reset();
    // TODO: add the filed of the LP mark table. Inoltre la tabella per il ric LP deve avere per ogni
    // feed tutte e quatto le polarizzazioni, quindi oltre a quelle circolari (L e R) anche quelle lineari
    // (V e H)
    if (!m_markTable->addField(error, "Feed", IRA::CDataField::LONGLONG)) {
        field = "Feed";
    }
    else if (!m_markTable->addField(error, "Polarization", IRA::CDataField::STRING)) {
        field="Polarization";
    }
    else if (!m_markTable->addField(error, "Coefficients", IRA::CDataField::STRING)) {
        field = "Coefficients";
    }
    if (!error.isNoError()) {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        dummy.setFieldName((const char *)field);
        throw dummy;
    }
    if (!m_markTable->openTable(error)) {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        throw dummy;
    }

    m_markTable->First();
    len = m_markTable->recordCount();
    if(len != m_feeds * 2) { // Two channels per feed
        _EXCPT(ComponentErrors::CDBAccessExImpl, dummy, "CConfiguration::init()");
        dummy.setFieldName("MarkCoefficients table size");
        throw dummy;
    }
    try {
        m_markVector = new TMarkValue[len];
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl, dummy, "CConfiguration::init()");
        throw dummy;
    }
    for (WORD i=0; i<len; i++) {
        m_markVector[i].feed = (*m_markTable)["Feed"]->asLongLong();
        m_markVector[i].polarization = \
            (*m_markTable)["Polarization"]->asString() == "LEFT" ? Receivers::RCV_LCP:Receivers::RCV_RCP;

        std::vector<std::string> marks_str = \
            split(std::string((*m_markTable)["Coefficients"]->asString()), std::string(","));

        // Vector of coefficients (double)
        for(std::vector<std::string>::iterator iter = marks_str.begin(); iter != marks_str.end(); iter++)
            (m_markVector[i].coefficients).push_back(str2double(*iter));
        m_markTable->Next();
    }
    m_markVectorLen = len;
    m_markTable->closeTable();
    delete m_markTable;
    m_markTable = NULL;

    // The feeds
    try {
        m_feedsTable = new IRA::CDBTable(Services, "Feed", FEEDTABLE_PATH);
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl, dummy, "CConfiguration::init()");
        throw dummy;
    }
    error.Reset();
    if (!m_feedsTable->addField(error, "feedCode", IRA::CDataField::LONGLONG)) {
        field = "feedCode";
    }
    else if (!m_feedsTable->addField(error, "xOffset", IRA::CDataField::DOUBLE)) {
        field = "xOffset";
    }
    else if (!m_feedsTable->addField(error, "yOffset", IRA::CDataField::DOUBLE)) {
        field = "yOffset";
    }
    else if (!m_feedsTable->addField(error, "relativePower", IRA::CDataField::DOUBLE)) {
        field = "relativePower";
    }
    if (!error.isNoError()) {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        dummy.setFieldName((const char *)field);
        throw dummy;
    }
    if (!m_feedsTable->openTable(error))    {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        throw dummy;
    }
    m_feedsTable->First();
    if (m_feeds != m_feedsTable->recordCount()) {
        _EXCPT(ComponentErrors::CDBAccessExImpl, dummy, "CConfiguration::init()");
        dummy.setFieldName("feed table size");
        throw dummy;
    }
    len = m_feeds;
    try {
        m_feedVector = new TFeedValue[len];
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl,dummy,"CConfiguration::init()");
        throw dummy;
    }
    for (WORD i=0; i<len; i++) {
        m_feedVector[i].xOffset = (*m_feedsTable)["xOffset"]->asDouble();
        m_feedVector[i].yOffset = (*m_feedsTable)["yOffset"]->asDouble();
        m_feedVector[i].relativePower = (*m_feedsTable)["relativePower"]->asDouble();
        m_feedVector[i].code = (WORD)(*m_feedsTable)["feedCode"]->asLongLong();
        ACS_LOG(LM_FULL_INFO,
                "CConfiguration::init()",
                (
                     LM_DEBUG, 
                     "FEED_VALUE_ENTRY: %d %lf %lf %lf", 
                     m_feedVector[i].code,
                     m_feedVector[i].xOffset,
                     m_feedVector[i].yOffset,
                     m_feedVector[i].relativePower
                 )
        );
        m_feedsTable->Next();
    }
    m_feedsTable->closeTable();
    delete m_feedsTable;
    m_feedsTable = NULL;

    //The taper.....
    try {
        m_taperTable = new IRA::CDBTable(Services, "TaperEntry", TAPERTABLE_PATH);
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl, dummy, "CConfiguration::init()");
        throw dummy;
    }
    error.Reset();
    if (!m_taperTable->addField(error, "Frequency", IRA::CDataField::DOUBLE)) {
        field = "Frequency";
    }
    else if (!m_taperTable->addField(error, "Taper", IRA::CDataField::DOUBLE)) {
        field = "OutputPower";
    }
    if (!error.isNoError()) {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        dummy.setFieldName((const char *)field);
        throw dummy;
    }
    if (!m_taperTable->openTable(error))    {
        _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
        throw dummy;
    }
    m_taperTable->First();
    len = m_taperTable->recordCount();
    try {
        m_taperVector = new TTaperValue[len];
    }
    catch (std::bad_alloc& ex) {
        _EXCPT(ComponentErrors::MemoryAllocationExImpl, dummy, "CConfiguration::init()");
        throw dummy;
    }
    ACS_LOG(LM_FULL_INFO, "CConfiguration::init()", (LM_DEBUG, "TAPER_ENTRY_NUMBER: %d", len));
    for (WORD i=0; i<len; i++) {
        m_taperVector[i].frequency = (*m_taperTable)["Frequency"]->asDouble();
        m_taperVector[i].taper = (*m_taperTable)["Taper"]->asDouble();
        ACS_LOG(LM_FULL_INFO, 
                "CConfiguration::init()",
                (
                     LM_DEBUG,
                     "SYNTH_VALUE_ENTRY: %lf %lf",
                     m_taperVector[i].frequency,
                     m_taperVector[i].taper
                 )
        );
        m_taperTable->Next();
    }
    m_taperVectorLen = len;
    m_taperTable->closeTable();
    delete m_taperTable;
    m_taperTable = NULL;

}

void CConfiguration::setMode(const char * mode) throw (
        ComponentErrors::CDBAccessExImpl, 
        ReceiversErrors::ModeErrorExImpl
        ) 
{
	IRA::CString cmdMode(mode);
	cmdMode.MakeUpper();

    CString MODE_PATH((std::string(CONFIG_PATH) + std::string("/Modes/") + std::string(cmdMode)).c_str());
    IRA::CString value, token;
    IRA::CError error;

    maci::ContainerServices* Services = m_services;
    m_mode = "";
    // TODO: E inoltre, sostituire gli asterisco con il valore attuale mode attuale

    // _GET_DWORD_ATTRIBUTE("Feeds","Number of feeds:", m_feeds, MODE_PATH);
    // _GET_DWORD_ATTRIBUTE("LBandIFs","Number of IFs per feed:", m_IFs, MODE_PATH);

    _GET_STRING_ATTRIBUTE("LBandPolarization", "LBand IF polarization:", value, MODE_PATH);
    int start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("LBandPolarization");
            throw dummy;
        }
        token.MakeUpper();
        if (token == "L") {
            m_LBandPolarizations[k] = Receivers::RCV_LCP;
        }
        else if (token == "R") {
            m_LBandPolarizations[k] = Receivers::RCV_RCP;
        }
        else if (token == "V") {
            m_LBandPolarizations[k] = Receivers::RCV_VLP;
        }
        else if (token == "H") {
            m_LBandPolarizations[k] = Receivers::RCV_HLP;
        }
        else {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("LBandPolarization");
            throw dummy;
        }
    }

    _GET_STRING_ATTRIBUTE("PBandPolarization", "PBand IF polarization:", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("PBandPolarization");
            throw dummy;
        }
        token.MakeUpper();
        if (token == "L") {
            m_PBandPolarizations[k] = Receivers::RCV_LCP;
        }
        else if (token == "R") {
            m_PBandPolarizations[k] = Receivers::RCV_RCP;
        }
        else if (token == "V") {
            m_PBandPolarizations[k] = Receivers::RCV_VLP;
        }
        else if (token == "H") {
            m_PBandPolarizations[k] = Receivers::RCV_HLP;
        }
        else {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("PBandPolarization");
            throw dummy;
        }
    }

    _GET_STRING_ATTRIBUTE("LBandRFMin", "L band RF lower limit (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("LBandRFMin");
            throw dummy;
        }
        m_LBandRFMin[k] = token.ToDouble();
    }
    _GET_STRING_ATTRIBUTE("LBandRFMax", "L band RF upper limit (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("LBandRFMax");
            throw dummy;
        }
        m_LBandRFMax[k] = token.ToDouble();
    }

    for (WORD k=0; k<m_IFs; k++)
        m_LBandIFBandwidth[k] = m_LBandRFMax[k] - m_LBandRFMin[k];
    

    _GET_STRING_ATTRIBUTE("PBandRFMin", "P band RF lower limit (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("PBandRFMin");
            throw dummy;
        }
        m_PBandRFMin[k] = token.ToDouble();
    }
    _GET_STRING_ATTRIBUTE("PBandRFMax", "P band RF upper limit (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("PBandRFMax");
            throw dummy;
        }
        m_PBandRFMax[k] = token.ToDouble();
    }

    for (WORD k=0; k<m_IFs; k++)
        m_PBandIFBandwidth[k] = m_PBandRFMax[k] - m_PBandRFMin[k];

    _GET_STRING_ATTRIBUTE("LBandIFMin","L band IF start frequency (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("LBandIFMin");
            throw dummy;
        }
        m_LBandIFMin[k] = token.ToDouble();
    }

    _GET_STRING_ATTRIBUTE("PBandIFMin","P band IF start frequency (MHz):", value, MODE_PATH);
    start = 0;
    for (WORD k=0; k<m_IFs; k++) {
        if (!IRA::CIRATools::getNextToken(value, start, ' ', token)) {
            _EXCPT_FROM_ERROR(ComponentErrors::CDBAccessExImpl, dummy, error);
            dummy.setFieldName("PBandIFMin");
            throw dummy;
        }
        m_PBandIFMin[k] = token.ToDouble();
    }

    for (WORD k=0; k<m_IFs; k++)
        m_LBandLO[k] = m_LBandRFMin[k] - m_LBandIFMin[k];

    for (WORD k=0; k<m_IFs; k++)
        m_PBandLO[k] = m_PBandRFMin[k] - m_PBandIFMin[k];

    // TODO: Impostare il filtro chiamando il corrispondente metodo della receiver Library
    // TODO: Impostare la polarizzazione, con metodo della ReceiverLibrary

    m_mode = cmdMode;
}


DWORD CConfiguration::getTaperTable(double * &freq,double *&taper) const
{
    freq = new double [m_taperVectorLen];
    taper = new double [m_taperVectorLen];
    for (DWORD j=0; j<m_taperVectorLen; j++) {
        freq[j] = m_taperVector[j].frequency;
        taper[j] = m_taperVector[j].taper;
    }
    return m_taperVectorLen;
}

DWORD CConfiguration::getFeedInfo(
        WORD *& code,
        double *& xOffset,
        double *& yOffset,
        double *& relativePower
        ) const
{
    code = new WORD[m_feeds];
    xOffset = new double [m_feeds];
    yOffset = new double [m_feeds];
    relativePower = new double [m_feeds];
    for (DWORD j=0; j<m_feeds; j++) {
        code[j] = m_feedVector[j].code;
        xOffset[j] = m_feedVector[j].xOffset;
        yOffset[j] = m_feedVector[j].yOffset;
        relativePower[j] = m_feedVector[j].relativePower;
    }
    return m_feeds;
}
