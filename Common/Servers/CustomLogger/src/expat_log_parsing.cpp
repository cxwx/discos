#include <expat_log_parsing.h>

LogRecord::LogRecord() : _initialized(false),
                         _parsing_message_cdata(false),
                         _parsing_data_cdata(false),
                         _finished(false),
                         _depth(0)
{
};

LogRecord::~LogRecord()
{
};
 
void
LogRecord::add_data(std::string k, std::string v)
{
    kwargs.insert(KVMap::value_type(k, v));
};

std::string 
LogRecord::get_data(std::string key)
{
    KVMap::iterator finder;
    finder = kwargs.find(key);
    if(finder != kwargs.end())
        return finder->second;
    else
        return NULL;
};

/*
* Gets a const char* timestamp as extracted from the xml log record
* and returns a string without '-' and ':' characters.
*/
std::string 
trim_date(const char *str)
{
    int pos;
    std::string res(str);
    pos = res.find('-');
    while(pos != -1)
    {   
        res = res.erase(pos, 1);
        pos = res.find('-');
    }
    pos = res.find(':');
    while(pos != -1)
    {   
        res = res.erase(pos, 1);
        pos = res.find(':');
    }
    return res;
};

/*
* Convert a TimpeStamp xml attribute to ACS::Time
* In the xml logs the timestamp is formatted as YYYY-mm-ddTHH:MM:SS.mmm
*/
ACS::Time 
parse_timestamp(const char *ts)
{
    EpochHelper eh; 
    //date becomes YYYYmmddTHHMMSS.mmm
    std::string timestamp = trim_date(ts);
    try{
        eh.year(boost::lexical_cast<unsigned long>(timestamp.substr(0, 4)));
        eh.month(boost::lexical_cast<long>(timestamp.substr(4, 2)));
        eh.day(boost::lexical_cast<long>(timestamp.substr(6, 2)));
        eh.hour(boost::lexical_cast<long>(timestamp.substr(9, 2)));
        eh.minute(boost::lexical_cast<long>(timestamp.substr(11, 2)));
        eh.second(boost::lexical_cast<long>(timestamp.substr(13, 2)));
        eh.microSecond(1000 * boost::lexical_cast<long>(timestamp.substr(16, 3)));
    }catch(boost::bad_lexical_cast& blc)
    {
        ACE_ERROR ((LM_ERROR, blc.what()));
    }
    return eh.value().value;
};

/*
* Convert a LogRecord to its stringified representation
* that will be used in the produced logfile
*/
std::string 
log_to_string(const LogRecord& log_record)
{
    std::stringstream res;
    IRA::CString timestamp;
    IRA::CIRATools::timeToStr(log_record.timestamp, timestamp);
    res << (const char*)timestamp << " ";
    res << log_record.log_level_name << " ";
    res << log_record.message;
    return res.str();
};

void
start_parsing_element(void *data, const char *el, const char **attr)
{
    LogRecord *lr = (LogRecord*) data;
    if(strcmp(el, "root") != 0)
    {
            lr->_depth++;
            lr->_element.assign(el);
	    if(!lr->_initialized)
	    {
		lr->log_level_name.assign(el);
                lr->log_level = IRA::CustomLoggerUtils::string2customLogLevel(el);
		for(int i = 0; attr[i]; i += 2)
	            if (strcmp(attr[i], "TimeStamp") == 0)
		    {
			lr->timestamp = parse_timestamp(attr[i+1]);   
		    }
		lr->_initialized = true;
            }else if(strcmp(el, "Data") == 0)
                for(int i = 0; attr[i]; i += 2)
                    if (strcmp(attr[i], "Name") == 0)
                        lr->_data_name.assign(attr[i+1]);
    }
};

void 
start_parsing_cdata_element(void * data)
{
    LogRecord *lr = (LogRecord*) data;
    if(lr->_initialized)
        if(lr->_element.compare(lr->log_level_name) == 0)
            lr->_parsing_message_cdata = true;
        if(lr->_element.compare("Data") == 0)
            lr->_parsing_data_cdata = true;
};

void 
end_parsing_cdata_element(void * data)
{
    LogRecord *lr = (LogRecord*) data;
    lr->_parsing_data_cdata = false;
    lr->_parsing_message_cdata = false;
};

void
end_parsing_element(void *data, const char *el)
{
    LogRecord *lr = (LogRecord*) data;
    lr->_depth--; 
    if(lr->_depth == 0)
        lr->_finished = true;
};

void 
parsing_character_data(void* data, const char *chars, int len)
{
    LogRecord *lr = (LogRecord*) data;
    if(lr->_parsing_message_cdata){
        lr->message.assign(chars, len);
    }else if(lr->_parsing_data_cdata){
        lr->add_data(lr->_data_name, std::string(chars, len));
    }
};

/*
* Create the expat XML parser adding the defined handlers.
*/
XML_Parser 
init_log_parsing(){
    XML_Parser log_parser = XML_ParserCreate(NULL);
    XML_SetElementHandler(log_parser, start_parsing_element, end_parsing_element);
    XML_SetCdataSectionHandler(log_parser, start_parsing_cdata_element, end_parsing_cdata_element);
    XML_SetCharacterDataHandler(log_parser, parsing_character_data);
    //Added to simulate single xml document
    if(!XML_Parse(log_parser, START_TAG, START_TAG_LEN, false))
        ACE_ERROR ((LM_ERROR, XML_ErrorString(XML_GetErrorCode(log_parser))));
    return log_parser;
};

void 
free_log_parsing(XML_Parser log_parser)
{
    if(!XML_Parse(log_parser, END_TAG, END_TAG_LEN, true))
        ACE_ERROR ((LM_ERROR, XML_ErrorString(XML_GetErrorCode(log_parser))));
    XML_ParserFree(log_parser);
};

boost::shared_ptr<LogRecord> 
get_log_record(XML_Parser log_parser, const char *xml_text)
{
    boost::shared_ptr<LogRecord> log_record(new LogRecord);
    log_record->xml_text.assign(xml_text);
    XML_SetUserData(log_parser, log_record.get());
    if(!XML_Parse(log_parser, xml_text, std::strlen(xml_text), false))
        ACE_ERROR ((LM_ERROR, XML_ErrorString(XML_GetErrorCode(log_parser))));
    return log_record;
};