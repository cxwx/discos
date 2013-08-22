#include "SchedulerTuiClient.h"

/*#define _GET_ACS_PROPERTY(TYPE,NAME) TYPE##_var NAME; \
{	\
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,"Trying to get property "#NAME"...")); \
	NAME=component->NAME(); \
	if (NAME.ptr()!=TYPE::_nil()) { \
		ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_DEBUG,"OK, reference is: %x",NAME.ptr())); \
	} \
	else { \
		_EXCPT(ClientErrors::CouldntAccessPropertyExImpl,impl,MODULE_NAME"::Main()"); \
		impl.setPropertyName(#NAME); \
		impl.log(); \
		goto ErrorLabel; \
	} \
}

#define _INSTALL_MONITOR(COMP,TRIGGERTIME) { \
	if (!COMP->installAutomaticMonitor(GUARDINTERVAL)) { \
		_EXCPT(ClientErrors::CouldntPerformActionExImpl,impl,MODULE_NAME"::Main()"); \
		impl.setAction("Install monitor"); \
		impl.setReason((const char*)COMP->getLastError()); \
		impl.log(); \
		ACE_OS::sleep(1); \
		goto ErrorLabel; \
	} \
	COMP->setTriggerTime(TRIGGERTIME); \
}

#define _GET_PROPERTY_VALUE_ONCE(TYPE,FORMAT,PROPERTY,PROPERTNAME,CONTROL) { \
	try { \
		TYPE temp; \
		ACSErr::Completion_var cmpl; \
		temp=PROPERTY->get_sync(cmpl.out()); \
		CompletionImpl cmplImpl(cmpl.in()); \
		if (cmplImpl.isErrorFree()) { \
			IRA::CString val(0,FORMAT,temp); \
			CONTROL->setValue(val); \
			CONTROL->Refresh(); \
		} \
		else { \
			_ADD_BACKTRACE(ClientErrors::CouldntAccessPropertyExImpl,impl,cmplImpl,"::Main()"); \
			impl.setPropertyName(PROPERTNAME); \
			impl.log(); \
			goto ErrorLabel; \
		} \
	} \
	catch (...) { \
		_EXCPT(ClientErrors::UnknownExImpl,impl,"::Main()"); \
		impl.log(); \
		goto ErrorLabel; \
	} \
}

#define _CATCH_ALL(OUTPUT,ROUTINE,COMP_EXCEPTION) \
	catch (COMP_EXCEPTION &E) { \
		_ADD_BACKTRACE(ClientErrors::CouldntPerformActionExImpl,impl,E,ROUTINE); \
		IRA::CString Message; \
		if (OUTPUT!=NULL) { \
			_EXCPT_TO_CSTRING(Message,impl); \
			OUTPUT->setValue(Message); \
			OUTPUT->Refresh(); \
		} \
		impl.log(); \
	} \
	catch (CORBA::SystemException &C) { \
		_EXCPT(ClientErrors::CORBAProblemExImpl,impl,ROUTINE); \
		impl.setName(C._name()); \
		impl.setMinor(C.minor()); \
		if (OUTPUT!=NULL) { \
			IRA::CString Message; \
			_EXCPT_TO_CSTRING(Message,impl); \
			OUTPUT->setValue(Message); \
			OUTPUT->Refresh(); \
		} \
		impl.log(); \
	} \
	catch (...) { \
		_EXCPT(ClientErrors::UnknownExImpl,impl,ROUTINE); \
		if (OUTPUT!=NULL) { \
			IRA::CString Message; \
			_EXCPT_TO_CSTRING(Message,impl); \
			OUTPUT->setValue(Message); \
			OUTPUT->Refresh(); \
		} \
		impl.log(); \
	}
*/
#define COMPONENT_INTERFACE COMPONENT_IDL_MODULE::COMPONENT_IDL_INTERFACE
#define COMPONENT_DECLARATION COMPONENT_IDL_MODULE::COMPONENT_SMARTPOINTER

using namespace TW;

static bool terminate;

void quintHandler(int sig)
{
	terminate=true;
}

int main(int argc, char *argv[]) {
	bool loggedIn=false;
	int fieldCounter;
	maci::SimpleClient client;
	ACE_Time_Value tv(1);
	IRA::CString fields[MAXFIELDNUMBER];
	char formatString[20];
	IRA::CString inputCommand;
	
	// Component declaration 
	COMPONENT_DECLARATION component;
		
	
	// Add frame controls declaration
	//TW::CPropertyText<_TW_PROPERTYCOMPONENT_T_RO(uLongLong)> *universalTime_field;
	// *******************************
	TW::CLabel *output_label;
	TW::CInputCommand *userInput;

	// Add callbacks object declaration
	//ACS::CBDescIn desc;
	//CCallbackVoid cbUnstow;
	//CCallbackVoid cbStow;
	// ********************************
	terminate=false;
	// mainframe 
	CFrame window(CPoint(0,0),CPoint(WINDOW_WIDTH,WINDOW_HEIGHT),'|','|','-','-');
	
	// disable ctrl+C
	signal(SIGINT,SIG_IGN);
	signal(SIGUSR1,quintHandler);
	
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::MANAGER_LOGGING"));
	/*try {
		if (client.init(argc,argv)==0) {
			_EXCPT(ClientErrors::CouldntInitExImpl,impl,MODULE_NAME"::Main()");
			impl.log();		
			goto ErrorLabel;
		}
		else {
			if (client.login()==0) {
				_EXCPT(ClientErrors::CouldntLoginExImpl,impl,MODULE_NAME"::Main()");
				impl.log();		
				goto ErrorLabel;
			}
			loggedIn=true;
		}	
	}
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,MODULE_NAME"::Main()");
		IRA::CString Message;
		impl.log();
		goto ErrorLabel;
	}*/
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::DONE"));	
	// Window initialization
	window.initFrame();
	window.setTitle(APPLICATION_TITLE);
	// Change the style of the main frame 
	window.setTitleStyle(CStyle(TITLE_COLOR_PAIR,TITLE_STYLE));
	
	if (window.colorReady()) {
		window.defineColor(GRAY_COLOR,255,255,255);		
		window.defineColorPair(BLUE_GRAY,CColor::BLUE,GRAY_COLOR);
		window.defineColorPair(GRAY_BLUE,GRAY_COLOR,CColor::BLUE);		
	}
	else {
		window.defineColorPair(BLUE_GRAY,CColor::BLUE,CColor::WHITE);
		window.defineColorPair(GRAY_BLUE,CColor::WHITE,CColor::BLUE);	
		window.defineColorPair(BLACK_RED,CColor::BLACK,CColor::RED);	
		window.defineColorPair(BLACK_GREEN,CColor::BLACK,CColor::GREEN);	
		window.defineColorPair(BLACK_YELLOW,CColor::BLACK,CColor::YELLOW);	
		window.defineColorPair(BLACK_BLUE,CColor::BLACK,CColor::BLUE);		
		window.defineColorPair(BLACK_MAGENTA,CColor::BLACK,CColor::MAGENTA);	
		window.defineColorPair(BLACK_WHITE,CColor::BLACK,CColor::WHITE	);			
	}
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::FRAME_INITIALIZED"));
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::GET_COMPONENENT: %s",COMPONENT_NAME));
	/*try {
		component=client.getComponent<COMPONENT_INTERFACE>(COMPONENT_NAME,0,true);
		if (CORBA::is_nil(component.in())==true) {
			_EXCPT(ClientErrors::CouldntAccessComponentExImpl,impl,MODULE_NAME"::Main()");
			impl.setComponentName(COMPONENT_NAME);
			impl.log();
			goto ErrorLabel;
		}
	}
	catch(CORBA::SystemException &E) {
		_EXCPT(ClientErrors::CORBAProblemExImpl,impl,MODULE_NAME"::Main()");
		impl.setName(E._name());
		impl.setMinor(E.minor());
		impl.log();
		goto ErrorLabel;
	}
	catch (maciErrType::CannotGetComponentExImpl& E) {
		_ADD_BACKTRACE(ClientErrors::CouldntAccessComponentExImpl,impl,E,MODULE_NAME"::Main()");
		impl.setComponentName(COMPONENT_NAME);
		impl.log();
		goto ErrorLabel;
	}	
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,MODULE_NAME"::Main()");
		impl.log();
		goto ErrorLabel;
	}		
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_DEBUG,"OK, reference is: %d",component.ptr()));
	ACE_OS::sleep(1);*/
	try {
		/*
		// Add all component properties here
		//_GET_ACS_PROPERTY(ACS::ROuLongLong,universalTime);
		// *********************************
		ACE_OS::sleep(1);
		
		// Frame controls creation
		//universalTime_field=new TW::CPropertyText<_TW_PROPERTYCOMPONENT_T_RO(uLongLong)>(universalTime.in());	
		//
		*/
		#if USE_OUTPUT_FIELD >=1 
			output_label=new TW::CLabel("");
		#else
			output_label=NULL;
		#endif
		userInput=new TW::CInputCommand();
		
		// setting up the properties of the components of the frame controls
		//_TW_SET_COMPONENT(universalTime_field,22,0,22,1,CColorPair::WHITE_BLACK,CStyle::BOLD,output_label);
		//universalTime_field->setFormatFunction(CFormatFunctions::dateTimeClockFormat,NULL);
		// ******************************************************************
		_TW_SET_COMPONENT(userInput,0,WINDOW_HEIGHT-6,WINDOW_WIDTH-1,1,USER_INPUT_COLOR_PAIR,USER_INPUT_STYLE,NULL);
		#if USE_OUTPUT_FIELD >=1 
			_TW_SET_COMPONENT(output_label,0,WINDOW_HEIGHT-(OUTPUT_FIELD_HEIGHT+1),WINDOW_WIDTH-1,OUTPUT_FIELD_HEIGHT,OUTPUT_FIELD_COLOR_PAIR,OUTPUT_FIELD_STYLE,NULL);	
		#endif 
		
		// Monitors
		ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::MONITORS_INSTALLATION"));
		// Add all required monitor installation here
		//_INSTALL_MONITOR(universalTime_field,200);
		// ******************************************
		ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::DONE"));
		
		// Add all the static labels
		//_TW_ADD_LABEL("Universal Time   :",0,0,18,1,CColorPair::WHITE_BLACK,CStyle::UNDERLINE,window);
		// *************************
		
		// Add all required association: components/Frame
		//window.addComponent((CFrameComponent*)universalTime_field);
		// **********************************************
		window.addComponent((CFrameComponent*)userInput);		
		#if USE_OUTPUT_FIELD >=1 
			window.addComponent((CFrameComponent*)output_label);
		#endif

	}
	catch(CORBA::SystemException &E) {		
		_EXCPT(ClientErrors::CORBAProblemExImpl,impl,MODULE_NAME"::Main()");
		impl.setName(E._name());
		impl.setMinor(E.minor());
		impl.log();
		goto ErrorLabel;
	}
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,MODULE_NAME"::Main()");
		IRA::CString Message;
		impl.log();
		goto ErrorLabel;
	}
	// now it is possible to show the frame
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::START"));
	window.showFrame();
	while(!terminate)	{
		if (userInput->readCommand(inputCommand)) {
			if (inputCommand=="exit") terminate=true;
			/*if (component->_is_a("IDL:alma/Management/CommandInterpreter:1.0")) {
				try {
					char * outputAnswer;
					component->command((const char *)inputCommand,outputAnswer);
					output_label->setValue(outputAnswer);
					CORBA::string_free(outputAnswer);
					output_label->Refresh();
				}
				catch (CORBA::SystemException& ex) {
					_EXCPT(ClientErrors::CORBAProblemExImpl,impl,"Main()");
					impl.setName(ex._name());
					impl.setMinor(ex.minor());
					IRA::CString Message;
					_EXCPT_TO_CSTRING(Message,impl);
					output_label->setValue(Message);
					output_label->Refresh();					
					impl.log(LM_ERROR);
				}
				catch(...) {
					_EXCPT(ClientErrors::CouldntPerformActionExImpl,impl,"Main()");
					impl.setAction("command()");
					impl.setReason("communication error to component server");
					IRA::CString Message;
					_EXCPT_TO_CSTRING(Message,impl);
					output_label->setValue(Message);
					output_label->Refresh();
	                impl.log(LM_ERROR); 					
				}
			}
			else {
				output_label->setValue("not enabled");
				output_label->Refresh();				
			}*/
		}
		/*client.run(tv);
		tv.set(0,MAINTHREADSLEEPTIME*1000);*/
		ACE_OS::sleep(1);
	}
		
	window.closeFrame();
	
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::RELEASING"));
	goto CloseLabel;
ErrorLabel:
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::ABORTING"));	
CloseLabel:
	window.Destroy();
	/*ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::FRAME_CLOSED"));
	ACE_OS::sleep(1);	
	try {
		client.releaseComponent(COMPONENT_NAME);
	}
	catch (maciErrType::CannotReleaseComponentExImpl& E) {
		E.log();
	}
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::COMPONENT_RELEASED"));
	if (loggedIn) client.logout();
	ACS_LOG(LM_FULL_INFO,MODULE_NAME"::Main()",(LM_INFO,MODULE_NAME"::SHUTTING_DOWN"));		
	signal(SIGINT,SIG_DFL);
	ACE_OS::sleep(1);*/
	return 0;
}

