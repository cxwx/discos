// $Id: coordinateGrabber.cpp,v 1.5 2011-06-21 16:38:02 a.orlati Exp $
#ifndef _COORDINATEGRABBER_CPP_
#define _COORDINATEGRABBER_CPP_

#include <AntennaBossC.h>
#include <ObservatoryC.h>
#include <IRA>
#include <maciSimpleClient.h>
#include <getopt.h>
#include <ClientErrors.h>
#include <iostream>
#include <fstream>
#include <Site.h>
#include <DateTime.h>
#include <slamac.h>

using namespace IRA;

bool done;

void printHelp() {
	printf("Saves observed Ra/Dec and Az/El in a file\n");
	printf("\n");
	printf("CoordinateGrabber [-h] [-s sample] -b bossInstance -o observatoryInstance FileName \n");
	printf("\n");
	printf("[-h] prints this help\n");
	printf("[-s sample] set the sample interval in milliseconds. It must be included between 10 and 999.\n");
	printf("-b  bossInstance select the name of the Antenna Boss instance to be used. This switch is mandatory.\n");
	printf("-o  observatoryInstance select the name of the Observatory instance to be used. This switch is mandatory.\n");
	printf("\n");
}

void handler(int p) 
{
	done=true;
}

int main(int argc, char *argv[]) 
{
	int out;
	maci::SimpleClient client;
	CString name,observ;
	ofstream file;
	ACE_Time_Value tv(3);
	Antenna::AntennaBoss_var boss;
	Antenna::Observatory_var observatory;
	long sampleInterval=200;
	CString fileName;
	char output[512];
	TIMEVALUE now,last,lst;
	double az,el;
	double ra,decl;
	double lon,lat;
	double dut1;
	CSite site;
	Antenna::TSiteInformation_var siteInfo;
	CDateTime ut;
	if (argc==1) {
		printf("Arguments are missing...\n\n");
		printHelp();
		exit(-1);
	}
	name="";
	while ((out=getopt(argc,argv,"hb:o:s:"))!=-1) {
		switch (out) {
			case 'h': {
				printHelp();
				exit(0);
			}
			case 'b': {
				name=CString(optarg);
				break;
			}
			case 'o': {
				observ=CString(optarg);
				break;
			}
			case 's': {
				if (sscanf(optarg,"%ld",&sampleInterval) !=1) {
					printf("Error in sample interval specification...\n\n");
					printHelp();
					exit(-1);
				}
				else {
					if ((sampleInterval<10) || (sampleInterval>1000)) {
						printf("Sample interval out of ranges...\n\n");
						printHelp();
						exit(-1);						
					}
					sampleInterval*=1000;
				}
				break;
			}
			default : {
				printf("Unknown switch.......\n");
				printHelp();
				exit(-1);
			}
		}
	}
	if (name=="") {
		printf("Antenna Boss instance must be provided...\n\n");
		printHelp();
		exit(-1);
	}
	if (observ=="") {
		printf("Observatory instance must be provided...\n\n");
		printHelp();
		exit(-1);		
	}
	fileName=CString(argv[optind]);
	try {
		if (client.init(argc,argv)==0) {
			_EXCPT(ClientErrors::CouldntInitExImpl,impl,"coordinateGrabber:main()");
			impl.log();
			exit(-1);
		}
		else {
			if (client.login()==0) {
				_EXCPT(ClientErrors::CouldntLoginExImpl,impl,"coordinateGrabber::main()");
				impl.log();
				exit(-1);
			}
		}	
	}
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,"coordinateGrabber::main()");
		impl.log();
		exit(-1);
	}
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"LOGGED_IN"));	
	try {
		boss=client.getComponent<Antenna::AntennaBoss>((const char *)name,0,true);
		if (CORBA::is_nil(boss.in())==true) {
			_EXCPT(ClientErrors::CouldntAccessComponentExImpl,impl,"coordinateGrabber::main()");
			impl.setComponentName((const char *)name);
			impl.log();
			exit(-1);
		}
	}
	catch(CORBA::SystemException &E) {
		_EXCPT(ClientErrors::CORBAProblemExImpl,impl,"coordinateGrabber::main()");
		impl.setName(E._name());
		impl.setMinor(E.minor());
		impl.log();
		exit(-1);
	}
	catch (maciErrType::CannotGetComponentExImpl& E) {
		_ADD_BACKTRACE(ClientErrors::CouldntAccessComponentExImpl,impl,E,"coordinateGrabber::main()");
		impl.setComponentName((const char *)name);
		impl.log();
		exit(-1);
	}	
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,"coordinateGrabber::main()");
		impl.log();
		exit(-1);
	}
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"GOT_COMPONENENT: %s",(const char *)name));
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_DEBUG,"Reference is: %d",boss.ptr()));
	try {
		observatory=client.getComponent<Antenna::Observatory>((const char*)observ,0,true);
		if (CORBA::is_nil(observatory.in())==true) {
			_EXCPT(ClientErrors::CouldntAccessComponentExImpl,impl,"coordinateGrabber::main()");
			impl.setComponentName((const char *)observ);
			impl.log();
			exit(-1);
		}
	}
	catch(CORBA::SystemException &E) {
		_EXCPT(ClientErrors::CORBAProblemExImpl,impl,"coordinateGrabber::main()");
		impl.setName(E._name());
		impl.setMinor(E.minor());
		impl.log();
		exit(-1);
	}
	catch (maciErrType::CannotGetComponentExImpl& E) {
		_ADD_BACKTRACE(ClientErrors::CouldntAccessComponentExImpl,impl,E,"coordinateGrabber::main()");
		impl.setComponentName((const char *)observ);
		impl.log();
		exit(-1);
	}	
	catch(...) {
		_EXCPT(ClientErrors::UnknownExImpl,impl,"coordinateGrabber::main()");
		impl.log();
		exit(-1);
	}	
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"GOT_COMPONENENT: %s",(const char *)observ));
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_DEBUG,"Reference is: %d",observatory.ptr()));
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"ALL_COMPONENTS_RETRIEVED"));
	ACE_OS::sleep(1);	
	try	{	
		siteInfo=observatory->getSiteSummary();
	}
	catch (CORBA::SystemException& ex)	{
		_EXCPT(ClientErrors::CORBAProblemExImpl,__dummy,"coordinateGrabber::main()");
		__dummy.setName(ex._name());
		__dummy.setMinor(ex.minor());
		__dummy.log();
		exit(-1);
	}
	site=CSite(siteInfo.out());
	dut1=siteInfo->DUT1;
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"FILE_CREATION"));
	file.open((const char *)fileName,ios_base::out|ios_base::trunc);
	if (!file.is_open()) {
		ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"FILE_CANNOT_BE_CREATED"));
		exit(-1);
	}
	file<<"UT\tSIDEREAL\tRA2000(deg)\tDEC2000(deg)\tAzimuth(deg)\tElevation(deg)\tGLongitude(deg)\tGLatitude(deg)\n";
	tv.set(0,sampleInterval);
	signal(SIGINT,handler);
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"RUNNING"));
	printf("Ctrl+c to terminate\n");
	CIRATools::getTime(last);
	done=false;
	while(!done) {
		IRA::CIRATools::getTime(now);
		ut.setDateTime(now,dut1);
		try {
			boss->getObservedHorizontal(now.value().value,0,az,el);
			boss->getObservedEquatorial(now.value().value,0,ra,decl);
			boss->getObservedGalactic(now.value().value,lon,lat);
		}
		catch (ACSErr::ACSbaseExImpl& ex) {
			_ADD_BACKTRACE(ClientErrors::CouldntPerformActionExImpl,impl,ex,"coordinateGrabber::main()"); 
			impl.log();
		}
		catch (CORBA::SystemException& ex) {
			_EXCPT(ClientErrors::CORBAProblemExImpl,impl,"coordinateGrabber::main()");
			impl.setName(ex._name());
			impl.setMinor(ex.minor());
			impl.log();
		}
		catch (...) {
			_EXCPT(ClientErrors::UnknownExImpl,impl,"coordinateGrabber::main()");
			impl.log();
		}
		ut.LST(site).getDateTime(lst);
		sprintf(output,"%02d:%02d:%02d.%03d\t%02d:%02d:%02d.%03d\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf\n",
				now.hour(),now.minute(),now.second(),now.microSecond()/1000,
				lst.hour(),lst.minute(),lst.second(),lst.microSecond()/1000,
				ra*DR2D,decl*DR2D,az*DR2D,el*DR2D,lon*DR2D,lat*DR2D);
		file << output;
		if (CIRATools::timeDifference(now,last)>1000000) {
			ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"FLUSHING"));
			file.flush();
			CIRATools::timeCopy(last,now);
		}
		tv.set(0,sampleInterval);
		client.run(tv);
	}
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"TERMINATING"));
	try {
		client.releaseComponent((const char *)name);
		client.releaseComponent((const char *)observ);
	}
	catch (maciErrType::CannotReleaseComponentExImpl& E) {
		E.log();
	}
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"COMPONENT_RELEASED"));
	CIRATools::Wait(1,0);
	client.logout();
	ACS_LOG(LM_FULL_INFO,"coordinateGrabber::main()",(LM_INFO,"LOGGED_OUT"));
	CIRATools::Wait(1,0);
	file.close();
	CIRATools::Wait(1,0);
	signal(SIGINT,SIG_DFL);
	return 0;
}

#endif /*COORDINATEGRABBER_CPP_*/
