#include <maciSimpleClient.h>
#include <RefractionC.h>
#include <RefractionImpl.h>
#include <iostream>
#include <stdio.h>
#include <slamac.h>

using namespace std;
using namespace maci;

int main (int argc, char *argv[])
{

	SimpleClient *client = new SimpleClient;
	ACSErr::Completion_var completion;
	Antenna::Refraction_var Refraction;
	ofstream m_file;
	IRA::CString out;
        IRA::CString fileName = "stringtest";

	m_file.open ((const char *) fileName, ios_base::out | ios_base::app);

  	/**
   	* Create the instance of Client and Init() it.
  	*/
    	ACS_SHORT_LOG((LM_INFO, "Activating client: %s","TestClient"));
  	if (client->init (argc, argv) == 0)
    	{
      		ACS_SHORT_LOG ((LM_ERROR, "TestClient::main: Cannot init client."));
      		return -1;
    	}
  	else
    	{
      		client->login();
    	}
    	ACS_SHORT_LOG((LM_INFO, "OK activation client: %s","TestClient"));
	ACE_OS::sleep (2); // wait two seconds


	//Get reference to the component

	ACS_SHORT_LOG((LM_INFO, "Getting component: %s","ANTENNA/Refraction"));
  	try
  	{
    		Refraction = client->get_object < Antenna::Refraction > ("ANTENNA/Refraction", 0, true);
    		if (CORBA::is_nil(Refraction.in()) == true)
      		{
			ACS_SHORT_LOG ((LM_ERROR, "TestClient::main: Failed to get a reference to Refraction component."));
			return -1;
      		}
   	}
   	catch(...)
	{
		ACS_SHORT_LOG((LM_ERROR,"Error!"));
		return -1;
   	}

	ACS_SHORT_LOG((LM_INFO, "OK activation component: %s","Refraction"));

	double zenithDist, zenithCorrected;
	//double temperature = 20.0;	// Metrology.getTemperature();
	//double humidity = 0.5;	// Metrology.getHumidity();
	//double pressure = 1013.0;	// Metrology.getPressure();
	//double i = 0.0;
	for (unsigned i = 10;i<= 90;i+=10) {
		zenithDist=90.0- i;
		//zenithCorrected = 90 - i*5.0;
		printf("Elevation = %f\n", (double)i);
//		Refraction->setMeteoParameters(20-i, 0.5+i/100., 1013+i);
		Refraction->getCorrection(zenithDist*DD2R,0.05,zenithCorrected); //@5 cm
		//printf ("zenithCorrected = %f arcseconds\n", zenithCorrected*DR2AS);
		//out.Format("zenithCorrected = %f arcseconds\n", zenithCorrected*DR2AS);
		//std::cout << (const char *) out << std::endl;
		//m_file << (const char *) out;
		printf ("elevation = %f, elevation corrected = %f\n",(double) i, zenithCorrected*DR2D);
		out.Format("elevation = %f, elevation corrected = %f\n", (double)i,zenithCorrected*DR2D);
		m_file << (const char *) out;
		std::cout << (const char *) out << std::endl;
		ACE_OS::sleep (1);
	}

	ACE_OS::sleep (3);
	ACS_SHORT_LOG ((LM_INFO, "Releasing..."));
  	client->manager ()->release_component (client->handle (), "ANTENNA/Refraction");
  	client->logout ();
  	delete client;

	//Sleep for 3 sec to allow everytihng to cleanup and stabilize
	ACE_OS::sleep (3);

	return 0;
}
