 // $Id: XDataThread.cpp,v 1.3 2010/06/21 11:04:51 bliliana Exp $

#include "XDataThread.h"

using namespace IRA;

XDataThread::XDataThread(const ACE_CString& name,DataParameter *par,
  const ACS::TimeInterval& responseTime,const ACS::TimeInterval& sleepTime)
{
	AUTO_TRACE("XDataThread::XDataThread()");
	m_sender=par->sender;
	groupS=par->group;
	ciclo=true;//Varialbile che mi mantiene il Loop del thread
}

XDataThread::~XDataThread()
{ 
	AUTO_TRACE("XDataThread::~XDataThread()");
}

void XDataThread::Start()
{
	AUTO_TRACE("XDataThread::onStart()");
	ciclo=true;   
	ACS_DEBUG("XDataThread::onStart()","Data thread starts");

}
	
void XDataThread::Stop()
{
	AUTO_TRACE("XDataThread::onStop()");
	ciclo=false;  
	ACS_DEBUG("XDataThread::onStop()","Data thread stops");
}

void XDataThread::run()
{
	while(ciclo){
		StartSem.acquire(1);
		groupS->GetData(groupS->Xspec.GetNCicli());
		groupS->countInt++;
		if(groupS->Xspec.GetNCicli()==0) StopSem.release(1);	
	}
}