#ifndef _BASECCONFIGURATION_H_
#define _BASECCONFIGURATION_H_

/*******************************************************************************\
 IRA Istituto di Radioastronomia                                 
 This code is under GNU General Public License (GPL).          
                                                              
 Author: Andrea Orlati (aorlati@ira.inaf.it)
 Modified by Marco Buttu (mbuttu@oa-cagliari.inaf.it)
\*******************************************************************************/


#include <IRA>
#include <maciContainerServices.h>
#include <ComponentErrors.h>
#include <ReceiversErrors.h>
#include <ReceiversDefinitionsC.h>
#include "utils.h"

// Dummy value for a board connection error
#define CEDUMMY 100000


/**
 * This class implements the component configuration. The data inside this class are initialized at the startup from the
 * configuration database and then are used (read) inside the component.
 * @author <a href=mailto:a.orlati@ira.cnr.it>Andrea Orlati</a>, Istituto di Radioastronomia, Italia
 * @author <a href=mailto:m.buttu@oa-cagliari.inaf.it>Maro Buttu</a>, Osservatorio Astronomico di Cagliari
 * <br> 
  */
class CConfiguration {
public:

	typedef struct {
		double skyFrequency;
		long long feed; 
		Receivers::TPolarization polarization;
		double markValue;
	} TMarkValue;

	typedef struct {
		double frequency;
		double outputPower;
	} TLOValue;

	typedef struct {
		long long feed; 
		double frequency;
		double taper;
	} TTaperValue;

	typedef struct {
		WORD code;
		double xOffset;
		double yOffset;
		double relativePower;
	} TFeedValue;

	typedef struct {
		double temperature;
        ACS::Time timestamp;
	} BoardValue;


	/**
	 * Default constructor
	 */
	CConfiguration();

	/**
	 * Destructor
	 */
	~CConfiguration();

	/**
	 * @return the IP address of the board that controls the dewar
	 */
	inline const IRA::CString& getDewarIPAddress() const { return m_dewarIPAddress; }

	/**
	 * @return the IP port of the board that controls the dewar
	 */
	inline const WORD& getDewarPort() const { return m_dewarPort; }

	/**
	 * @return the IP address of the board that controls the receiver LNA
	 */
	inline const IRA::CString& getLNAIPAddress() const { return m_LNAIPAddress; }

	/**
	 * @return the IP port of the board that controls the receiver LNA
	 */
	inline const WORD& getLNAPort() const { return m_LNAPort; }

	/**
	 * @return the IP address of the board that controls the filter and the polarization
	 */
	inline const IRA::CString& getSwitchIPAddress() const { return m_switchIPAddress; }

	/**
	 * @return the IP port of the board that controls the receiver filter and polarization
	 */
	inline const WORD& getSwitchPort() const { return m_switchPort; }

	/**
	 * @return the time allowed to the watch dog thread to complete an iteration and 
     * respond to the thread manager (microseconds)
	 */
	inline const DDWORD& getWarchDogResponseTime() const { return m_watchDogResponseTime; }

	/**
	 * @return the time between two execution of the watch dog thread (microseconds)
	 */
	inline const DDWORD& getWatchDogSleepTime() const { return m_watchDogSleepTime; }

	/**
	 * @return the time required by the board to sample the values of LNA (microseconds)
	 */
	inline const DDWORD& getLNASamplingTime() const { return m_LNASamplingTime; }

	/**
	 * @return the time the repetition log guard will cache new log messages before sending 
     * to the central logger (microseconds)
	 */
	inline const DDWORD& getRepetitionCacheTime() const { return m_repetitionCacheTime; }

	/**
	 * @return the time of expiration of a log sent to the repetition log guard (microseconds)
	 */
	inline const DDWORD& getRepetitionExpireTime() const { return m_repetitionExpireTime; }

	/**
	 * Allows to get the table of mark values relative to left polarization
	 * @param freq vector containing the frequency value of the mark table. It must be freed by caller.
	 * @param markValue vector of the value of the calibration diode. It must be freed by caller.
	 * @param len used to return the length of the mark values array
	 * @return the size of the output vectors
	 */
	DWORD getLeftMarkTable(double *& freq,double *& markValuel, short feed) const;

	/**
	 * Allows to get the table of mark values relative to left polarization
	 * @param freq vector containing the frequency value of the mark table. It must be freed by caller.
	 * @param markValue vector of the value of the calibration diode. It must be freed by caller.
	 * @param len used to return the length of the mark values array
	 * @return the size of the output vectors
	 */
	DWORD getRightMarkTable(double *& freq,double *& markValue, short feed) const;

	/**
	 * @return the instance of the local oscillator component  that the receiver will use to drive the set its LO
	 */
	inline const IRA::CString& getLocalOscillatorInstance() const { return m_localOscillatorInstance; }

	/**
	 * @param freq vector with the synthesizer frequencies. It must be freed by caller.
	 * @param power corresponding powers for the frequencies vector. It must be freed by caller.
	 * @return the size of the output vectors
	 */
	DWORD getSynthesizerTable(double *& freq,double *& power) const;

	/**
	 * @param freq vector with the frequencies. It must be freed by caller.
	 * @param taper corresponding taper  for the frequencies vector. It must be freed by caller.
	 * @return the size of the output vectors
	 */
	DWORD getTaperTable(double * &freq,double *&taper, short feed) const;

	/**
	 * @param code feed identification codes. It must be freed by caller.
	 * @param xOffset displacement of the feed with respect to the central one along x axis. 
     * It must be freed by caller.
	 * @param yOffset displacement of the feed with respect to the central one along y axis. 
     * It must be freed by caller.
	 * @param relativePower expected percentage of variation of gain with respect to the central one. 
     * It must be freed by caller.
	 * @return the size of the output vectors
	 */
	DWORD getFeedInfo(WORD *& code,double *& xOffset,double *& yOffset,double *& relativePower) const;


	/**
	 * @return mnemonic of the working mode of the receiver
	 */
	inline const IRA::CString& getActualMode() const { return m_mode; }


	/**
	 * @return mnemonic of the default mode of the receiver
	 */
	inline const IRA::CString& getDefaultMode() const { return m_defaultMode; }


    /**
     * @return the markVector
     */
    inline const TMarkValue * getMarkVector() const { return m_markVector; }


    inline const DWORD getMarkVectorLen() const { return m_markVectorLen; }


    void setMode(const char * mode) throw (ComponentErrors::CDBAccessExImpl, ReceiversErrors::ModeErrorExImpl);


	/**
	 * @return the lower limit of the RF coming from the L band feed (MHz)
	 */
	inline double const * const  getLBandRFMin() const { return m_LBandRFMin; }


	/**
	 * @return the upper limit of the RF coming from the L band feed (MHz)
	 */
	inline double const * const  getLBandRFMax() const { return m_LBandRFMax; }


	/**
	 * @return the lower limit of the RF coming from the P band feed (MHz)
	 */
	inline double const * const  getPBandRFMin() const { return m_PBandRFMin; }


	/**
	 * @return the upper limit of the RF coming from the P band feed (MHz)
	 */
	inline double const * const  getPBandRFMax() const { return m_PBandRFMax; }


	/**
	 * @return the start frequency of the IF coming from the L band feed (MHz)
	 */
	inline double const * const  getLBandIFMin() const { return m_LBandIFMin; }


	/**
	 * @return the start frequency of the IF coming from the P band feed (MHz)
	 */
	inline double const * const  getPBandIFMin() const { return m_PBandIFMin; }


	/**
	 * @return the bandwidth of the IF coming from the L band (MHz)
	 */
	inline double const * const  getLBandIFBandwidth() const{ 
        return  m_LBandIFBandwidth;
    }


	/**
	 * @return the bandwidth of the IF coming from the P band (MHz)
	 */
	inline double const * const  getPBandIFBandwidth() const{ 
        return  m_PBandIFBandwidth;
    }


	/**
	 * @return the number of IF chains for each feed
	 */
	inline const  DWORD& getIFs() const { return m_IFs; }


	/**
	 * @return the L band filter ID
	 */
	inline const  DWORD& getLBandFilterID() const { return m_LBandFilterID; }


	/**
	 * @return the P band filter ID
	 */
	inline const  DWORD& getPBandFilterID() const { return m_PBandFilterID; }

	/**
	 * @return the pointer to the L band polarizations array; the first element corresponds 
     * to the first IF and so on.
	 */
	inline Receivers::TPolarization const * const getLBandPolarizations() const { return m_LBandPolarizations; }


	/**
	 * @return the pointer to the P band polarizations array; the first element corresponds 
     * to the first IF and so on.
	 */
	inline Receivers::TPolarization const * const getPBandPolarizations() const { return m_PBandPolarizations; }

	inline IRA::CString getLBandPolarization() const { return m_LBandPolarization; }

	inline IRA::CString getPBandPolarization() const { return m_PBandPolarization; }


	/**
	 * @return the number of feeds
	 */
	inline const  DWORD& getFeeds() const { return m_feeds; }


	/**
	 * @return the default frequency for the synthesizer  (MHz)
	 */
	inline double const * const  getDefaultLO()  const { return m_DefaultLO; }


	/**
	 * @return the value of the fixed synthesizer used for the second conversion (MHz)
	 */
	inline double const * const  getFixedLO2() const { return m_FixedLO2; }


	/**
	 * @return lower limit for the synthesizer tuning (MHz)
	 */
	inline double const * const  getLOMin() const { return  m_LOMin; }


	/**
	 * @return upper limit for  the  synthesizer tuning (MHz)
	 */
	inline double const * const  getLOMax() const { return  m_LOMax; }


	inline double const * const  getLowpassFilterMin() const { return  m_LowpassFilterMin; }


	inline double const * const  getLowpassFilterMax() const { return  m_LowpassFilterMax; }


	/**
     * This member function is used to configure component by reading the configuration parameter from the CDB.
	 * This must be the first call before using any other function of this class.
	 * @throw ComponentErrors::CDBAccess
	 * @param Services pointer to the container services object
	*/
	void init(maci::ContainerServices *Services)  throw (
            ComponentErrors::CDBAccessExImpl,
            ComponentErrors::MemoryAllocationExImpl,
            ReceiversErrors::ModeErrorExImpl
    );
	

    std::vector<double> getLBandRFMaxFromMode(IRA::CString cmdMode) throw (
            ComponentErrors::CDBAccessExImpl,
            ComponentErrors::MemoryAllocationExImpl, 
            ReceiversErrors::ModeErrorExImpl);
 

    std::vector<double> getLBandRFMinFromMode(IRA::CString cmdMode) throw (
            ComponentErrors::CDBAccessExImpl,
            ComponentErrors::MemoryAllocationExImpl, 
            ReceiversErrors::ModeErrorExImpl);

    
private:
	IRA::CString m_dewarIPAddress;
	WORD m_dewarPort;
	IRA::CString m_LNAIPAddress;
	WORD m_LNAPort;
	IRA::CString m_switchIPAddress;
	WORD m_switchPort;
	DDWORD m_watchDogResponseTime;
	DDWORD m_watchDogSleepTime;
	DDWORD m_LNASamplingTime;
	DDWORD m_repetitionCacheTime;
	DDWORD m_repetitionExpireTime;
	IRA::CString m_localOscillatorInstance;
    maci::ContainerServices* m_services;
	IRA::CString m_mode;
	IRA::CString m_defaultMode;
	double *m_LBandRFMin;
	double *m_LBandRFMax;
	double *m_PBandRFMin;
	double *m_PBandRFMax;
	double *m_LBandIFMin;
	double *m_PBandIFMin;
	double *m_LBandIFBandwidth;
	double *m_PBandIFBandwidth;
    double *m_DefaultLO;
    double *m_FixedLO2;
    double *m_LOMin;
    double *m_LOMax;
    double *m_LowpassFilterMin;
    double *m_LowpassFilterMax;
	DWORD m_IFs;
	DWORD m_LBandFilterID;
	DWORD m_PBandFilterID;
	Receivers::TPolarization *m_LBandPolarizations;
	Receivers::TPolarization *m_PBandPolarizations;
    IRA::CString m_LBandPolarization;
    IRA::CString m_PBandPolarization;
	DWORD m_feeds;
	double *m_LBandLO;
	double *m_PBandLO;

	IRA::CDBTable *m_markTable;
	IRA::CDBTable *m_loTable;
	IRA::CDBTable *m_taperTable;
	IRA::CDBTable *m_feedsTable;

	TMarkValue *m_markVector;
	TLOValue * m_loVector;
	DWORD m_loVectorLen;
	DWORD m_markVectorLen;
	TTaperValue * m_taperVector;
	DWORD m_taperVectorLen;
	TFeedValue * m_feedVector; // length given by m_feeds
};


#endif
