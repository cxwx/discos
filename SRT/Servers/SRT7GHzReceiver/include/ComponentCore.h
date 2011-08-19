#ifndef _COMPONENTCORE_H_
#define _COMPONENTCORE_H_

/* **************************************************************************************************** */
/* IRA Istituto di Radioastronomia                                                                      */
/*                                                                                                      */
/* This code is under GNU General Public License (GPL).                                                 */
/*                                                                                                      */
/* Who                                when            What                                              */
/* Andrea Orlati(aorlati@ira.inaf.it) 03/08/2011     Creation                                         */

#include "Configuration.h"
#include <ReceiverControl.h>
#include <LocalOscillatorInterfaceC.h>
#include <ReceiversErrors.h>

/**
 * This class contains the code of almost all the features  of the component
 * @author <a href=mailto:a.orlati@ira.cnr.it>Andrea Orlati</a>,
 * Istituto di Radioastronomia, Italia
 * <br>
  */
class CComponentCore {
public:
	/**
	 * Constructor
	 */
	CComponentCore();

	/**
	 * Destructor
	 */
	virtual ~CComponentCore();

	/**
	 * This method initializes the object
	 * @param service pointer to container services object provided by the container
	 */
	virtual void initialize(maci::ContainerServices* services);

	/**
	 * This method prepares the object for execution.
	 * @return the pointer to the configuration class
	 */
	virtual CConfiguration const * const execute() throw (ComponentErrors::CDBAccessExImpl,ComponentErrors::MemoryAllocationExImpl,ComponentErrors::SocketErrorExImpl);

	/**
	 * This function is responsible to free all allocated resources
	 */
	virtual void cleanup();

	/*
	 * It sets the local oscillator. Only the first value is considered in this case, since the receiver has just one synthesizer. Before commanding the new value some check are done. The the correspnding signal
	 * amplitude is computed.
	 * @param lo lists of values for the local oscillator (MHz), one for each IF. In that case just the first one is significant. In a -1 is passed the present value is kept
	 * @throw  ComponentErrors::ValidationErrorExImpl
	 * @throw ComponentErrors::ValueOutofRangeExImpl
	 * @throw ComponentErrors::CouldntGetComponentExImpl
	 * @throw ComponentErrors::CORBAProblemExImpl
	 * @thorw ReceiversErrors::LocalOscillatorErrorExImpl
	 */
	void setLO(const ACS::doubleSeq& lo) throw (ComponentErrors::ValidationErrorExImpl,ComponentErrors::ValueOutofRangeExImpl,ComponentErrors::CouldntGetComponentExImpl,
			ComponentErrors::CORBAProblemExImpl,ReceiversErrors::LocalOscillatorErrorExImpl);

	/**
	 * It allows to change the operating mode of the receiver. If the mode does not correspond to a valid mode an error is thrown.
	 * @param  mode mode code as a string
	 */
	void setMode(const char * mode) throw  (ReceiversErrors::ModeErrorExImpl,ComponentErrors::ValidationErrorExImpl,ComponentErrors::ValueOutofRangeExImpl,
			ComponentErrors::CouldntGetComponentExImpl,ComponentErrors::CORBAProblemExImpl,ReceiversErrors::LocalOscillatorErrorExImpl);

	/**
	 * It activate the receiver, in other words it allows to setup the default configuration, and to check is the LNA are turned on.
	 */
	void activate() throw (ReceiversErrors::ModeErrorExImpl,ComponentErrors::ValidationErrorExImpl,ComponentErrors::ValueOutofRangeExImpl,
			ComponentErrors::CouldntGetComponentExImpl,ComponentErrors::CORBAProblemExImpl,ReceiversErrors::LocalOscillatorErrorExImpl);

	/**
	 * It allows to compute the value of the calibration mark for any given sub bands in the IF space.
	 * @param result this the sequence of computed mark values, first entry correspond to first sub band and so on....
	 * @param freqs  list of start frequencies (MHz)
	 * @param bandwidth list of the band widths (MHz)
	 * @param feeds list of feed identifier, it allows to specifies form which feed the sub band comes from. In that case it is neglected since the receiver is a single feed
	 * @param ifs list of IF identifier, it allows to specifies from which receiver IF the sub band comes from.
	 */
	void getCalibrationMark(ACS::doubleSeq& result,const ACS::doubleSeq& freqs,const ACS::doubleSeq& bandwidths,const ACS::longSeq& feeds,
			const ACS::longSeq& ifs) throw (ComponentErrors::ValidationErrorExImpl,ComponentErrors::ValueOutofRangeExImpl);

	/**
	 * It computes the taper given a reference band.
	 * @param freq start frequency of the reference band
	 * @param bw width of the reference band
	 * @param feed feed number
	 * @param ifNumber IF chain identifier
	 * @param waveLen wave length of the reference band, the band is transformed in a real sky observed band and the the central frequency is taken
	 */
	double getTaper(const double& freq,const double& bw,const long& feed,const long& ifNumber,double& waveLen) throw (
			ComponentErrors::ValidationErrorExImpl,ComponentErrors::ValueOutofRangeExImpl);

	/**
	 * It turns the calibration diode on.
	 */
	void calOn() throw (ComponentErrors::ValidationErrorExImpl,ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It turns the calibration diode off
	 */
	void calOff() throw (ComponentErrors::ValidationErrorExImpl,ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the control board the current value of the vacuum
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateVacuum() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It check if the vacuum pump is on
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateVacuumPump() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It checks if the vacuum valve is opened or not
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateVacuumValve() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the control board the current cryo temperature measured near the cool head
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateCryoCoolHead() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the control board the current cryo temperature measured near the cool head window
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateCryoCoolHeadWin() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the control board the current cryo temperature measured near the LNA
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateCryoLNA() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the LNA control board the current values of current and voltage of gate and drain of the transistors
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateLNAControls() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It reads and updates from the control board the current cryo temperature measured near the LNA window
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateCryoLNAWin() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It checks if the Dewar power box is in remote or not
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateIsRemote() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * It checks if the cool head is turned on or not
	 * @throw ReceiversErrors::ReceiverControlBoardErrorExImpl
	 */
	void updateCoolHead() throw (ReceiversErrors::ReceiverControlBoardErrorExImpl);

	/**
	 * This is getter method. No need to make it thread safe......
	 * @return the current value of the vacuum in mbar
	 */
	double  getVacuum() const { return m_vacuum; }

	/**
	 * This is getter method. No need to make it thread safe......
	 * @return the current value of the cryogenic temperature at cool head in °K
	 */
	double getCryoCoolHead() const { return m_cryoCoolHead; }

	/**
	 * This is getter method. No need to make it thread safe......
	  * @return the current value of the cryogenic temperature at cool head window in °K
	 */
	double getCryoCoolHeadWin() const { return m_cryoCoolHeadWin; }

	/**
	 * This is getter method. No need to make it thread safe......
	  * @return the current value of the cryogenic temperature at LNA in °K
	 */
	double getCryoLNA() const { return m_cryoLNA; }

	/**
	 * This is getter method. No need to make it thread safe......
	 * @return the current value of the cryogenic temperature at LNA  window in °K
	 */
	double getCryoLNAWin() const { return m_cryoLNAWin; }

	/**
	 * This is getter method. No need to make it thread safe......
	 * @return the current status word
	 */
	DWORD getStatusWord() const  { return  m_statusWord; }

	/**
	 * This is getter method. In this case, since it makes use of some class members that could be changed by other methods it is advisable to protect this method with the class mutex.
	 * @param control name of the parameter that must be returned
	 * @param ifs Intermediate frequency identifier, it permits to select which amplification chain we are interested in
	 * @return a specific value of from the transistor control parameters
	 */
	double getFetValue(const IRA::FetValue& control,const DWORD& ifs);

	/**
	 * It returns the feed geometry of the receiver with respect to the central one. For this implementation it is just a placeholder since there is just one feed.
	 */
	long getFeeds(ACS::doubleSeq& X,ACS::doubleSeq& Y,ACS::doubleSeq& power);

	/**
	 * It returns back the current local oscillator frequency settings.
	 * @param lo output sequence
	 */
	void getLO(ACS::doubleSeq& lo);

	/**
	 * It returns back the current bandwidth for each IF.
	 * @param bw output sequence
	 */
	void getBandwidth(ACS::doubleSeq& bw);

	/**
	 * It returns back the current start frequency for each IF.
	 * @param sf output sequence
	 */
	void getStartFrequency(ACS::doubleSeq& sf);

	/**
	 * It returns back the current polarization for each IF.
	 * @param pol output sequence
	 */
	void getPolarization(ACS::longSeq& pol);

	/**
	 * It returns the current operating mode of the receiver.
	 * @return output string
	 */
	const IRA::CString& getSetupMode();

	/**
	 * It returns the number of IF chains
	 * @return output value
	 */
	const DWORD& getIFs();

	/**
	 * It returns the number of feeds
	 * @return output value
	 */
	const DWORD& getFeeds();

protected:
	/**
	 * Obtain a valid reference to the local oscillator device
	 */
	void loadLocalOscillator()  throw (ComponentErrors::CouldntGetComponentExImpl);

	/**
	 * used to free the reference to the local oscillator device
	 */
	void unloadLocalOscillator();
private:
	enum TStatusBit {
		LOCAL=0,
		VACUUMSENSOR=1,
		VACUUMPUMPSTATUS=2,
		VACUUMPUMPFAULT=3, //**/
		VACUUMVALVEOPEN=4,
		COOLHEADON=5,
		NOISEMARK=6,
		NOISEMARKERROR=7, //**/
		EXTNOISEMARK=8, //**/
		CONNECTIONERROR=9, /**/
		UNLOCKED=10   //**/
	};

	CConfiguration m_configuration;
	maci::ContainerServices* m_services;
	BACIMutex m_mutex;
	IRA::ReceiverControl *m_control; // this object is thread safe
	Receivers::LocalOscillator_var m_localOscillatorDevice;
	bool m_localOscillatorFault;
	double m_localOscillatorValue;
	ACS::doubleSeq m_startFreq;
	ACS::doubleSeq m_bandwidth;
	ACS::longSeq m_polarization;
	IRA::CString m_setupMode;
	double m_vacuum;
	double m_cryoCoolHead;
	double m_cryoCoolHeadWin;
	double m_cryoLNA;
	double m_cryoLNAWin;
	IRA::FetValues m_fetValues;
	DWORD m_statusWord;

	/**
	 * This function will set the a status bit
	 */
	inline void setStatusBit(TStatusBit bit) { m_statusWord |= 1 << bit; }

	/**
	 * This function will unset (clear) a status bit
	 */
	inline void clearStatusBit(TStatusBit  bit) { m_statusWord &= ~(1 << bit); }

	/**
	 * This function check is a bit is set or not
	 */
	inline bool checkStatusBit(TStatusBit bit) { return m_statusWord & (1 << bit); }

	double linearFit(double *X,double *Y,const WORD& size,double x);
	/************************ CONVERSION FUNCTIONS **************************/
	// Convert the voltage value of the vacuum to mbar
	static double voltage2mbar(double voltage) { return(pow(10, 1.5 * voltage - 12)); }
	// Convert the voltage value of the temperatures to Kelvin
	static double voltage2Kelvin(double voltage) {
	    return voltage < 1.12 ?  (660.549422889947 * pow(voltage, 6)) - (2552.334255456860 * pow(voltage, 5)) + (3742.529989384570 * pow(voltage, 4))
	        - (2672.656926956470 * pow(voltage, 3)) + (947.905578508975 * pow(voltage, 2)) - 558.351002849576 * voltage + 519.607622398508  :
	        (865.747519105672 * pow(voltage, 6)) - (7271.931957100480 * pow(voltage, 5)) + (24930.666241800500 * pow(voltage, 4))
	        - (44623.988512320400 * pow(voltage, 3)) + (43962.922216886600 * pow(voltage, 2)) - 22642.245121997700 * voltage + 4808.631312836750;
	}
	// Convert the ID voltage value to the mA value
	static double currentConverter(double voltage) { return(10 * voltage); }
	// Convert the VD and VG voltage values using a right scale factor
	static double voltageConverter(double voltage) { return(voltage); }
};


#endif
