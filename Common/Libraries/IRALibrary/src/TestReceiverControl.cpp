#include "ReceiverControl.h"
#include <cstdlib>

template<class T>
std::string any2string(T i) {
    std::ostringstream buffer;
    buffer << i;
    return buffer.str();
}


int main(int argc, char *argv[])
{
    if(argc != 5) {
        cerr << "\nWrong number of parameters. ";
        cerr << "Usage:\n\tTestReceiverControl dewar_IP dewar_port LNA_IP LNA_port\n" << endl;
        return 1;
    }

    std::string dewar_IP = std::string(argv[1]);
    unsigned int dewar_port = atoi(argv[2]);
    std::string lna_IP = std::string(argv[3]);
    unsigned int lna_port = atoi(argv[4]);
    std::vector<BYTE> data;


    // Test the constructor 
    try {
        cout << "\nConnecting to the boards (LNA and deward)..." << endl;
        ReceiverControl rc = ReceiverControl(dewar_IP, dewar_port, lna_IP, lna_port);
        cout << "Connection" << " done!" << endl << endl;
    
        // Test the setCalibrationOn()
        cout << "Test setCalibrationOn() with a reliable communication" << endl;
        rc.setCalibrationOn();
        cout << "Done!\n" << endl;

        // Test the setReliableCommOff()
        cout << "Test setCalibrationOff()" << endl;
        rc.setReliableCommOff();
        cout << "Done!\n" << endl;
    
        // Test the setCalibrationOn()
        cout << "Test setCalibrationOn()" << endl;
        rc.setCalibrationOn();
        cout << "Done!\n" << endl;

        // Test the setReliableCommOn()
        cout << "Test setReliableCommOn()" << endl;
        rc.setReliableCommOn();
        cout << "Done!\n" << endl;

        // Test the setCalibrationOff()
        cout << "Test setCalibrationOff() with a reliable communication" << endl;
        rc.setCalibrationOff();
        cout << "Done!\n" << endl;

        // Test the vacuum()
        cout << "Test vacuum() with a reliable communication" << endl;
        cout << "Vacuum value: " << rc.vacuum() << endl;
        cout << "Done!\n" << endl;

        // Test lna()
        cout << "Test lna() with a reliable communication" << endl;
        FetValues values = rc.lna(1, LEFT, 4);
        cout << "LNA values: (VD=" << values.VD << ", ID=" << values.ID;
        cout << ", VG=" << values.VG << ")" << endl;
        cout << "Done!\n" << endl;

        // Test the isCalibrationOn()
        cout << "Test isCalibrationOn() with a reliable communication" << endl;
        cout << "It the noise mark generator ON? " << (rc.isCalibrationOn() == true ? "yes" : "no") << endl;
        cout << "Done!\n" << endl;
    }
    catch(ReceiverControlEx& ex) {
        cout << ex.what() << endl;
        return 1;
    }


    return 0;
}