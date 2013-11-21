__all__ = ['antenna', 'recorder', 'ato', 'ANT_J2000', 'ACU_NEUTRAL']

import datetime
import logging
import time
import sys
from os.path import join, exists, curdir, abspath
from os import mkdir
from Acspy.Common.TimeHelper import getTimeStamp


class Recorder(object):
    def __init__(self, backend_name='BACKENDS/XBackends', writer_name='MANAGEMENT/FitsZilla',\
            projectName='ison', observerName='Gavino', file_name=''):
        wclient = PySimpleClient()
        bclient = PySimpleClient()
        self.writer = wclient.getComponent(writer_name)
        self.backend = bclient.getComponent(backend_name)
        # Data stuffs
        self.projectNickname = '_'.join(projectName.split())
        datadir = '%s_data/' %self.projectNickname
        if not exists(datadir):
            mkdir(datadir)
        full_path = join(abspath(curdir), datadir)
        self.scanSetup = Management.TScanSetup(
            scanTag=0, 
            scanId=1, 
            projectName=projectName,
            observerName=observerName,
            path=full_path, 
            extraPath='',
            baseName='',
            scanLayout='', 
            schedule='noSchedule', 
            device=0)
        self.subScanSetup = Management.TSubScanSetup(
            startUt = getTimeStamp().value,
            subScanId=0,
            axis=Management.MNG_NO_AXIS,
            targetID='',
            extraPath='',
            baseName='')
        _setStartTime(datetime.datetime.utcnow())
        _setStopTime(datetime.datetime.utcnow())

    def start(self):
        self.writer.reset();
        time.sleep(1)
        self.backend.connect(self.writer);
        self._updateScanSetup()
        self.writer.startScan(self.scanSetup)
        time.sleep(1)
        self.backend.sendHeader()
        self._updateSubScanSetup()
        self.writer.startSubScan(self.subScanSetup);
        self.backend.sendData(self.subScanSetup.startUt)
        _setStartTime(self.subScanSetup)
        while getTimeStamp().value < self.subScanSetup.startUt + 110000000: 
            time.sleep(1)

    def _updateSubScanSetup(self):
        self.subScanSetup.startUt = getTimeStamp().value + 10000000 # Add one second
        self.subScanSetup.subScanId += 1 

    def stop(self):
        self.backend.sendStop();
        time.sleep(2)
        while self.writer.isRecording():
            time.sleep(1)
        _setStopTime(datetime.datetime.utcnow())
        self.writer.stopScan() 
        time.sleep(1)
        self.backend.terminate()
        self.backend.disconnect()
        self.writer.closeReceiver()

    def _updateScanSetup(self):
        self.scanSetup.scanTag += 1
        utc = datetime.datetime.utcnow()
        self.scanSetup.extraPath = '%04d-%s' %(self.scanSetup.scanTag, utc.strftime("%Y%m%d-%H%M%S"))

    def _setStartTime(self, value):
        self.startTime = value

    def _setStartTime(self, value):
        self.stopTime = value

try:
    from Acspy.Clients.SimpleClient import PySimpleClient
    from maciErrTypeImpl import CannotGetComponentExImpl
    from Antenna import ANT_J2000, ACU_NEUTRAL
    from AntennaErrors import AntennaErrorsEx
    import Management
    client = PySimpleClient()
    antenna = client.getComponent('ANTENNA/Boss')
    recorder = Recorder()
    # If the backend rounds the time to a next upper value, set a negative `ato` 
    ato = -2 # Offset to add to the backend acquisition time
except Exception, e:
    logging.exception('Cannot get the ACS modules and components')
    print('ERROR: Cannot get the ACS modules and components')
    raise


