import math
import time
import os
from datetime import datetime

import unittest2 # https://pypi.python.org/pypi/unittest2
import Management
import MinorServo
import Antenna
from Acspy.Common.TimeHelper import getTimeStamp
from Acspy.Clients.SimpleClient import PySimpleClient


__author__ = "Marco Buttu <mbuttu@oa-cagliari.inaf.it>"

class CheckScanTest(unittest2.TestCase):
    """Test the MinorServoBoss.checkScan() behavior"""

    def setUp(self):    
        code = 'KKG' if os.getenv('TARGETSYS') == 'SRT' else 'KKC'
        if hasattr(self, 'boss') and self.boss.isReady():
            pass
        else:
            client = PySimpleClient()
            self.boss = client.getComponent('MINORSERVO/Boss')
            if self.boss.getActualSetup() != code:
                self.boss.setup(code)
            t0 = datetime.now()
            while not self.boss.isReady() and (datetime.now() - t0).seconds < 60*5:
                time.sleep(2)
            self.assertTrue(self.boss.isReady()) # Timeout expired?

        self.assertEqual(self.boss.getActualSetup(), code)

        self.scan = MinorServo.MinorServoScan(
            range=50,
            total_time=500000000, # 50 seconds
            axis_code='SRP_TZ',
            is_empty_scan=False)

        self.antennaInfo = Antenna.TRunTimeParameters(
            targetName='dummy',
            azimuth=math.pi,
            elevation=math.pi/8,
            startEpoch=getTimeStamp().value + 100000000,
            onTheFly=False,
            slewingTime=100000000,
            section=Antenna.ANT_SOUTH,
            axis=Management.MNG_TRACK,
            timeToStop=0)

        self.msInfo = MinorServo.TRunTimeParameters( # TODO: remove me
            startEpoch=0,
            onTheFly=False,
            centerScan=0,
            scanAxis='SRP_TZ',
            timeToStop=0
        )

        self.boss.setElevationTracking('OFF')
        self.boss.setASConfiguration('OFF')
        axes, units = self.boss.getAxesInfo()
        self.idx = axes.index('SRP_TZ')

    def testStartASAP(self):
        """Check the behavior when the scan must start as soon as possible"""
        center = self.boss.getAxesPosition(getTimeStamp().value)[self.idx]
        startTime = 0
        res, msInfo = self.boss.checkScan(startTime, self.scan, self.antennaInfo)

        self.assertTrue(res)
        # Cannot know the startEpoch time...
        self.assertGreater(msInfo.startEpoch, getTimeStamp.now().value)
        self.assertTrue(msInfo.onTheFly)
        self.assertEqual(msInfo.centerScan, center)
        self.assertEqual(msInfo.scanAxis, 'SRP_TZ')
        # Cannot know the timeToStop time...
        self.assertGreater(msInfo.timeToStop, getTimeStamp.now().value+50000000)

    def _testStartAtGivenTime(self):
        """Check the behavior when the scan must start at a given time"""
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = getTimeStamp().value + 60*10**7 # Start in a minute from now
        msInfo = MinorServo.TRunTimeParameters
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)

        self.assertTrue(res)
        self.assertEqual(msInfo.startEpoch, startTime)
        self.assertTrue(msInfo.onTheFly)
        self.assertEqual(msInfo.centerScan, center)
        self.assertEqual(msInfo.scanAxis, 'SRP_TZ')
        # Cannot know the timeToStop time...
        self.assertGreater(msInfo.timeToStop, getTimeStamp.now().value+50000000)

    def _testTooFast(self):
        """The result must be False when the total time is not enough"""
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = getTimeStamp().value + 60*10**7 # Start in a minute from now
        msInfo = MinorServo.TRunTimeParameters
        self.scan.total_time = 500000 # 50 ms
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)
        self.assertFalse(res)

    def _testOutOfRange(self):
        """The result must be False in case of scan out of range"""
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = getTimeStamp().value + 60*10**7 # Start in a minute from now
        msInfo = MinorServo.TRunTimeParameters
        self.scan.range = 5000 # 5 meters...
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)
        self.assertFalse(res) # does it currently throw an exception?

    def _testTooCloseToNow(self):
        """The result must be False when the starting time is too close to now"""
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = getTimeStamp().value + 2 # Start in 2 seconds from now
        msInfo = MinorServo.TRunTimeParameters
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)
        self.assertFalse(res) # does it currently throw an exception?

    def _testStartASAPEmptyScan(self):
        """An empty scan when that must start as soon as possible"""
        self.scan = MinorServo.MinorServoScan(
            range=0,
            total_time=0,
            axis_code='',
            is_empty_scan=True)
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = 0
        msInfo = MinorServo.TRunTimeParameters
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)

        self.assertTrue(res)
        # Cannot know the startEpoch time...
        self.assertGreater(msInfo.startEpoch, getTimeStamp.now().value)
        self.assertFalse(msInfo.onTheFly)
        self.assertEqual(msInfo.centerScan, center)
        self.assertEqual(msInfo.scanAxis, '')
        self.assertEqual(msInfo.timeToStop, 0)

    def _testStartAtGivenTimeEmptyScan(self):
        """An empty scan that must start at a given time"""
        self.scan = MinorServo.MinorServoScan(
            range=0,
            total_time=0,
            axis_code='',
            is_empty_scan=True)
        center = self.boss.getAxesPosition(getTimeStamp().value)['SRP_TZ']
        startTime = getTimeStamp().value + 60*10**7 # Start in a minute from now
        msInfo = MinorServo.TRunTimeParameters
        res = self.boss.checkScan(startTime, self.scan, self.antennaInfo, msInfo)

        self.assertTrue(res)
        self.assertEqual(msInfo.startEpoch, startTime)
        self.assertFalse(msInfo.onTheFly)
        self.assertEqual(msInfo.centerScan, center)
        self.assertEqual(msInfo.scanAxis, '')
        self.assertGreater(msInfo.timeToStop, 0)


if __name__ == '__main__':
    unittest2.main()


