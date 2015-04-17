from __future__ import with_statement
import random 
import math
import time
import os
from datetime import datetime

import unittest2 # https://pypi.python.org/pypi/unittest2
import Management
import MinorServo
import Antenna

from MinorServoErrors import MinorServoErrorsEx
from Acspy.Common.TimeHelper import getTimeStamp
from Acspy.Clients.SimpleClient import PySimpleClient
from Acspy.Util import ACSCorba

__author__ = "Marco Buttu <mbuttu@oa-cagliari.inaf.it>"

class PositionTest(unittest2.TestCase):

    telescope = os.getenv('TARGETSYS')

    @classmethod
    def setUpClass(cls):
        cls.client = PySimpleClient()
        cls.boss = cls.client.getComponent('MINORSERVO/Boss')
        
    @classmethod
    def tearDownClass(cls):
        cls.client.releaseComponent('MINORSERVO/Boss')

    def setUp(self):
        self.axis_code='SRP_TZ' if self.telescope == 'SRT' else 'Z'
        setupCode = 'KKG' if self.telescope == 'SRT' else 'CCC'
        # Wait (maximum one minute) in case the boss is parking
        if self.boss.isParking():
            t0 = datetime.now()
            while self.boss.isParking() and (datetime.now() - t0).seconds < 60:
                time.sleep(2)
            if self.boss.isParking():
                self.fail('The system can not exit form a parking state')

        if self.boss.getActualSetup() != setupCode:
            self.boss.setup(setupCode)
            # Wait (maximum 5 minutes) in case the boss is starting
            t0 = datetime.now()
            while not self.boss.isReady() and (datetime.now() - t0).seconds < 60*5:
                time.sleep(2)
            if not self.boss.isReady():
                self.fail('The system is not ready for executing the tests')
        self.boss.setElevationTracking('OFF')
        self.boss.setASConfiguration('OFF')
        axes, units = self.boss.getAxesInfo()
        self.idx = axes.index(self.axis_code)

    def tearDown(self):
        # self.boss.clearUserOffset(self.axis_code)
        self.boss.setUserOffset(self.axis_code, 0)
        self.wait_tracking()

    def test_get_current_position(self):
        timestamp = getTimeStamp().value
        position = self.get_position()
        position_now = self.get_position(timestamp)
        self.assertAlmostEqual(position, position_now, delta=0.1)

    def test_get_offset_position(self):
        position = self.get_position()
        self.boss.setUserOffset(self.axis_code, 10)
        self.wait_tracking()
        position_now = self.get_position()
        self.assertAlmostEqual(position_now, position + 10, delta=0.1)

    def test_get_past_position(self):
        timestamp = getTimeStamp().value
        position = self.get_position()
        self.boss.setUserOffset(self.axis_code, 10)
        self.wait_tracking()
        position_past = self.get_position(timestamp)
        self.assertAlmostEqual(position, position_past, delta=0.1)

    def test_get_past_position_with_sleep(self):
        timestamp = getTimeStamp().value
        position = self.get_position()
        self.boss.setUserOffset(self.axis_code, delta=10)
        self.wait_tracking()
        time.sleep(10)
        position_past = self.get_position(timestamp)
        self.assertAlmostEqual(position, position_past, delta=0.1)

    def wait_tracking(self):
        while not self.boss.isTracking():
            time.sleep(0.1)

    def get_position(self, timestamp=0):
        return self.boss.getAxesPosition(timestamp)[self.idx]
        

if __name__ == '__main__':
    if 'Configuration' in os.getenv('ACS_CDB'):
        unittest2.main(verbosity=2, failfast=True) # Real test using the antenna CDB
    else:
        from PyMinorServoTest import simunittest
        simunittest.run(PositionTest)