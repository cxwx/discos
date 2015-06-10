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

class SetupTest(unittest2.TestCase):

    telescope = os.getenv('TARGETSYS')
    
    def setUp(self):
        self.client = PySimpleClient()
        self.boss = self.client.getComponent('MINORSERVO/Boss')
        
    def tearDown(self):
        self.client.releaseComponent('MINORSERVO/Boss')
        self.client.disconnect()

    def test_elevation_tracking_ON(self):
        """The setup turns the elevation tracking on"""
        self.boss.setup('KKG')
        self.wait_until_ready()
        self.assertTrue(self.boss.isElevationTrackingEn())

    def test_as_configuration_ON(self):
        """The setup turns the AS configuration on"""
        self.boss.setup('KKG')
        self.wait_until_ready()
        self.assertEqual(self.boss.getActualSetup(), 'KKG_ASACTIVE')

    def wait_until_ready(self, timeout=20):
        t0 = datetime.now()
        while not self.boss.isReady():
            if (datetime.now() - t0).seconds > timeout:
                assert False, 'The minor servo boss is not ready'
            else:
                time.sleep(1)


if __name__ == '__main__':
    if 'Configuration' in os.getenv('ACS_CDB'):
        unittest2.main(verbosity=2, failfast=True) # Real test using the antenna CDB
    else:
        from PyMinorServoTest import simunittest
        simunittest.run(SetupTest)
