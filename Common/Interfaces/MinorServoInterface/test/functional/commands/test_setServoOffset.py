from __future__ import with_statement

import math
import time
import unittest2

import MinorServo
import Management
import Antenna

from Acspy.Clients.SimpleClient import PySimpleClient
from MinorServoErrors import MinorServoErrorsEx
from Acspy.Common.TimeHelper import getTimeStamp

__author__ = "Marco Buttu <mbuttu@oa-cagliari.inaf.it>"


class TestSetServoOffsetCmd(unittest2.TestCase):
    """Test the setServoOffset command"""

    def setUp(self):
        client = PySimpleClient()
        self.boss = client.getComponent('MINORSERVO/Boss')

    def test_wrong_axis_code(self):
        success, answer = self.boss.command('setServoOffset=FOO_TX,0')
        self.assertFalse(success)


if __name__ == '__main__':
    unittest2.main()