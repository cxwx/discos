from __future__ import with_statement

import os
import time
import datetime
import subprocess

import unittest

import Management
import MinorServo

from Acspy.Clients.SimpleClient import PySimpleClient
from acswrapper.system import acs
from acswrapper.containers import (
    Container, ContainerError, start_containers_and_wait,
    stop_containers_and_wait
)

__author__ = "Marco Buttu <mbuttu@oa-cagliari.inaf.it>"


FNULL = open(os.devnull, 'w')


class TestFailure(unittest.TestCase):

    telescope = os.getenv('STATION')

    @classmethod
    def setUpClass(cls):
        if not acs.is_running():
            acs.start()
        cls.containers = [
            Container('MinorServoContainer', 'cpp'),
            Container('MinorServoBossContainer', 'cpp'),
        ]
        try:
            start_containers_and_wait(cls.containers)
        except ContainerError, ex:
            cls.fail(ex.message)

        cls.client = PySimpleClient()
        cls.boss = cls.client.getComponent('MINORSERVO/Boss')
        cls.srp = cls.client.getComponent('MINORSERVO/SRP')
        
    @classmethod
    def tearDownClass(cls):
        cls.client.releaseComponent('MINORSERVO/Boss')
        cls.client.releaseComponent('MINORSERVO/SRP')
        stop_containers_and_wait(cls.containers)

    def setUp(self):
        self.setup_code = "CCB" if self.telescope == "SRT" else "CCC"
        self.boss.setup(self.setup_code)

    def tearDown(self):
        if self.srp.isReady():
            self.boss.park()
            time.sleep(0.2)
            self.wait_parked()

    def test_failure_during_tracking(self):
        """Make sure a failure during tracking gets notified"""
        counter = 0 # Seconds
        now = time_ref = datetime.datetime.now()
        while not self.boss.isReady() or (time_ref - now).seconds < 20:
            time.sleep(1)
            now = datetime.datetime.now()
        time.sleep(4)
        self.srp.disable(0)
        time.sleep(5)
        status_obj = self.boss._get_status()
        status_value, comp = status_obj.get_sync()
        self.assertEqual(status_value, Management.MNG_FAILURE)

    def test_failure_during_positioning(self):
        """Make sure a failure during SRP positioning gets notified"""
        counter = 0 # Seconds
        now = time_ref = datetime.datetime.now()
        while not self.boss.isReady() or (time_ref - now).seconds < 20:
            time.sleep(1)
            now = datetime.datetime.now()
        time.sleep(4)
        try:
            subprocess.Popen(['srt-mscu-sim', 'setpos_NAK'], stdout=FNULL, stderr=FNULL)
            time.sleep(5)
            status_obj = self.boss._get_status()
            status_value, comp = status_obj.get_sync()
            self.assertEqual(status_value, Management.MNG_FAILURE)
        finally:
            subprocess.Popen(['srt-mscu-sim', 'setpos_ACK'], stdout=FNULL, stderr=FNULL)
            time.sleep(2)

    def wait_parked(self):
        while self.boss.isParking():
            time.sleep(0.1)


if __name__ == '__main__':
    if 'Configuration' in os.getenv('ACS_CDB'):
        unittest.main() # Real test using the antenna CDB
    else:
        from testing import simulator
        simulator.run(TestFailure, 'srt-mscu-sim')
