#!/usr/bin/env python
# Author: Marco Buttu <m.buttu@oa-cagliari.inaf.it>
# Copyright: This module has been placed in the public domain.
"""Use this module to test the minor servo boss startScan() method"""

import unittest
import threading
import os
import time
from random import randrange as rr
from socket import *       
from Acspy.Clients.SimpleClient import PySimpleClient
from Acspy.Common import TimeHelper
from maciErrTypeImpl import CannotGetComponentExImpl
import datetime

class NackEx(Exception):
    pass

class TestScan(unittest.TestCase):
    """Test the minor servo properties."""

    def setUp(self):
        self.pyclient = PySimpleClient()
        self.cmd_num = 0
        self.srp = self.pyclient.getComponent("MINORSERVO/SRP")
        self.boss = self.pyclient.getComponent("MINORSERVO/Boss")

    def tearDownClass():
        print "in tearDown"
        self.sockobj.close()
        self.pyclient.releaseComponent(self.srp)
        self.pyclient.releaseComponent(self.boss)

    def test_startScan(self):
        """Test the startScan() method."""
        setup_code = 'CCB'
        self.boss.setup(setup_code)
        time.sleep(1) # Wait a bit, until the boss begins the configuration process
        
        counter = 0
        print "\nExecuting the minor servo %s setup. Wait a bit ..." %setup_code
        delay_ready = 2
        while not self.boss.isReady(): # Wait until the minor servo boss is ready
            time.sleep(delay_ready) # Wait a bit, until the boss is active
            counter += delay_ready
            if counter > 240:
                self.assertTrue(counter > 240)
                return

        print "\nThe MinorServoBoss is ready."

        time.sleep(10)
        delay = 15 * 10 ** 7 # 15 seconds
        starting_time = TimeHelper.getTimeStamp().value + delay # Start the scan in `delay` seconds from now
        range_ = 8.0 # mm 
        total_time = 4 * 10 ** 7 # 4 seconds
        mm_per_sec = range_ / total_time
        actPos_obj = self.srp._get_actPos()
        actPos, cmp = actPos_obj.get_sync()
        print "\nCalling the startScan() method, starting in %d s from now." %(delay/10**7)
        self.boss.startScan(starting_time, range_, total_time, 2, "SRP")
        expected_positions = {}
        sampling_time = 2000000 # 200ms
        exe_times = range(starting_time, starting_time + total_time + sampling_time, sampling_time)
        expected_position = actPos[:]
        for i, t in enumerate(exe_times):
            expected_position[2] = actPos[2] -range_/2 + (t - starting_time) * mm_per_sec
            expected_positions[t] = expected_position[:]
        for i in range(1, 5):
            expected_positions[t + sampling_time * i] = expected_position[:]

        actual_time = TimeHelper.getTimeStamp().value 
        print "\nWaiting for the scan..."
        while actual_time < starting_time:
            actual_time = TimeHelper.getTimeStamp().value 
            time.sleep(0.2)

        print "\nScanning..."
        while actual_time < starting_time + total_time:
            actual_time = TimeHelper.getTimeStamp().value 
            time.sleep(0.2)
        print "\nDone!"

        print "Writing the file..."
        out_file = open('../data/scan.data', 'w')
        out_file.write('# Running time: %s\n\n' %str(datetime.datetime.now()))
        out_file.write('# Starting time: %d\n' %starting_time)
        out_file.write('# Total time: %d\n\n s' %(total_time/10**7))
        out_file.write('# Range: %d mm\n\n' %range_)
        out_file.write('# History sampling time: %d\n' %sampling_time)

        hpositions = {}
        for t in sorted(expected_positions.keys()):
            print "time: ", t
            ep = ' '.join(['%.4f' %item for item in expected_positions[t]])
            hp = self.srp.getPositionFromHistory(t)
            out_file.write("\n#" + "-"*68)
            out_file.write("\nTime: starting time + %0*d ms" %(5, (t - starting_time) / (10**4)))
            out_file.write("\nExpected: %s" %ep)
            out_file.write("\nHistory:  %s" %' '.join(['%.4f' %item for item in hp]))
            hpositions[t] = hp
            print "-"*40
            print "Timestamp: ", t
            print "History : ", hp
            print "Expected: ", ep
            print
        pos = [expected_positions[t] for t in sorted(expected_positions.keys())]
        hpos = [hpositions[t] for t in sorted(hpositions.keys())]
        max_diff = max([abs(p[2] - h[2]) for p, h in zip(pos, hpos)])
        out_file.write("\n\n# Max diff: %.4f" %max_diff)
        print "Max diff: ", max_diff
        print "Done!"

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestScan)
    unittest.TextTestRunner(verbosity=2).run(suite)
    print "\n" + "="*70

