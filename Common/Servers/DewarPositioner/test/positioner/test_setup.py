import unittest2
import random
import time
from maciErrType import CannotGetComponentEx
from DewarPositioner.positioner import Positioner, PositionerError
from Acspy.Clients.SimpleClient import PySimpleClient


class PositionerSetupTest(unittest2.TestCase):

    def setUp(self):
        self.cdb_info = {
                'updating_time': 0.1,
                'rewinding_timeout': 1.5,
                'rewinding_sleep_time': 1
        }
        try:
            client = PySimpleClient()
            self.device = client.getComponent('RECEIVERS/SRTKBandDerotator')
        except CannotGetComponentEx:
            print '\nINFO -> component not available: we will use a mock device'
            from DewarPositionerTest.mock_components import MockDevice
            self. device = MockDevice()

    def test_starting_pos(self):
        """Verify it sets properly the starting position"""
        p = Positioner(self.cdb_info)
        p.setup(site_info={}, source=None, device=self.device, starting_position=2)
        self.assertEqual(p.getStartingPosition(), 2)

    def test_device_name(self):
        """Verify it sets properly the device_name"""
        p = Positioner(self.cdb_info)
        p.setup(site_info={}, source=None, device=self.device)
        self.assertEqual(p.getDeviceName(), self.device._get_name())


if __name__ == '__main__':
    unittest2.main()
