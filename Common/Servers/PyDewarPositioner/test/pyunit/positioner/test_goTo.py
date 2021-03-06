import random
import time
import unittest
from maciErrType import CannotGetComponentEx
from DewarPositioner.positioner import Positioner, NotAllowedError
from DewarPositioner.cdbconf import CDBConf
from Acspy.Clients.SimpleClient import PySimpleClient


class GoTo(unittest.TestCase):

    def setUp(self):
        self.cdbconf = CDBConf()
        self.p = Positioner(self.cdbconf)
        try:
            client = PySimpleClient()
            self.device = client.getComponent('RECEIVERS/SRTKBandDerotator')
            self.using_mock = False
        except CannotGetComponentEx:
            print '\nINFO -> component not available: we will use a mock device'
            from DewarPositionerMockers.mock_components import MockDevice
            self.device = MockDevice()
            self.using_mock = True

    def tearDown(self):
        self.p.park()
        time.sleep(0.2)

    def test_goTo(self):
        """Verify the get method returns the position we set"""
        # Not allowed when the system is not yet configured
        self.assertRaises(NotAllowedError, self.p.goTo, 1.5)
        self.p.setup(siteInfo={}, source=None, device=self.device)
        self.cdbconf.setup('KKG') # Default configuration: FIXED
        self.p.goTo(1.5)
        time.sleep(2)
        self.assertAlmostEqual(self.p.getPosition(), 1.5, places=1)
        self.cdbconf.setConfiguration('BSC')


if __name__ == '__main__':
    unittest.main()
