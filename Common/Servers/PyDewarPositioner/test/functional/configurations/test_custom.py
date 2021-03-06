import unittest
import time
from ComponentErrors import ValidationErrorEx, NotAllowedEx
from Acspy.Clients.SimpleClient import PySimpleClient


class CustomTest(unittest.TestCase):
    """Test the CUSTOM Configuration"""

    def setUp(self):
        client = PySimpleClient()
        self.dp = client.getComponent('RECEIVERS/DewarPositioner')
        self.dp.setup('KKG')
    
    def tearDown(self):
        self.dp.park()
        time.sleep(0.5)

    def test_setConfiguration(self):
        self.dp.setConfiguration('CUSTOM')
        self.assertEqual(self.dp.getConfiguration(), 'CUSTOM')

    def test_setPosition(self):
        """Verify we can set the position"""
        self.dp.setConfiguration('CUSTOM')
        self.dp.setPosition(1)


if __name__ == '__main__':
    unittest.main()
