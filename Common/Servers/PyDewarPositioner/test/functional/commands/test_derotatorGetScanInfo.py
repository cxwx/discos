import unittest
import mocker
import time
from DewarPositioner.DewarPositionerImpl import DewarPositionerImpl


class DerotatorGetScanInfo(unittest.TestCase):

    def test_default(self):
        dp = DewarPositionerImpl()
        success, answer = dp.command('derotatorGetScanInfo')
        self.assertEqual(success, True)
        self.assertTrue('axis: MNG_NO_AXIS' in answer)

if __name__ == '__main__':
    unittest.main()
