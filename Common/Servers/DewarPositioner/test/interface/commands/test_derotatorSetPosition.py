import unittest2
import mocker
import time
from DewarPositioner.DewarPositionerImpl import DewarPositionerImpl


class DerotatorSetPositionTest(unittest2.TestCase):
    """Test the Set[Get]Position commands"""

    def test_setupcmd(self):
        dp = DewarPositionerImpl()
        dp.command('derotatorSetup=KKG') # Default conf: FIXED
        success, answer = dp.command('derotatorSetPosition=1.5d')
        self.assertEqual(success, True)
        time.sleep(2)
        success, answer = dp.command('derotatorGetPosition')
        self.assertEqual(success, True)
        value = float(answer.rstrip('d'))
        self.assertAlmostEqual(value, 1.5, places=1)
        success, answer = dp.command('derotatorSetPosition=0')
        self.assertEqual(success, True)
        time.sleep(2)
        success, answer = dp.command('derotatorGetPosition')
        self.assertEqual(success, True)
        value = float(answer.rstrip('d'))
        self.assertAlmostEqual(value, 0, places=1)
        success, answer = dp.command('derotatorSetPosition=foo')
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'invalid')

if __name__ == '__main__':
    unittest2.main()
