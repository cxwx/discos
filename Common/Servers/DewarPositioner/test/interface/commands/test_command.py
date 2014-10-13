import unittest2
import mocker
from DewarPositioner.DewarPositionerImpl import DewarPositionerImpl


class CommandTest(unittest2.TestCase):
    """Test the DewarPositioner.command() method"""

    def setUp(self):
        self.dp = DewarPositionerImpl()

    def test_unknown_command(self):
        """Verify the answer is the exception message and success is False"""
        cmd = 'pippo'
        success, answer = self.dp.command(cmd)
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'command %s does not exist' %cmd)
        cmd = 'pippo=='
        success, answer = self.dp.command(cmd)
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'invalid command: maybe there are')

    def test_wrong_parameters(self):
        """Verify the answer is the exception message and success is False"""
        success, answer = self.dp.command('derotatorSetRewindingMode=a,b')
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'mode A unknown')
        success, answer = self.dp.command('derotatorSetRewindingMode=')
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'mode  unknown')
        success, answer = self.dp.command('derotatorPark=a,b')
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'positioner not ready: a setup')

    def test_command_fails(self):
        """Verify the answer is the exception message and success is False"""
        success, answer = self.dp.command('derotatorSetConfiguration=FIXED')
        self.assertEqual(success, False)
        self.assertRegexpMatches(answer, 'cannot set the configuration')

    def test_command_returns_None(self):
        """Verify the answer is an empty string and success is True"""
        success, answer = self.dp.command('derotatorSetRewindingMode=AUTO')
        self.assertEqual(success, True)
        self.assertRegexpMatches(answer, '')


if __name__ == '__main__':
    unittest2.main()
