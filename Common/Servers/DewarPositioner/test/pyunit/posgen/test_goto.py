import unittest2
from DewarPositioner.posgenerator import PosGenerator, PosGeneratorError


class PosGeneratorGoToTest(unittest2.TestCase):

    def test_goto(self):
        """The goto take a position and yield it"""
        posgen = PosGenerator()
        position = 10
        for pos in posgen.goto(position):
            self.assertEqual(pos, position)

if __name__ == '__main__':
    unittest2.main()
