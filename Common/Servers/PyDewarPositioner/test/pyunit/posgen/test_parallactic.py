import unittest
import mocker
from math import sin, cos, tan, atan, radians, degrees, pi, atan2
from DewarPositioner.posgenerator import PosGenerator, PosGeneratorError


class PosGeneratorParallacticTest(unittest.TestCase):

    def setUp(self):
        self.m = mocker.Mocker()
        self.source = self.m.mock()
        self.posgen = PosGenerator(zdtimeout=2)

    def tearDown(self):
        self.m.restore()
        self.m.verify()

    def test_wrong_site_info(self):
        """Raise a PosGeneratorError in case of wrong site_info"""
        gen = self.posgen.parallactic(self.source, siteInfo={'LATITUDE': 39})
        self.assertRaises(PosGeneratorError, gen.next)
        gen = self.posgen.parallactic(self.source, siteInfo=11)
        self.assertRaises(PosGeneratorError, gen.next)

    def test_az_el_values(self):
        """Raise a PosGeneratorError if cannot get the (az, el) values"""
        source = self.m.mock()
        mocker.expect(source.getApparentCoordinates(mocker.ANY)).result((None, None))
        mocker.expect(source._get_name()).result('mocker')
        self.m.replay()
        gen = self.posgen.parallactic(source, siteInfo={'latitude': 39})
        self.assertRaises(PosGeneratorError, gen.next)

    def test_right_behavior(self):
        """Generate 3 positions and do not stop in case of isolated zero div"""
        azimuths = (45, 45, 90, 45)
        elevations = (45, 45, 0, 45) # The value 0 will cause a zero division
        latitude = 0
        source = self.m.mock()
        for (az, el) in zip(azimuths, elevations):
            mocker.expect(
                    source.getApparentCoordinates(mocker.ANY)
                    ).result((radians(az), radians(el)) + (None,)*5)
            self.m.count(1)
        self.m.replay()
        gen = self.posgen.parallactic(source, siteInfo={'latitude': latitude})
        latitude = radians(latitude)
        for azd, eld in zip(azimuths, elevations):
            az = radians(azd)
            el = radians(eld)
            p = atan2(-sin(az), (tan(latitude)*cos(el) - sin(el)*cos(az)))
            self.assertEqual(degrees(p), gen.next())

    def test_parallactic_greater_than_90(self):
        """Verify the (abs) parallactic angle can be greater than 90 degrees"""
        lat = 0.68940505453776013 # 39.5 degrees
        az = 6.1086523819801535 # 350 degrees
        el = 1.2217304763960306 # 70 degrees
        # The related parallactic angle must be 165 degrees
        parallactic = PosGenerator.getParallacticAngle(lat, az, el)
        self.assertGreater(abs(parallactic), 90)

if __name__ == '__main__':
    unittest.main()
