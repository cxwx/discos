from ACSImpl.DevIO import DevIO

class GenericDevIO(DevIO):
    def __init__(self, value=0):
        DevIO.__init__(self, value)

    def read(self):
        return self.value

    def write(self, value):
        self.value = value


class amplitudeDevIO(DevIO):
    def __init__(self, value=0):
        DevIO.__init__(self, value)

    def read(self):
        return self.value

    def write(self, value):
        self.value = value
                                    
class frequencyDevIO(DevIO):
    def __init__(self, cl):
        DevIO.__init__(self, cl)
        self.cl=cl
    def read(self):
        self.cl.getFrequency()
        return self.value

    def write(self, value):
        self.value = value
class isLockedDevIO(DevIO):
    def __init__(self, value=0):
        DevIO.__init__(self, value)

    def read(self):
        return self.value

    def write(self, value):
        self.value = value
