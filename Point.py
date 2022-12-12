class Point:

    def __init__(self, a, b):
        self.latitude = a
        self.longitude = b

    def getLatitude(self):
        return self.latitude

    def getLongitude(self):
        return self.longitude

    def setLatitude(self, x):
        self.latitude = x

    def setLongitude(self, y):
        self.longitude = y

    def asArray(self):
        return [self.latitude, self.longitude]

    def toString(self):
        return "{} latitude {} longitude".format(self.latitude, self.longitude)
