from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import *

class CoordinateSystem(QLabel):
    def __init__(self, parent=None):
        super(CoordinateSystem, self).__init__(parent)
        self.windowBL = (-4,-4)
        self.windowTR = (4,4)
        self.resize(200,200)

    def inWindow(self,z):
        (x,y) = (z.real,z.imag)
        return (self.windowBL[0] < x and x < self.windowTR[0]) and (self.windowBL[1] < y and y < self.windowTR[1])

    def complexToWindowCoords(self,z):
        return ((z.real-self.windowBL[0]) * self.width()/(self.windowTR[0]-self.windowBL[0]), self.height()*(1-(z.imag-self.windowBL[1])/(self.windowTR[1]-self.windowBL[1])))

    def paintPoints(self, points):
        pixmap = QPixmap(self.size())
        pixmap.fill()
        painter = QPainter(pixmap)
        painter.setPen(QPen(QColor(0), 3))
        #painter.scale(self.width() / 100.0, self.height() / 100.0)

        for point in points:
            if self.inWindow(point):
                coords = self.complexToWindowCoords(point)
                painter.drawPoint(int(coords[0]),int(coords[1]))

        painter.end()
        self.setPixmap(pixmap)
