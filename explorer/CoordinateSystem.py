from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class CoordinateSystem(QLabel):
    selected_changed = pyqtSignal(complex, name='selectedPointChanged')
    def __init__(self, parent=None):
        super(CoordinateSystem, self).__init__(parent)
        self.windowBL = (-4,-2)
        self.windowTR = (4,2)
        self.ratio = (self.windowTR[1]-self.windowBL[1])/(self.windowTR[0]-self.windowBL[0])
        self.resize(400,int(400*self.ratio))
        self.points = None
        self.setMinimumSize(1, 1);
        self.selectedPoint = 0

    def inWindow(self,z):
        (x,y) = (z.real,z.imag)
        return (self.windowBL[0] < x and x < self.windowTR[0]) and (self.windowBL[1] < y and y < self.windowTR[1])

    def complexToWindowCoords(self,z):
        return ((z.real-self.windowBL[0]) * self.width()/(self.windowTR[0]-self.windowBL[0]), self.height()*(1-(z.imag-self.windowBL[1])/(self.windowTR[1]-self.windowBL[1])))
    def windowCoordsToComplex(self,x,y):
        re = (x * (self.windowTR[0]-self.windowBL[0])/self.width()) + self.windowBL[0]
        im = (1-y/self.height())*(self.windowTR[1]-self.windowBL[1]) + self.windowBL[1]
        return re + 1j * im

    def paintPoints(self, points):
        self.points = points
        pixmap = QPixmap(self.size())
        pixmap.fill()
        painter = QPainter(pixmap)
        painter.setPen(QPen(QColor(0,0,0,100), int(self.width()/400)))
        #painter.scale(self.width() / 100.0, self.height() / 100.0)

        for point in points:
            if self.inWindow(point):
                coords = self.complexToWindowCoords(point)
                painter.drawPoint(int(coords[0]),int(coords[1]))

        painter.end()
        self.setPixmap(pixmap)

        if self.selectedPoint != None:
            self.selectPoint(self.selectedPoint)

    def selectPoint(self,z):
        if self.selectedPoint != z:
            self.selectedPoint = z
            self.selected_changed.emit(complex(z))
        if self.inWindow(z):
            pixmap = self.pixmap()
            painter = QPainter(pixmap)
            painter.setPen(QPen(QColor(255,0,0), int(self.width()/100)))
            coords = self.complexToWindowCoords(self.selectedPoint)
            painter.drawPoint(int(coords[0]),int(coords[1]))
            painter.end()
            self.setPixmap(pixmap)

    def resizeEvent(self,e):
        super().resizeEvent(e)
        if self.points != None:
            self.resize(self.width(),int(self.width()*self.ratio))
            self.paintPoints(self.points)

    def mouseMoveEvent(self,e):
        #super().mouseMoveEvent(e) # Disable this, otherwise the drag event gets sent to KDE and the whole window moves
        self.selectPoint(self.windowCoordsToComplex(e.x(),e.y()))

