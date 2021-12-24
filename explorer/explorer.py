""" A better version of graphical_limits.py, a GUI program to display limit sets of groups in the Riley slices.

    Inspired by the schottky software written by Danny Calegari and Alden Walker (https://github.com/dannycalegari/schottky).
"""

import riley
import mpmath as mp

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
from PyQt5.QtGui import QPixmap, QPainter
from ui_explorer import Ui_MainWindow

scale = 100
riley_bounds = (-4,4,-4,4) # -x,x,-y,y
limit_bounds = (-4,4,-4,4) # -x,x,-y,y

def canvas_to_usual_coords():
    y = scale*(riley_bounds[3]-riley_bounds[2]) - y
    x = x/scale + riley_bounds[0]
    y = y/scale + riley_bounds[2]
    return x,y

def usual_coords_to_canvas(x,y):
    x = scale*(x - riley_bounds[0])
    y = scale*(y - riley_bounds[2])
    y = scale*(riley_bounds[3]-riley_bounds[2]) - y
    return x,y

def show_about_dialog():
    text = "<center>" \
           "<h1>Riley Slice Explorer</h1>" \
           "</center" \
           "<p><a href=\"https://github.com/aelzenaar/riley\">https://github.com/aelzenaar/riley</a><br/>" \
           "Copyright &copy; Alex Elzenaar.</p>"
    QMessageBox.about(window, "About Explorer", text)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(window)

    ui.actionAbout.triggered.connect(show_about_dialog)
    ui.sliceView.paintPoints(riley.riley_slice(3,4,10))


    window.show()
    sys.exit(app.exec_())
