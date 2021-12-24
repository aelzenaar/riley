""" A better version of graphical_limits.py, a GUI program to display limit sets of groups in the Riley slices.

    Inspired by the schottky software written by Danny Calegari and Alden Walker (https://github.com/dannycalegari/schottky).
"""

import riley
import mpmath as mp

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
from PyQt5.QtGui import QPixmap, QPainter, QIntValidator
from ui_explorer import Ui_MainWindow

scale = 100
riley_bounds = (-4,4,-4,4) # -x,x,-y,y
limit_bounds = (-4,4,-4,4) # -x,x,-y,y

app = QApplication(sys.argv)
window = QMainWindow()
ui = Ui_MainWindow()
ui.setupUi(window)


def show_about_dialog():
    text = "<center>" \
           "<h1>Riley Slice Explorer</h1>" \
           "</center" \
           "<p><a href=\"https://github.com/aelzenaar/riley\">https://github.com/aelzenaar/riley</a><br/>" \
           "Copyright &copy; Alex Elzenaar.</p>"
    QMessageBox.about(window, "About Explorer", text)

def show_validation_failed_dialog():
    QMessageBox.critical(window, "Invalid input", "Only enter positive numbers")

def slice_parameters_changed(_=None):
    if ui.pInfRadioButton.isChecked():
        pOrder = mp.inf
    else:
        if ui.pOrderEdit.text() == '':
            return
        pOrder = int(ui.pOrderEdit.text())

    if ui.qInfRadioButton.isChecked():
        qOrder = mp.inf
    else:
        if ui.qOrderEdit.text() == '':
            return
        qOrder = int(ui.pOrderEdit.text())

    fareyDenom = int(ui.fareyDenomEdit.text())

    ui.sliceView.paintPoints(riley.riley_slice(pOrder,qOrder,fareyDenom))


if __name__ == '__main__':
    ui.actionAbout.triggered.connect(show_about_dialog)

    ui.pOrderRadios.buttonClicked.connect(slice_parameters_changed)
    ui.qOrderRadios.buttonClicked.connect(slice_parameters_changed)

    order_validator = QIntValidator(0,100000)
    def setup_edit(edit):
        edit.editingFinished.connect(slice_parameters_changed)
        edit.setValidator(order_validator)
        edit.inputRejected.connect(show_validation_failed_dialog)
    setup_edit(ui.pOrderEdit)
    setup_edit(ui.qOrderEdit)
    setup_edit(ui.fareyDenomEdit)

    slice_parameters_changed()

    ui.limitView.paintPoints(riley.riley_slice(4,3,40))


    window.show()
    sys.exit(app.exec_())
