""" A better version of graphical_limits.py, a GUI program to display limit sets of groups in the Riley slices.

    Inspired by the schottky software written by Danny Calegari and Alden Walker (https://github.com/dannycalegari/schottky).
"""

import riley
import farey
import kleinian
import mpmath as mp

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
from PyQt5.QtGui import QPixmap, QPainter, QIntValidator, QDoubleValidator, QImage
from ui_explorer import Ui_MainWindow

app = QApplication(sys.argv)
window = QMainWindow()
ui = Ui_MainWindow()
ui.setupUi(window)
pOrder = None
qOrder = None

def show_about_dialog():
    text = "<center>" \
           "<h1>Riley Slice Explorer</h1>" \
           "</center" \
           "<p><a href=\"https://github.com/aelzenaar/riley\">https://github.com/aelzenaar/riley</a><br/>" \
           "Copyright &copy; Alex Elzenaar.</p>"
    QMessageBox.about(window, "About Explorer", text)

def show_posint_validation_failed_dialog():
    QMessageBox.critical(window, "Invalid input", "Only enter positive numbers")
def show_real_validation_failed_dialog():
    QMessageBox.critical(window, "Invalid input", "Only enter real numbers")

def slice_parameters_changed(_=None):
    global pOrder, qOrder
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
        qOrder = int(ui.qOrderEdit.text())

    fareyDenom = int(ui.fareyDenomEdit.text())

    ui.sliceView.paintPoints(riley.riley_slice(pOrder,qOrder,fareyDenom))

def slice_point_changed_via_edit():
    ui.sliceView.selectPoint(float(ui.muReEdit.text()) + 1j * float(ui.muImEdit.text()))

def slice_point_changed_via_click(z):
    ui.muReEdit.setText(str(z.real)[:8])
    ui.muImEdit.setText(str(z.imag)[:8])

def limit_parameters_changed(z):
    ui.limitView.redrawLimitSet(pOrder,qOrder,z)


if __name__ == '__main__':
    ui.actionAbout.triggered.connect(show_about_dialog)

    ui.pOrderRadios.buttonClicked.connect(slice_parameters_changed)
    ui.qOrderRadios.buttonClicked.connect(slice_parameters_changed)

    order_validator = QIntValidator(0,100000)
    mu_validator = QDoubleValidator()

    def setup_edit(edit,changed,validator,fail_function):
        edit.editingFinished.connect(changed)
        edit.setValidator(validator)
        edit.inputRejected.connect(fail_function)

    setup_edit(ui.pOrderEdit,slice_parameters_changed,order_validator,show_posint_validation_failed_dialog)
    setup_edit(ui.qOrderEdit,slice_parameters_changed,order_validator,show_posint_validation_failed_dialog)
    setup_edit(ui.fareyDenomEdit,slice_parameters_changed,order_validator,show_posint_validation_failed_dialog)
    setup_edit(ui.muReEdit,slice_point_changed_via_edit,mu_validator,show_real_validation_failed_dialog)
    setup_edit(ui.muImEdit,slice_point_changed_via_edit,mu_validator,show_real_validation_failed_dialog)
    ui.sliceView.selected_changed.connect(slice_point_changed_via_click)
    ui.sliceView.selected_changed.connect(limit_parameters_changed)

    slice_parameters_changed()
    limit_parameters_changed(2j)



    window.show()
    sys.exit(app.exec_())
