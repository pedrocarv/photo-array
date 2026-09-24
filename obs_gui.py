'''
GUI to view images from all fits files in a given directory

By JC

FIXME: loading all images at once may not be feasible for nights with many observations
(Very future) TODO: Would be nice to make this a web-app instead that can be accessed whenever
'''
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QPushButton, QHBoxLayout, QVBoxLayout, QStackedLayout, QLabel
import sys

from scip.data_processing.read_fits_files import *

from pathlib import Path
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('agg')

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure

# show_imgs func from SCIP_operation/util.py 
class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(20,30))
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122)
        super().__init__(self.fig)


# Get data from specified directory
# FIXME: Better way to do this than manually setting directory
# path = Path('/run/media/jc/Seagate Portable Drive/SCIP Data/4_26_25/')
# path = Path('/run/media/jc/Seagate Portable Drive/SCIP Data/20250602/')
path = input('enter file path to folder of interest: ')
start = (int)(input('enter image index to start at: '))

on_band, off_band = get_fits_from_folder(directory=path, open_image=True)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        
        self.max_ind = len(on_band) - 1

        if start < self.max_ind:
            self.ind = start
        else:
            self.ind = 0

        self.setWindowTitle('PATHS Fits Viewer')

        # Init buttons to navigate images
        btn_layout = QHBoxLayout()

        self.nxt = QPushButton('Next')
        self.nxt.clicked.connect(self.next_click)

        self.prv = QPushButton('Prev.')
        self.prv.clicked.connect(self.prev_click)
        if self.ind == 0:
            self.prv.setEnabled(False) # Init at start index, shouldn't be able to click

        btn_layout.addWidget(self.prv)
        btn_layout.addWidget(self.nxt)
        
        btn_wid = QWidget()
        btn_wid.setLayout(btn_layout)

        # Init labels to display image info (NOTE: On and off band should have same exact date and exp info, so only need to pull from one)
        self.t_obs = QLabel(on_band[self.ind]['DATE-OBS'])
        self.t_obs.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.t_exp = QLabel('Exp. Time: %.1f s'%(on_band[self.ind]['EXP-TIME']))
        self.t_exp.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.count = QLabel('Image index %d of %d'%(self.ind, self.max_ind))
        self.count.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # self.ims = self.show_ims()
        self.canv = MplCanvas(self)

        # Create GUI widget
        layout = QVBoxLayout()
        # layout.addWidget(self.ims)
        layout.addWidget(self.canv)
        layout.addWidget(self.count)
        layout.addWidget(self.t_obs)
        layout.addWidget(self.t_exp)
        layout.addWidget(btn_wid)

        wid = QWidget()
        wid.setLayout(layout)

        self.setCentralWidget(wid)
        self.update()
    
    def next_click(self):
        # Re-enable from min
        if self.ind == 0:
            self.prv.setEnabled(True)
    
        self.ind += 1
        # print(self.ind)
        self.t_obs.setText(on_band[self.ind]['DATE-OBS'])
        self.t_exp.setText('Exp. Time: %.1f s'%(on_band[self.ind]['EXP-TIME']))
        self.count.setText('Image index %d of %d'%(self.ind, self.max_ind))
        # self.ims = self.show_ims()
        # self.ims.draw()
        self.update()

        self.show()

        # Disable if max is reached
        if self.ind >= self.max_ind:
            self.nxt.setEnabled(False)
    
    def prev_click(self):
        # Re-enable if decrementing from max
        if self.ind == self.max_ind:
            self.nxt.setEnabled(True)

        self.ind -= 1
        # print(self.ind)
        self.t_obs.setText(on_band[self.ind]['DATE-OBS'])
        self.t_exp.setText('Exp. Time: %.1f s'%(on_band[self.ind]['EXP-TIME']))
        self.count.setText('Image index %d of %d'%(self.ind, self.max_ind))
        # self.ims = self.show_ims()
        # self.ims.draw()
        self.update()

        # Disable if this reaches 0
        if self.ind <= 0:
            self.prv.setEnabled(False)

    def show_ims(self):
        plt.clf()
        plt.cla()
        show_imgs(on_band[self.ind]['image'], off_band[self.ind]['image'])
        return plt.gcf().canvas
    
    def update(self):
        self.canv.ax1.imshow(on_band[self.ind]['image'], interpolation='none', vmin=np.percentile(on_band[self.ind]['image'], 1), vmax=np.percentile(on_band[self.ind]['image'], 99), cmap='gray')
        # self.canv.fig.colorbar(cax=self.canv.ax1)
        self.canv.ax1.set_title("Camera C1")

        self.canv.ax2.imshow(off_band[self.ind]['image'], interpolation='none', vmin=np.percentile(off_band[self.ind]['image'], 1), vmax=np.percentile(off_band[self.ind]['image'], 99), cmap='gray')
        # self.canv.fig.colorbar(cax=self.canv.ax2)
        self.canv.ax2.set_title("Camera C2")
        self.canv.draw_idle()


app = QApplication(sys.argv) #<- use commandline args to pass in directory of interest

window = MainWindow()
window.show()

app.exec()