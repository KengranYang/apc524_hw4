####### INSTRUCTION #######
# to use it, simply input the file name in the variable "file_name", and type "python contour_plot.py" in the terminal


from matplotlib.pyplot import *
import getData as gd
from numpy import *

file_name = ['/Users/hty/Google Drive/grad course work/apc524 software engineering/hw4/apc524_hw4/example.txt']

Readfiles=[file(file_name[0])]

interm_syn = [array(Readfiles[0].readlines()) for i in range (shape(Readfiles)[0])]

''' remember to set start line number here '''
interm1_syn =[ [str.split(interm_syn[j][i]) for i in range(1,shape(interm_syn[j])[0])] for j in range (shape(interm_syn)[0])]

xy= [[[float(interm1_syn[k][j][i]) for j in range (0,shape(interm1_syn[k])[0])]for i in range (shape(interm1_syn[k])[0]) ] for k in range (shape(interm_syn)[0])]
xy = transpose(xy[0])
close('all')
imshow(xy)
colorbar()

show()
