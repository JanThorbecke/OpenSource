from scitools.std import *

x = y = linspace(1, 100, 100)

xv, yv = ndgrid(x, y)

#values = sin(sqrt(xv**2 + yv**2))
values = fromfile("corr2d_scl.bin", dtype=single, count=-1, sep='')
values.shape=(100,100)

#mesh(xv, yv, values)
surf(xv, yv, values,
     shading='flat',
     colorbar='on',
     colormap=gray(m=1),
	axis=[1,100,1,100,0,1],
     view=[35,45])

hardcopy('tmp0.eps')
