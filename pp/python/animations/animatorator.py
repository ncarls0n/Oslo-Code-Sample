import matplotlib.pyplot as plt
import matplotlib.image as mgimg
from matplotlib import animation
from matplotlib.animation import PillowWriter

fig = plt.figure(figsize=(15,15))
plt.axis('off')

myimages = []

for p in range(1,23):
    fname = 'Figure_'+str(p)+'.png'
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img, aspect='auto')
    imgplot.axes.get_xaxis().set_visible(False)
    imgplot.axes.get_xaxis().set_visible(False)
    myimages.append([imgplot])

my_anim = animation.ArtistAnimation(fig,myimages,interval=1000,blit=True,repeat_delay=2000)

# plt.savefig('pict.png', bbox_inches='tight', pad_inches = 0)

writer = PillowWriter(fps=2)
my_anim.save('test.gif', writer=writer)

#plt.show()
