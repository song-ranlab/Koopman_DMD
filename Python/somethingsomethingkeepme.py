#%% md

# Not my animator but im throwing a Hail Mary something finally works
x = np.sin(y0[:,0])
y = -np.cos(y0[:,0])

from matplotlib import animation
from IPython.display import HTML
FPS=30
plt.style.use('default')

# Set up the figure
fig = plt.figure(figsize=(4, 4), dpi=60)
ax = plt.axes(xlim=(-7.0, 7.0), ylim=(-7.0, 7.0))
ax.set_aspect('equal')
ax.axis('off')

# Define the different elements in the animation
rod, = ax.plot([], [], color="grey", linewidth=2)
ball = plt.Circle((x[0], y[0]), 0.1, fc="grey")
ax.add_patch(ball)

# Calculates the number of frames
framesNum = int(FPS*t_span[-1])

# Animation function. This is called sequentially.
def animate(j):
    i = j*int(len(t_span)/framesNum)
    ball.center = (x[i], y[i])
    rod.set_data([0, x[i]], [0, y[i]])

# Create animation
anim = animation.FuncAnimation(fig, animate, frames=framesNum, interval=1000/FPS)

plt.close(anim._fig)

# Display the animation
#HTML(anim.to_html5_video())