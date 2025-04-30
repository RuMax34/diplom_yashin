import numpy as np
import matplotlib.pyplot as plt

# Enable LaTeX rendering in Matplotlib
#from matplotlib.pyplot import rc
#rc('text', usetex=True)
#rc('font', family='Times')
#rc('text.latex',unicode=True)
#rc('text.latex', preamble=r'\usepackage[utf8]{inputenc}')
#rc('text.latex', preamble=r'\usepackage[russian]{babel}')

# Enable regular rendering in Matplotlib (without LaTeX)
from matplotlib import rcParams
# Set the font globally
rcParams['font.family'] = 'Times New Roman'  # Replace with your font name
rcParams['mathtext.fontset'] = 'stix'

# Colors
c = [
    (68, 119, 170),  # good blue
    (238, 119, 51),  # orange
    (34, 136, 51),   # green
    (51, 34, 136),   # indigo, dark purple
    (204, 51, 17),   # red
    (0, 153, 136),   # teal
    (221, 170, 51),  # good yellow
    (80, 80, 80),    # dark gray
    (170, 51, 119),  # purple, cherry
    (102, 37, 6),    # brown
    (102, 204, 238), # sky blue
    (204, 187, 68),  # lemon
    (187, 187, 187), # gray
    (0, 68, 136),    # darker blue
    (153, 153, 51),  # olive
    (136, 34, 85),   # wine
    (0, 0, 0),       # black
]

# Dash types
d = [
    'solid',
    'dashed',
    'dotted',
    'dashdot',
    (5, (10, 3)),
    'solid',
    'dashed',
    'dotted',
    'dashdot',
    (5, (10, 3)),
    'solid',
    'dashed',
    'dotted',
    'dashdot',
    (5, (10, 3)),
]

def col(i):
    return tuple(rgb / 255 for rgb in c[i])

# Read the file
data = np.loadtxt('results/droplet.dat')
x = data[:, 0]
y1 = data[:, 1]

# Desired figure size in centimeters
def inch(width, height):
    return (width / 2.54, height / 2.54)

# =========================
# Rectangular (full page) plot
# =========================

width = 15
height = 10
lwd = 1.5
bwd = 0.75
msz = 2
fsz = 12
scale = 1.5

# Create a figure with custom size and layout
fig, ax = plt.subplots(figsize=inch(width, height))  # Figure size in inches (width, height)

# Plot the data
ax.plot(x, y1, label=r'$R_d(t)/R_d(0)$', color=col(0), linestyle=d[0], linewidth=lwd, markersize=msz)

# Customize the axes
ax.set_xlabel(r'$t$, наносекунды', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_ylabel(r'Функции, отн.ед.', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_title('Зависимость радиуса капли и потока Ga от времени', fontsize=fsz, fontweight='normal', pad=5*scale)
ax.tick_params(axis='both', which='major', labelsize=fsz-2, length=2*scale, width=bwd, direction='out')
ax.grid(True, which='both', linestyle='dotted', linewidth=bwd, alpha=1)

# Change the box (spines) line width
for spine in ax.spines.values():
    spine.set_linewidth(bwd)
plt.tick_params(axis='y', labelrotation=90)

# Add a legend
ax.legend(
    fontsize=fsz,
    loc='lower left',
    frameon=True,
    framealpha=1,
    fancybox=True,
    edgecolor='black',
    borderaxespad=0.5*scale,
    handleheight=0.25*scale,
    labelspacing=0.25*scale,
).get_frame().set_linewidth(bwd)

# Adjust subplot parameters for better spacing
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)

# Save the plot to a PDF file
plt.savefig("droplet.pdf", format="pdf", bbox_inches="tight")

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Read the file
data = np.loadtxt('results/height.dat')
num_columns = data.shape[1]
x = data[:, num_columns-1]
y1 = data[:, 0]
y2 = data[:, 1]
y3 = data[:, 2]
y4 = data[:, 3]
y5 = data[:, 4]
y6 = data[:, 5]
y7 = data[:, 6]
y8 = data[:, 7]
y9 = data[:, 8]
y10 = data[:, 9]

# Desired figure size in centimeters
def inch(width, height):
    return (width / 2.54, height / 2.54)

# =========================
# Rectangular (full page) plot
# =========================

width = 15
height = 10
lwd = 1.5
bwd = 0.75
msz = 2
fsz = 12
scale = 1.5

# Create a figure with custom size and layout
fig, ax = plt.subplots(figsize=inch(width, height))  # Figure size in inches (width, height)

# Plot the data
ax.plot(x, y1, label=r't=1', color=col(0), linestyle=d[0], linewidth=lwd, markersize=msz)
ax.plot(x, y2, label=r't=2', color=col(1), linestyle=d[1], linewidth=lwd, markersize=msz)
ax.plot(x, y3, label=r't=3', color=col(2), linestyle=d[2], linewidth=lwd, markersize=msz)
ax.plot(x, y4, label=r't=4', color=col(3), linestyle=d[3], linewidth=lwd, markersize=msz)
ax.plot(x, y5, label=r't=5', color=col(4), linestyle=d[4], linewidth=lwd, markersize=msz)
ax.plot(x, y6, label=r't=6', color=col(5), linestyle=d[5], linewidth=lwd, markersize=msz)
ax.plot(x, y7, label=r't=7', color=col(6), linestyle=d[6], linewidth=lwd, markersize=msz)
ax.plot(x, y8, label=r't=8', color=col(7), linestyle=d[7], linewidth=lwd, markersize=msz)
ax.plot(x, y9, label=r't=9', color=col(8), linestyle=d[8], linewidth=lwd, markersize=msz)
ax.plot(x, y10, label=r't=10', color=col(9), linestyle=d[9], linewidth=lwd, markersize=msz)

# Customize the axes
ax.set_xlabel(r'$r$, нм', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_ylabel(r'$h(r)$, нм', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_title('Профиль КК в разные моменты времени', fontsize=fsz, fontweight='normal', pad=5*scale)
ax.tick_params(axis='both', which='major', labelsize=fsz-2, length=2*scale, width=bwd, direction='out')
ax.grid(True, which='both', linestyle='dotted', linewidth=bwd, alpha=1)

# Change the box (spines) line width
for spine in ax.spines.values():
    spine.set_linewidth(bwd)
plt.tick_params(axis='y', labelrotation=90)

# Add a legend
ax.legend(
    fontsize=fsz,
    bbox_to_anchor=(1.05, 0),
    loc='lower left',
    frameon=True,
    framealpha=1,
    fancybox=True,
    edgecolor='black',
    borderaxespad=0.5*scale,
    handleheight=0.25*scale,
    labelspacing=0.25*scale,
).get_frame().set_linewidth(bwd)

# Adjust subplot parameters for better spacing
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)

# Save the plot to a PDF file
plt.savefig("height_profile.pdf", format="pdf", bbox_inches="tight")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Read the file
data = np.loadtxt('results/c_ga.dat')
num_columns = data.shape[1]
x = data[:, num_columns-1]
y1 = data[:, 0]
y2 = data[:, 1]
y3 = data[:, 2]
y4 = data[:, 3]
y5 = data[:, 4]
y6 = data[:, 5]
y7 = data[:, 6]
y8 = data[:, 7]
y9 = data[:, 8]
y10 = data[:, 9]

# Desired figure size in centimeters
def inch(width, height):
    return (width / 2.54, height / 2.54)

# =========================
# Rectangular (full page) plot
# =========================

width = 15
height = 10
lwd = 1.5
bwd = 0.75
msz = 2
fsz = 12
scale = 1.5

# Create a figure with custom size and layout
fig, ax = plt.subplots(figsize=inch(width, height))  # Figure size in inches (width, height)

# Plot the data
ax.plot(x, y1, label=r't=1', color=col(0), linestyle=d[0], linewidth=lwd, markersize=msz)
ax.plot(x, y2, label=r't=2', color=col(1), linestyle=d[1], linewidth=lwd, markersize=msz)
ax.plot(x, y3, label=r't=3', color=col(2), linestyle=d[2], linewidth=lwd, markersize=msz)
ax.plot(x, y4, label=r't=4', color=col(3), linestyle=d[3], linewidth=lwd, markersize=msz)
ax.plot(x, y5, label=r't=5', color=col(4), linestyle=d[4], linewidth=lwd, markersize=msz)
ax.plot(x, y6, label=r't=6', color=col(5), linestyle=d[5], linewidth=lwd, markersize=msz)
ax.plot(x, y7, label=r't=7', color=col(6), linestyle=d[6], linewidth=lwd, markersize=msz)
ax.plot(x, y8, label=r't=8', color=col(7), linestyle=d[7], linewidth=lwd, markersize=msz)
ax.plot(x, y9, label=r't=9', color=col(8), linestyle=d[8], linewidth=lwd, markersize=msz)
ax.plot(x, y10, label=r't=10', color=col(9), linestyle=d[9], linewidth=lwd, markersize=msz)

# Customize the axes
ax.set_xlabel(r'$r$, нм', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_ylabel(r'$С(r)$, усл.ед.', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_title('Концентрация Ga в разные моменты времени', fontsize=fsz, fontweight='normal', pad=5*scale)
ax.tick_params(axis='both', which='major', labelsize=fsz-2, length=2*scale, width=bwd, direction='out')
ax.grid(True, which='both', linestyle='dotted', linewidth=bwd, alpha=1)

# Change the box (spines) line width
for spine in ax.spines.values():
    spine.set_linewidth(bwd)
plt.tick_params(axis='y', labelrotation=90)

# Add a legend
ax.legend(
    fontsize=fsz,
    bbox_to_anchor=(1.05, 0),
    loc='lower left',
    frameon=True,
    framealpha=1,
    fancybox=True,
    edgecolor='black',
    borderaxespad=0.5*scale,
    handleheight=0.25*scale,
    labelspacing=0.25*scale,
).get_frame().set_linewidth(bwd)

# Adjust subplot parameters for better spacing
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)

# Save the plot to a PDF file
plt.savefig("ga_concentration.pdf", format="pdf", bbox_inches="tight")


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Read the file
data = np.loadtxt('results/c_as.dat')
num_columns = data.shape[1]
x = data[:, num_columns-1]
y1 = data[:, 0]
y2 = data[:, 1]
y3 = data[:, 2]
y4 = data[:, 3]
y5 = data[:, 4]
y6 = data[:, 5]
y7 = data[:, 6]
y8 = data[:, 7]
y9 = data[:, 8]
y10 = data[:, 9]

# Desired figure size in centimeters
def inch(width, height):
    return (width / 2.54, height / 2.54)

# =========================
# Rectangular (full page) plot
# =========================

width = 15
height = 10
lwd = 1.5
bwd = 0.75
msz = 2
fsz = 12
scale = 1.5

# Create a figure with custom size and layout
fig, ax = plt.subplots(figsize=inch(width, height))  # Figure size in inches (width, height)

# Plot the data
ax.plot(x, y1, label=r't=1', color=col(0), linestyle=d[0], linewidth=lwd, markersize=msz)
ax.plot(x, y2, label=r't=2', color=col(1), linestyle=d[1], linewidth=lwd, markersize=msz)
ax.plot(x, y3, label=r't=3', color=col(2), linestyle=d[2], linewidth=lwd, markersize=msz)
ax.plot(x, y4, label=r't=4', color=col(3), linestyle=d[3], linewidth=lwd, markersize=msz)
ax.plot(x, y5, label=r't=5', color=col(4), linestyle=d[4], linewidth=lwd, markersize=msz)
ax.plot(x, y6, label=r't=6', color=col(5), linestyle=d[5], linewidth=lwd, markersize=msz)
ax.plot(x, y7, label=r't=7', color=col(6), linestyle=d[6], linewidth=lwd, markersize=msz)
ax.plot(x, y8, label=r't=8', color=col(7), linestyle=d[7], linewidth=lwd, markersize=msz)
ax.plot(x, y9, label=r't=9', color=col(8), linestyle=d[8], linewidth=lwd, markersize=msz)
ax.plot(x, y10, label=r't=10', color=col(9), linestyle=d[9], linewidth=lwd, markersize=msz)

# Customize the axes
ax.set_xlabel(r'$r$, нм', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_ylabel(r'$С(r)$, усл.ед.', fontsize=fsz, fontweight='normal', labelpad=5*scale)
ax.set_title('Концентрация As в разные моменты времени', fontsize=fsz, fontweight='normal', pad=5*scale)
ax.tick_params(axis='both', which='major', labelsize=fsz-2, length=2*scale, width=bwd, direction='out')
ax.grid(True, which='both', linestyle='dotted', linewidth=bwd, alpha=1)

# Change the box (spines) line width
for spine in ax.spines.values():
    spine.set_linewidth(bwd)
plt.tick_params(axis='y', labelrotation=90)

# Add a legend
ax.legend(
    fontsize=fsz,
    bbox_to_anchor=(1.05, 0),
    loc='lower left',
    frameon=True,
    framealpha=1,
    fancybox=True,
    edgecolor='black',
    borderaxespad=0.5*scale,
    handleheight=0.25*scale,
    labelspacing=0.25*scale,
).get_frame().set_linewidth(bwd)

# Adjust subplot parameters for better spacing
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)

# Save the plot to a PDF file
plt.savefig("as_concentration.pdf", format="pdf", bbox_inches="tight")