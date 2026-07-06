import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')
colors =  plt.colormaps["tab10"]
x = np.linspace(-4, 4, 100000)


zeta_1 = np.minimum(abs(1/(x+1e-20)), 1)
zeta_2 = 1/(np.sqrt(1+x**2))
zeta_3 = 1/(np.sqrt(0.25+x**2))
zeta_4 = np.tanh(1/(abs(x)))

labels = [r'max($1/|x|,1$)',r'$1/\sqrt{1+x^2}$',r'$1/\sqrt{0.25+x^2}$',r'$\tanh(1/|x|)$']


fig = plt.figure()

ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(x,zeta_1, color = colors(0), label = labels[0])
plt.plot(x,zeta_2, color = colors(1), label = labels[1])
plt.plot(x,zeta_3, color = colors(2), label = labels[2])
plt.plot(x,zeta_4, color = colors(3), label = labels[3])

plt.axvspan(-1,1,color = 'royalblue',alpha = 0.25, lw = 0, edgecolor = None)
plt.axvspan(-4,-1,color = 'crimson',alpha = 0.25, lw = 0, edgecolor = None)
plt.axvspan(1,4,color = 'crimson',alpha = 0.25, lw = 0, edgecolor = None)
# ax.spines['bottom'].set_alpha(0.3)

plt.xlabel('Strain rate',fontsize = 14, alpha = 0.5)
ax.annotate("Viscosity", xy=(-0.34,1), xytext=(0,3), 
            xycoords="axes fraction", textcoords="offset points", 
            ha="left", alpha = 0.5, fontsize = 14) 

# Change all spines at once
for spine in ax.spines.values():
    spine.set_alpha(0.5)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_alpha(0.5)
 
# for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
#     tick.tick1line.set_alpha(0.5)
#     tick.tick2line.set_alpha(0.5)
    

# plt.subplots_adjust(left=0.25, top=0.85)
fig.legend(labelcolor = 'linecolor',handlelength=0, bbox_to_anchor = (1.18,0.72))
# plt.tight_layout()

plt.savefig('regularizations.png', bbox_inches='tight')

def make_figure(n_lines):
    fig = plt.figure(facecolor='none')
    ax = plt.axes()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    all_lines = [
        (x, zeta_1, colors(1), labels[0]),
        (x, zeta_2, colors(2), labels[1]),
        (x, zeta_3, colors(3), labels[2]),
        (x, zeta_4, colors(4), labels[3]),
    ]

    for i, (xd, yd, c, l) in enumerate(all_lines):
        if i < n_lines:
            ax.plot(xd, yd, color=c, label=l)        # visible, in legend
        else:
            ax.plot(xd, yd, color=c, alpha=0)         # invisible, no label → not in legend

    ax.axvspan(-1, 1, color='crimson', alpha=0.25, lw=0, edgecolor=None)
    ax.axvspan(-4, -1, color='royalblue', alpha=0.25, lw=0, edgecolor=None)
    ax.axvspan(1, 4, color='royalblue', alpha=0.25, lw=0, edgecolor=None)

    plt.xlabel('Strain rate', fontsize=14, alpha=0.5)
    ax.annotate("Viscosity", xy=(-0.34, 1), xytext=(0, 3),
                xycoords="axes fraction", textcoords="offset points",
                ha="left", alpha=0.5, fontsize=14)

    for spine in ax.spines.values():
        spine.set_alpha(0.5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_alpha(0.5)

    # Always plot all lines with alpha=0 on a twin axis to fix the axes limits
    ax2 = ax.twinx()
    ax2.set_visible(False)
    for (xd, yd, c, l) in all_lines:
        ax2.plot(xd, yd, alpha=0)

    fig.legend(labelcolor='linecolor', handlelength=0, bbox_to_anchor=(0.9, 0.92))
    plt.savefig(f'regularizations_{n_lines}.png', bbox_inches='tight')
    plt.close()

for n in range(1, 5):
    make_figure(n)