"""Consistent manuscript-inspection style for candidate 693168 figures."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

COLORS={"highT_reference":"#CC79A7","success":"#009E73","lower_failure":"#0072B2","upper_failure":"#D55E00","lowT_isothermal":"#666666"}
RATES={1:"#000000",20:"#0072B2",50:"#009E73",100:"#D55E00"}
CLASSES={"DENSIFICATION_EXHAUSTION_FAILURE":"#3B75AF","SUCCESS":"#2CA25F","GRAIN_GROWTH_FAILURE":"#E6550D","MIXED_FAILURE":"#984EA3","UNATTAINABLE_FIRST_STEP":"#BDBDBD","TARGET_REACHED_DURING_FIRST_STEP":"#636363","NUMERICAL_CENSOR":"#FDD0A2"}

def apply():
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.labelsize":11,"axes.titlesize":11,"legend.fontsize":8,"xtick.labelsize":9,"ytick.labelsize":9,"lines.linewidth":2,"axes.linewidth":.8,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"savefig.dpi":600})
def clean(ax):ax.grid(True,alpha=.18);ax.tick_params(direction="out",length=3)
def letters(axes):
    for i,ax in enumerate(getattr(axes,"flat",[axes])):ax.text(-.13,1.04,chr(65+i),transform=ax.transAxes,fontweight="bold",fontsize=13)
def save(fig,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(path.with_suffix('.png'),dpi=600,bbox_inches='tight');plt.close(fig)
