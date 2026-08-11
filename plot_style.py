"""Shared journal-style plotting utilities for the manuscript figure package."""
from pathlib import Path
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt

COLORS={"slow":"#0072B2","fast":"#D55E00","highT":"#CC79A7","two_step":"#009E73","positive":"#009E73","neutral":"#999999","negative":"#D55E00","unattainable":"#6A3D9A","q0":"#56B4E9","q1":"#E69F00"}
TJ_COLORS={"current_all_TJ_multihit":"#4C78A8","pore_relaxed_TJ":"#72B7B2","pore_pinned_drag":"#F58518","mixed_relaxed_pinned":"#B279A2"}

LABELS={"rho":r"Relative density, $\rho$","G":r"Grain size, $G$ [nm]","radius":r"Connected mean pore radius, $\bar{r}_{p}^{\,c}$ [nm]","fine":r"Connected fine-pore fraction, $f_{\mathrm{fine}}^{c}$","WPR":r"PR/de-sintering work, $W_{\mathrm{PR}}$ [model units]","XJ":r"Persistent junction state, $X_J$","lambdaK":r"TJ activity ratio, $\Lambda_{\mathrm{TJ}}/K_{\mathrm{TJ}}$"}

def apply_style():
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.labelsize":11,"axes.titlesize":11,"xtick.labelsize":9,"ytick.labelsize":9,"legend.fontsize":9,"figure.dpi":120,"savefig.dpi":600,"pdf.fonttype":42,"ps.fonttype":42,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.8,"lines.linewidth":2.2,"lines.markersize":6,"grid.alpha":.2,"grid.linewidth":.6})

def panel_labels(axes):
    axes=list(getattr(axes,"flat",[axes]))
    for i,ax in enumerate(axes):ax.text(-.12,1.04,f"({chr(97+i)})",transform=ax.transAxes,fontweight="bold",fontsize=13,va="bottom")

def finish(fig,outbase):
    outbase=Path(outbase);outbase.parent.mkdir(parents=True,exist_ok=True);fig.savefig(outbase.with_suffix(".pdf"),bbox_inches="tight");fig.savefig(outbase.with_suffix(".png"),dpi=600,bbox_inches="tight");plt.close(fig)

def write_inventory(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);fields=("figure_id","filename_pdf","filename_png","short_title","source_table_or_script","purpose","manuscript_location")
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def inventory_row(fid,stem,title,source,purpose,location):
    return dict(figure_id=fid,filename_pdf=f"{stem}.pdf",filename_png=f"{stem}.png",short_title=title,source_table_or_script=source,purpose=purpose,manuscript_location=location)

def clean(ax,grid=True):
    if grid:ax.grid(True)
    ax.tick_params(direction="out",length=3.5,width=.8)
