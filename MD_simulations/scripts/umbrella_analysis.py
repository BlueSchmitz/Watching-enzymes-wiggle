### This script plots histograms from .xvg files ###

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from tqdm import tqdm

# ---------------------------
# Thesis plotting style
# ---------------------------

cmap = plt.get_cmap("viridis")
main_color = cmap(0.5)  # default color for single-line plots

plt.rcParams.update({
    "figure.dpi": 200,
    "image.cmap": "viridis",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "lines.linewidth": 2,
})

# ---------------------------


def load_xvg_columns(path, cols=(0,1)):
    """Fast XVG loader reading only needed columns."""
    return np.loadtxt(path, comments=['@','#'], usecols=cols, unpack=True)


def find_file(folder, suffix):
    for f in os.listdir(folder):
        if f.endswith(suffix) and "#" not in f:
            return os.path.join(folder, f)
    raise FileNotFoundError(f"{suffix} not found in {folder}")


def plot_line(x, y, xlabel, ylabel, title, path):
    plt.figure(figsize=(8,5))
    plt.plot(x,y,color=main_color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def plot_hist(data, xlabel, title, path):
    plt.figure(figsize=(8,5))
    plt.hist(data, bins=50, color=main_color)
    plt.xlabel(xlabel)
    plt.ylabel("Counts")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def process_window(folder):

    center = float(folder.split("_")[-1])

    pullx = find_file(folder, "_pullx.xvg")
    pullf = find_file(folder, "_pullf.xvg")

    time, pos = load_xvg_columns(pullx,(0,1))
    _, force = load_xvg_columns(pullf,(0,1))

    plot_line(time,pos,"Time (ps)","COM (nm)",
              f"Position vs Time (COM {center})",
              os.path.join(folder,"position_vs_time.png"))

    plot_line(time,force,"Time (ps)","Pull force (kJ/mol/nm)",
              f"Pull Force vs Time (COM {center})",
              os.path.join(folder,"pull_force_vs_time.png"))

    plot_hist(pos,"COM (nm)",
              f"Position Distribution (COM {center})",
              os.path.join(folder,"position_distribution.png"))

    return center,pos


def plot_pmf(folder):

    bsres = os.path.join(folder,"profile_bootstrap.xvg")
    profile = os.path.join(folder,"profile.xvg")

    plt.figure(figsize=(8,5))

    if os.path.exists(bsres):

        x, pmf, err = np.loadtxt(bsres, comments=['@','#'], unpack=True)
        pmf -= np.min(pmf)  # Shift PMF to zero minimum

        plt.plot(x, pmf, color=cmap(0.7))
        plt.fill_between(x, pmf-err, pmf+err,
                         color=cmap(0.7), alpha=0.3)

    else:
        x, pmf = load_xvg_columns(profile,(0,1))
        plt.plot(x, pmf, color=cmap(0.7))

    plt.xlabel("COM distance Lys167-Tyr259 (nm)")
    plt.ylabel("PMF (kJ/mol)")
    plt.title("Potential of Mean Force")
    plt.tight_layout()
    plt.savefig(os.path.join(folder,"PMF.png"))
    plt.close()


def plot_histogram_overlap(folder):

    histo = os.path.join(folder,"histo.xvg")

    data = np.loadtxt(histo,comments=['@','#']).T
    x = data[0]
    series = data[1:]

    plt.figure(figsize=(8,5))

    colors = cmap(np.linspace(0,1,len(series)))

    for h,c in zip(series,colors):
        plt.plot(x,h,color=c)

    plt.xlabel("COM distance Lys167-Tyr259 (nm)")
    plt.ylabel("Counts")
    plt.title("Histogram Overlap")
    plt.tight_layout()
    plt.savefig(os.path.join(folder,"histogram_overlap.png"))
    plt.close()


def pmf_convergence(folder):

    profile = os.path.join(folder,"profile.xvg")
    x,pmf = load_xvg_columns(profile,(0,1))

    plt.figure(figsize=(8,5))

    n=len(pmf)

    colors=cmap([0.2,0.5,0.8])

    plt.plot(x,pmf,color=colors[2],label="100%")
    plt.plot(x[:int(n*0.75)],pmf[:int(n*0.75)],color=colors[1],label="75%")
    plt.plot(x[:int(n*0.5)],pmf[:int(n*0.5)],color=colors[0],label="50%")

    plt.xlabel("COM distance Lys167-Tyr259 (nm)")
    plt.ylabel("PMF (kJ/mol)")
    plt.title("PMF Convergence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(folder,"PMF_convergence.png"))
    plt.close()


def center_vs_mean(centers,means,folder):

    plt.figure(figsize=(6,6))

    sc = plt.scatter(centers,means,c=centers,cmap="viridis")

    plt.plot(centers,centers,'--',color="black")

    plt.xlabel("Umbrella center (nm)")
    plt.ylabel("Sampled mean (nm)")
    plt.title("Center vs Sampled Mean")

    plt.colorbar(sc,label="Window center (nm)")

    plt.tight_layout()
    plt.savefig(os.path.join(folder,"center_vs_sampled_mean.png"))
    plt.close()


def sampling_width(centers,stds,folder):

    plt.figure(figsize=(8,5))

    plt.plot(centers,stds,'o-',color=cmap(0.6))

    plt.xlabel("Umbrella center (nm)")
    plt.ylabel("Sampling width (nm)")
    plt.title("Window Sampling Width")

    plt.tight_layout()
    plt.savefig(os.path.join(folder,"window_sampling_width.png"))
    plt.close()


def overlap_matrix(all_pos,folder):

    n=len(all_pos)
    overlap=np.zeros((n,n))

    bins=np.linspace(min(np.concatenate(all_pos)),
                     max(np.concatenate(all_pos)),100)

    hists=[np.histogram(p,bins,density=True)[0] for p in all_pos]

    for i in range(n):
        for j in range(n):
            overlap[i,j]=np.sum(np.minimum(hists[i],hists[j]))

    plt.figure(figsize=(6,6))

    im=plt.imshow(overlap,origin="lower",cmap="viridis")

    plt.colorbar(im,label="Overlap")

    plt.xlabel("Window")
    plt.ylabel("Window")
    plt.title("Window Overlap Matrix")

    plt.tight_layout()
    plt.savefig(os.path.join(folder,"window_overlap_matrix.png"))
    plt.close()


def main(folder):

    windows = sorted(
        [f for f in os.listdir(folder) if f.startswith("COM_")],
        key=lambda x: float(x.split("_")[-1])
    )

    centers=[]
    all_pos=[]
    means=[]
    stds=[]

    for w in tqdm(windows):

        path=os.path.join(folder,w)

        if os.path.isdir(path):

            c,pos=process_window(path)

            centers.append(c)
            all_pos.append(pos)
            means.append(np.mean(pos))
            stds.append(np.std(pos))

    plot_pmf(folder)
    plot_histogram_overlap(folder)
    pmf_convergence(folder)

    center_vs_mean(centers,means,folder)
    sampling_width(centers,stds,folder)
    overlap_matrix(all_pos,folder)


if __name__=="__main__":

    parser=argparse.ArgumentParser()
    parser.add_argument("-f","--folder",default=os.getcwd())

    args=parser.parse_args()

    main(args.folder)