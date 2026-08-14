def save_basic_history(frame, path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(); ax.plot(frame["t_s"]/60, frame["rho"]); ax.set(xlabel="time (min)", ylabel="relative density"); fig.tight_layout(); fig.savefig(path); plt.close(fig)
