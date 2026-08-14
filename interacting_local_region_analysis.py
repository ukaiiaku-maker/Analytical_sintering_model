#!/usr/bin/env python3
import numpy as np
def distribution_metrics(state):
 g=np.sort(state.G);return dict(G50=float(np.median(g)),G90=float(np.quantile(g,.9)),connected_mean=float(np.average(state.connected_removable_fraction,weights=state.weights)),topology_variance=float(np.var(state.connected_removable_fraction)))
