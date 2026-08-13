#!/usr/bin/env python3
"""Interacting local pore regions with conservative, schedule-local kinetics.

The local-region parameters are deliberately explicit: every decoder block enters
either initialization, a conservative pore transfer, closed-pore densification,
or the migration-only resistance.  Density is always reconstructed from the
four local pore-volume stores.
"""
from dataclasses import dataclass
import numpy as np
from massive_latent_topology_models import R, sigmoid


@dataclass
class NetworkState:
    rho: np.ndarray
    G: np.ndarray
    phi_GBseg: np.ndarray
    phi_TJ: np.ndarray
    phi_iso: np.ndarray
    phi_closed: np.ndarray
    N_GBseg: np.ndarray
    N_TJ: np.ndarray
    N_iso: np.ndarray
    N_closed: np.ndarray
    connected_removable_fraction: np.ndarray
    damaged_connected_fraction: np.ndarray
    sweep_coalescence_seed: np.ndarray
    large_attached_fraction: np.ndarray
    large_TJ_fraction: np.ndarray
    isolated_fraction: np.ndarray
    closed_fraction: np.ndarray
    X_J: np.ndarray
    C_TJ: np.ndarray
    C_GBseg: np.ndarray
    f_clean_GB: np.ndarray
    residual_stress: np.ndarray
    PR_damage_memory: np.ndarray
    sweep_memory: np.ndarray
    closed_accommodation: np.ndarray
    migration_factor: np.ndarray
    densification_eligibility: np.ndarray
    weights: np.ndarray


def clone_state(state):
    return NetworkState(**{
        name: (getattr(state, name).copy() if isinstance(getattr(state, name), np.ndarray)
               else getattr(state, name))
        for name in state.__dataclass_fields__
    })


def _thermal(T_C, Q, reference_K=1473.15):
    T = T_C + 273.15
    exponent = -Q / R * (1.0 / T - 1.0 / reference_K)
    return np.exp(np.clip(exponent, -60.0, 60.0))


def network_adjacency(n, p, seed=1):
    """Deterministic ring-plus-shortcuts graph controlled by degree/clustering."""
    degree = int(np.clip(p.get("degree", 2), 1, max(n - 1, 1)))
    cluster = float(np.clip(p.get("cluster", 0.0), 0.0, 1.0))
    adjacency = np.zeros((n, n), float)
    for offset in range(1, degree + 1):
        for i in range(n):
            local_j = (i + offset) % n
            long_j = (i + 1 + ((i * 7 + offset * 11 + seed) % max(n - 1, 1))) % n
            j = local_j if ((i + offset) % 100) / 99.0 < cluster else long_j
            if j != i:
                adjacency[i, j] = adjacency[j, i] = 1.0
    return adjacency


def initial(n=8, rho0=.70, G0=100., seed=1, p=None):
    p = {} if p is None else p
    rng = np.random.default_rng(seed)
    w = rng.lognormal(0.0, p.get("weight_sigma", .25), n)
    w /= w.sum()
    pore = (1.0 - rho0) * np.clip(
        rng.normal(1.0, p.get("rho_sigma", .12), n), .25, 2.0
    )
    pore *= (1.0 - rho0) / (w @ pore)
    z = np.linspace(-1.0, 1.0, n)
    local_connected = np.clip(
        p.get("connected_init", .6) * (1.0 + p.get("cluster", 0.0) * .5 * z),
        .02,
        .98,
    )
    closed = np.clip(p.get("closed_init", 0.0), 0.0, .4) * pore
    available = pore - closed
    gb = local_connected * available
    tj = .25 * (available - gb)
    iso = available - gb - tj
    zeros = lambda: np.zeros(n)
    ones = lambda: np.ones(n)
    G = np.maximum(G0 * rng.lognormal(0.0, p.get("G_sigma", .05), n), 1.0)
    state = NetworkState(
        1.0 - pore, G, gb, tj, iso, closed,
        gb.copy(), tj.copy(), iso.copy(), closed.copy(),
        gb / np.maximum(pore, 1e-30), .05 * ones(), zeros(), zeros(), zeros(),
        iso / np.maximum(pore, 1e-30), closed / np.maximum(pore, 1e-30),
        .1 * ones(), tj / np.maximum(pore, 1e-30), gb / np.maximum(pore, 1e-30),
        1.0 - (gb + tj) / np.maximum(pore, 1e-30), zeros(), zeros(), zeros(),
        p.get("closed_capacity", 1.0) * ones(), ones(), ones(), w,
    )
    return state


def global_observables(s):
    pore = s.phi_GBseg + s.phi_TJ + s.phi_iso + s.phi_closed
    connected_mean = np.average(s.connected_removable_fraction, weights=s.weights)
    return dict(
        rho_global=1.0 - float(s.weights @ pore),
        G_mean=float(s.weights @ s.G),
        topology_variance=float(np.average(
            (s.connected_removable_fraction - connected_mean) ** 2,
            weights=s.weights,
        )),
        closed_fraction=float(s.weights @ s.phi_closed / max(s.weights @ pore, 1e-30)),
    )


def local_fluxes(s, T_C, p):
    """Instantaneous local laws using only temperature and the current state."""
    total_pore = np.maximum(1.0 - s.rho, 1e-30)
    base_activity = sigmoid((T_C - 1180.0) / 70.0)
    stress_inhibition = np.exp(-np.clip(p.get("stress_nucleation", 0.0) * s.residual_stress, 0, 60))
    activity = np.clip(base_activity * stress_inhibition, 0.0, 1.0)
    open_flux = (
        p["k_open"] * _thermal(T_C, p["Q_density"]) * s.phi_GBseg
        * np.clip(s.connected_removable_fraction, 0, 1)
        * np.clip(s.densification_eligibility, 0, 1)
    )
    radius_proxy = np.clip(
        s.phi_closed / np.maximum(s.N_closed, 1e-30), 1e-6, 1e6
    )
    radius_factor = radius_proxy ** (-p.get("closed_radius_exp", 0.0) / 3.0)
    closed_flux = (
        p["k_closed"] * _thermal(T_C, p["Q_closed"]) * s.phi_closed
        * np.clip(s.closed_accommodation, 0, p.get("closed_capacity", 1.0))
        * max(1.0 - p["gas_ratio"], 0.0) * radius_factor
    )
    low_activity = sigmoid(
        (p.get("activity_mid", .2) - activity) / max(p.get("activity_width", .1), 1e-6)
    )
    pr = p["k_PR"] * _thermal(T_C, p["Q_PR"]) * low_activity * s.phi_GBseg
    detach = p.get("detachment", 0.0) * _thermal(T_C, p["Q_PR"]) * s.phi_GBseg * np.clip(s.large_attached_fraction, 0, 1)
    recapture = p.get("recapture", 0.0) * _thermal(T_C, p["Q_PR"]) * s.phi_iso * np.clip(s.connected_removable_fraction, 0, 1)
    close_source = s.phi_iso + .5 * s.phi_TJ
    close_transition = p.get("closed_transition", 0.0) * _thermal(T_C, p["Q_closed"]) * close_source * np.clip(s.rho, 0, 1) ** 3

    C_pore = np.clip(s.C_TJ, 0, 1)
    C_constraint = np.clip(s.X_J * (1.0 - p.get("pore_relax", 0.0) * C_pore), 0, 1)
    thermal_events = _thermal(T_C, min(p["Q_growth"], p.get("Q_PR", p["Q_growth"])))
    lambda_eff = p.get("lambda_TJ", 1.0) * thermal_events * np.maximum(C_constraint, 1e-8)
    K_eff = p.get("K_TJ", 1.0) * np.maximum(s.G / 100.0, 1e-6) ** p.get("q_TJ", 0)
    P_comp = sigmoid((lambda_eff - K_eff) / np.sqrt(lambda_eff + K_eff + 1e-12))
    pore_drag = p.get("pore_drag_fraction", 0.0) * C_pore
    drag = (
        p["attached_drag"] * (s.large_attached_fraction + s.large_TJ_fraction + pore_drag)
        + p["junction_drag"] * s.X_J
        + p.get("stress_migration", p.get("stress_drag", 1.0)) * s.residual_stress
    )
    migration_factor = np.clip(P_comp / (1.0 + np.maximum(drag, 0.0)), 0.0, 1.0)
    growth = p["k_growth"] * _thermal(T_C, p["Q_growth"]) * migration_factor / np.maximum(s.G, 1.0)
    migration_rate = growth / np.maximum(s.G, 1.0)
    sweep = (
        p.get("k_sweep_damaged", p.get("k_sweep", 20.0))
        * migration_rate ** p["sweep_exp"] * np.clip(s.damaged_connected_fraction, 0, 1)
        + p.get("k_sweep_connected", 0.0)
        * migration_rate ** p["sweep_exp"] * s.phi_GBseg
    )
    return dict(
        rho_dot_open=np.maximum(open_flux, 0),
        rho_dot_closed=np.maximum(closed_flux, 0),
        PR_damage=np.maximum(pr, 0),
        sweep=np.maximum(sweep, 0),
        detachment=np.maximum(detach, 0),
        recapture=np.maximum(recapture, 0),
        closed_transition=np.maximum(close_transition, 0),
        G_dot=np.maximum(growth, 0),
        migration_factor=migration_factor,
        activity=activity,
        Lambda_TJ=lambda_eff,
        K_TJ=K_eff,
        P_comp_TJ=P_comp,
    )


def _bounded_amount(rate, available, dt):
    return np.minimum(np.maximum(rate, 0.0) * dt, np.maximum(available, 0.0))


def advance(s, T_C, p, dt, adj=None):
    f = local_fluxes(s, T_C, p)
    pore_before = s.phi_GBseg + s.phi_TJ + s.phi_iso + s.phi_closed

    pr = _bounded_amount(f["PR_damage"], s.phi_GBseg, dt)
    gb_after_pr = s.phi_GBseg - pr
    sweep = _bounded_amount(f["sweep"], gb_after_pr, dt)
    gb_after_sweep = gb_after_pr - sweep
    detach = _bounded_amount(f["detachment"], gb_after_sweep, dt)
    recapture = _bounded_amount(f["recapture"], s.phi_iso, dt)

    # PR redistribution is conservative and uses every declared partition.
    s.phi_GBseg = gb_after_sweep - detach + recapture
    s.phi_GBseg += (p["PR_damaged"] + p["PR_large"]) * pr + .70 * sweep
    s.phi_TJ += p["PR_TJ"] * pr + .20 * sweep
    s.phi_iso += p["PR_iso"] * pr + detach + .10 * sweep - recapture
    s.phi_closed += p["PR_closed"] * pr

    # Isolated/TJ-to-closed transition is conservative and cannot densify.
    transition = _bounded_amount(f["closed_transition"], s.phi_iso + .5 * s.phi_TJ, dt)
    from_iso = np.minimum(transition, s.phi_iso)
    from_tj = np.minimum(transition - from_iso, s.phi_TJ)
    s.phi_iso -= from_iso
    s.phi_TJ -= from_tj
    s.phi_closed += from_iso + from_tj

    open_loss = _bounded_amount(f["rho_dot_open"], s.phi_GBseg, dt)
    closed_loss = _bounded_amount(f["rho_dot_closed"], s.phi_closed, dt)
    s.phi_GBseg -= open_loss
    s.phi_closed -= closed_loss

    coalescence_strength = np.clip(
        p.get("number_loss", 0.0)
        * (1.0 + p.get("coalescence_exp", 1.0) * sweep / np.maximum(pore_before, 1e-30)),
        0.0, 10.0,
    )
    s.N_GBseg = np.maximum(s.N_GBseg - open_loss - coalescence_strength * sweep - detach + recapture, 0)
    s.N_TJ = np.maximum(s.N_TJ + p["PR_TJ"] * pr + .10 * sweep - from_tj, 0)
    s.N_iso = np.maximum(s.N_iso + p["PR_iso"] * pr + detach + .05 * sweep - recapture - from_iso, 0)
    s.N_closed = np.maximum(s.N_closed + p["PR_closed"] * pr + from_iso + from_tj, 0)

    total_pore = np.maximum(s.phi_GBseg + s.phi_TJ + s.phi_iso + s.phi_closed, 1e-30)
    memory_relax = np.exp(-dt / max(p.get("PR_tau", 1e30), 1e-30))
    s.PR_damage_memory = np.clip(s.PR_damage_memory * memory_relax + pr / total_pore, 0, 1)
    s.damaged_connected_fraction = np.clip(
        s.damaged_connected_fraction * memory_relax + pr / total_pore - sweep / total_pore, 0, 1
    )
    s.sweep_memory = np.clip(s.sweep_memory + sweep / total_pore, 0, 1)
    radius_gain = np.clip(p.get("coalescence_exp", 1.0) * sweep / total_pore, 0, 1)
    s.large_attached_fraction = np.clip(s.large_attached_fraction + p["PR_large"] * pr / total_pore + radius_gain, 0, 1)
    s.large_TJ_fraction = np.clip(s.large_TJ_fraction + p["PR_TJ"] * pr / total_pore + .5 * radius_gain, 0, 1)

    xj_relax = np.exp(-dt / max(p.get("XJ_tau", 1e30), 1e-30))
    xj_source = p.get("XJ_prod", 0.0) * (p["PR_TJ"] * pr + .2 * sweep) / total_pore
    s.X_J = np.clip(s.X_J * xj_relax + xj_source, 0, 1)
    stress_relax = np.exp(-dt / max(p.get("stress_tau", 1e30), 1e-30))
    stress_source = (
        p.get("stress_PR", 0.0) * pr / total_pore
        + p.get("stress_shear", 0.0) * f["G_dot"] * dt / np.maximum(s.G, 1.0)
    )
    s.residual_stress = np.clip(s.residual_stress * stress_relax + stress_source, 0, 1e3)

    cap = p.get("closed_capacity", 1.0)
    cap_relax = 1.0 - np.exp(-dt / max(p.get("capacity_tau", 1e30), 1e-30))
    s.closed_accommodation += cap_relax * (cap - s.closed_accommodation)
    s.closed_accommodation = np.clip(s.closed_accommodation - closed_loss / total_pore, 0, cap)

    # Exact integration of dG/dt = coefficient/G avoids Euler runaway.
    s.G = np.sqrt(np.maximum(s.G * s.G + 2.0 * s.G * f["G_dot"] * dt, 1.0))
    s.rho = 1.0 - (s.phi_GBseg + s.phi_TJ + s.phi_iso + s.phi_closed)
    s.connected_removable_fraction = np.clip(s.phi_GBseg / total_pore, 0, 1)
    s.isolated_fraction = np.clip(s.phi_iso / total_pore, 0, 1)
    s.closed_fraction = np.clip(s.phi_closed / total_pore, 0, 1)
    s.C_TJ = np.clip(s.phi_TJ / total_pore, 0, 1)
    s.C_GBseg = np.clip(s.phi_GBseg / total_pore, 0, 1)
    s.f_clean_GB = np.clip(1.0 - s.C_TJ - s.C_GBseg, 0, 1)
    s.migration_factor = f["migration_factor"]

    if adj is not None:
        neighbor = adj @ s.connected_removable_fraction / np.maximum(adj.sum(axis=1), 1.0)
        alpha = 1.0 - np.exp(-max(p.get("exchange_rate", 0.0), 0.0) * dt)
        s.connected_removable_fraction = np.clip(
            s.connected_removable_fraction + alpha * (neighbor - s.connected_removable_fraction), 0, 1
        )
    return f


def defaults():
    return dict(
        k_open=1.5e-5, Q_density=475e3, k_closed=2e-6, Q_closed=475e3,
        k_PR=1e-5, Q_PR=250e3, k_sweep=20., k_sweep_damaged=20.,
        k_sweep_connected=0., sweep_exp=1., coalescence_exp=1., k_growth=9e3,
        Q_growth=500e3, activity_mid=1180., activity_width=70., gas_ratio=.25,
        attached_drag=30., junction_drag=10., stress_drag=1., stress_migration=1.,
        stress_nucleation=0., stress_PR=0., stress_shear=0., stress_tau=1e30,
        number_loss=2., exchange_rate=1e-7, detachment=0., recapture=0.,
        closed_transition=0., closed_capacity=1., capacity_tau=1e30,
        closed_radius_exp=0., PR_damaged=.4, PR_large=.2, PR_TJ=.15,
        PR_iso=.15, PR_closed=.1, PR_tau=1e30, XJ_prod=0., XJ_tau=1e30,
        lambda_TJ=1e3, K_TJ=1., q_TJ=0, pore_relax=0., pore_drag_fraction=0.,
        degree=2, cluster=0., weight_sigma=.25, rho_sigma=.12, G_sigma=.05,
        connected_init=.6, closed_init=0.,
    )


LOCAL_FUNCTIONS = (local_fluxes,)
