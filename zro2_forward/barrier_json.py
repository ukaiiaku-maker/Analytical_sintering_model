from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np


class BarrierInputError(RuntimeError):
    pass


def _find(d: dict, names: tuple[str, ...]):
    for name in names:
        if name in d:
            return d[name]
    raise BarrierInputError(f"barrier slice lacks one of {names}")


@dataclass
class BarrierModel:
    temperatures_K: np.ndarray
    G0_J: np.ndarray
    Gfloor_J: np.ndarray
    a: np.ndarray
    sigmahat_Pa: np.ndarray
    n: np.ndarray
    schema: str

    @classmethod
    def load(cls, path: str | Path) -> "BarrierModel":
        path = Path(path)
        if not path.is_file():
            raise BarrierInputError(f"required fitted barrier JSON is missing: {path}")
        raw = json.loads(path.read_text())
        schema = str(raw.get("schema", raw.get("model_schema", "")))
        if schema != "bicrystal_surface_triple_line_arrhenius_EXP_floor_v1":
            raise BarrierInputError(f"unsupported barrier schema: {schema!r}")
        slices = raw.get("temperature_slices", raw.get("slices"))
        if not isinstance(slices, list) or len(slices) < 2:
            raise BarrierInputError("barrier JSON requires at least two fitted temperature slices")
        eV_J = 1.602176634e-19
        vals = []
        for s in slices:
            T = float(_find(s, ("T_K", "temperature_K", "temperature")))
            def energy(j_names, ev_names):
                try:
                    return float(_find(s, j_names))
                except BarrierInputError:
                    return float(_find(s, ev_names)) * eV_J
            vals.append((T, energy(("G0_J",), ("G0_eV",)),
                         energy(("Gfloor_J", "G_floor_J"), ("Gfloor_eV", "G_floor_eV")),
                         float(_find(s, ("a", "a_fit"))),
                         float(_find(s, ("sigmahat_Pa", "sigma_hat_Pa"))),
                         float(_find(s, ("n", "stress_exponent")))))
        vals.sort()
        a = np.asarray(vals, float)
        return cls(a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5], schema)

    def _p(self, values: np.ndarray, T_K: float) -> float:
        # Shape-preserving interpolation is preferred; linear is deliberately not
        # substituted. scipy is optional, so fail clearly when prediction is asked.
        try:
            from scipy.interpolate import PchipInterpolator
        except ImportError as exc:
            raise BarrierInputError("scipy is required for PCHIP barrier interpolation") from exc
        return float(PchipInterpolator(self.temperatures_K, values, extrapolate=False)(T_K))

    def Gstar(self, sigma_Pa: float, T_K: float) -> float:
        G0, floor = self._p(self.G0_J, T_K), self._p(self.Gfloor_J, T_K)
        a, sh, n = self._p(self.a, T_K), self._p(self.sigmahat_Pa, T_K), self._p(self.n, T_K)
        return floor + (G0 - floor) * np.exp(-a * (max(sigma_Pa, 0.0) / sh) ** n)

    def diagnostics(self, sigma_Pa: float, T_K: float) -> dict[str, float]:
        dT = max(0.1, T_K * 1e-4); ds = max(1e3, sigma_Pa * 1e-4)
        return {"Gstar_J": self.Gstar(sigma_Pa, T_K),
                "Sstar_J_K": -(self.Gstar(sigma_Pa, T_K+dT)-self.Gstar(sigma_Pa, T_K-dT))/(2*dT),
                "Vstar_m3": -(self.Gstar(sigma_Pa+ds, T_K)-self.Gstar(max(0., sigma_Pa-ds), T_K))/(2*ds)}
