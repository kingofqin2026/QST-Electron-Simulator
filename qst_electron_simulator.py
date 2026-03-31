from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

try:
    from scipy.sparse import csr_matrix, diags, kron, identity
    from scipy.sparse.linalg import eigsh
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


PHI = (1.0 + math.sqrt(5.0)) / 2.0


@dataclass
class SimulatorConfig:
    # hard-locked / declared parameters
    kappa: float = 1.0e-2
    g_s: float = 0.1
    D0: float = PHI ** 3
    sigma_init: float = 0.45

    # radial domain
    r_max: float = 20.0
    n_r: int = 512
    r0_shell: float = 0.25
    n_shells: int = 12

    # channel / matrix side
    n_channels: int = 12
    e8_mode: str = "adjacency"  # or "closure"
    sigma_tol: float = 1.0e-4
    sigma_max_iter: int = 20

    # radial BVP lowest-order knobs
    core_weight: float = 1.0
    lambda_D: float = 0.10
    lambda_J: float = 0.12
    lambda_reg: float = 0.06
    reg_eps: float = 1.0e-4
    closure_strength: float = 1.0
    potential_mix: float = 1.0

    # shell / audit
    shell_payload_metric: str = "energy"
    machine_cutoff: float = 1.0e-10


@dataclass
class ChannelInfo:
    index: int
    degree: int
    closeness: float
    vortex_class: str
    chi_a: int
    winding_tag: int


@dataclass
class RadialSolution:
    r: np.ndarray
    u: np.ndarray
    D: np.ndarray
    J_sc: np.ndarray
    rho: np.ndarray
    energy_density: np.ndarray
    eigenvalue: float


@dataclass
class SigmaAudit:
    sigma_in: float
    sigma_out: float
    c_lock: float
    c_n: List[float]
    converged: bool
    iterations: int


@dataclass
class ShellPayload:
    n_star: int
    shell_radius: float
    sigma_out: float
    alpha_coeffs: Dict[str, float]
    delta_psi: Dict[str, float]
    delta_J: Dict[str, float]
    sigma_psi: Dict[str, float]
    sigma_J: Dict[str, float]
    shell_weights: List[float]


@dataclass
class BranchAuditResult:
    branch_name: str
    score: float
    passed_regularity: bool
    passed_finite_core: bool
    passed_monotonic_decay: bool
    notes: List[str] = field(default_factory=list)


@dataclass
class FinalReport:
    config: Dict
    sigma_audit: Dict
    dominant_shell_payload: Dict
    branch_audit: Dict
    closure_status: Dict
    observables: Dict


class E8SelectionGraph:
    """
    Minimal zero-calibration stand-in:
    generates a sparse graph from a deterministic discrete rule.
    This is a working closure, not the final unique E8 basis.
    """

    def __init__(self, n_channels: int, mode: str = "adjacency") -> None:
        self.n_channels = int(n_channels)
        self.mode = mode
        self.adj = self._build_graph()

    def _build_graph(self) -> np.ndarray:
        n = self.n_channels
        A = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                if self.mode == "adjacency":
                    cond = (abs(i - j) == 1) or ((i + j) % 5 == 0)
                else:  # "closure"
                    cond = ((i + j) % 4 == 0) or ((i * j) % 7 == 0 and i != 0 and j != 0)
                if cond:
                    A[i, j] = 1.0
                    A[j, i] = 1.0
        return A

    def classify_channels(self) -> List[ChannelInfo]:
        degrees = self.adj.sum(axis=1)
        # closeness surrogate on unweighted graph
        closeness = np.zeros(self.n_channels, dtype=float)
        for i in range(self.n_channels):
            dists = self._bfs_distances(i)
            finite = dists[np.isfinite(dists)]
            closeness[i] = 0.0 if finite.size <= 1 else (finite.size - 1) / finite.sum()

        order_deg = np.argsort(degrees)
        order_closeness = np.argsort(closeness)

        degree_threshold_low = np.quantile(degrees, 0.33)
        degree_threshold_high = np.quantile(degrees, 0.67)

        infos: List[ChannelInfo] = []
        for a in range(self.n_channels):
            deg = int(degrees[a])
            cl = float(closeness[a])

            if deg <= degree_threshold_low:
                vortex_class = "core-like"
                chi_a = +1
                winding = 0
            elif deg >= degree_threshold_high:
                vortex_class = "flow-arm/chiral"
                chi_a = -1 if (a % 2) else +1
                winding = 1 if (a % 2 == 0) else -1
            else:
                vortex_class = "edge-like"
                chi_a = +1 if (a % 2 == 0) else -1
                winding = 0

            infos.append(
                ChannelInfo(
                    index=a,
                    degree=deg,
                    closeness=cl,
                    vortex_class=vortex_class,
                    chi_a=chi_a,
                    winding_tag=winding,
                )
            )
        return infos

    def _bfs_distances(self, source: int) -> np.ndarray:
        n = self.n_channels
        dist = np.full(n, np.inf)
        dist[source] = 0
        q = [source]
        head = 0
        while head < len(q):
            u = q[head]
            head += 1
            nbrs = np.where(self.adj[u] > 0)[0]
            for v in nbrs:
                if not np.isfinite(dist[v]):
                    dist[v] = dist[u] + 1
                    q.append(int(v))
        return dist


class RadialBVPSolver:
    """
    Minimal stationary radial envelope solver.

    Working closure:
        D(r) = d/dr ln|u(r)| = u'(r)/(u(r)+eps)
        J_sc(r) = kappa * C_core * u(r)^2

    The PDE/BVP is reduced to a Schrödinger-type radial eigenproblem:
        -u'' + U_eff[u] u = lambda u
    with iterative self-consistency.
    """

    def __init__(self, cfg: SimulatorConfig):
        self.cfg = cfg

    def solve(self, sigma: float, core_channel_count: int) -> RadialSolution:
        r = np.linspace(0.0, self.cfg.r_max, self.cfg.n_r)
        dr = r[1] - r[0]

        # finite regular centered initial guess
        u = np.exp(-(r / 2.0) ** 2) * (1.0 + 0.15 * np.exp(-(r / 0.8) ** 2))
        u /= math.sqrt(np.trapz(u ** 2, r))

        c_core = max(1, core_channel_count)

        for _ in range(14):
            D = self._compute_D(u, r)
            J = self._compute_J(u, c_core)
            U = self._effective_potential(r, D, J, sigma)

            eigval, eigvec = self._ground_state(U, dr)
            u_new = np.abs(eigvec)
            u_new /= math.sqrt(max(np.trapz(u_new ** 2, r), 1e-18))

            # under-relaxation for stability
            u = 0.65 * u + 0.35 * u_new
            u /= math.sqrt(max(np.trapz(u ** 2, r), 1e-18))

        D = self._compute_D(u, r)
        J = self._compute_J(u, c_core)
        rho = u ** 2
        energy_density = 0.5 * np.gradient(u, r) ** 2 + 0.5 * U * rho

        return RadialSolution(
            r=r,
            u=u,
            D=D,
            J_sc=J,
            rho=rho,
            energy_density=energy_density,
            eigenvalue=float(eigval),
        )

    def _compute_D(self, u: np.ndarray, r: np.ndarray) -> np.ndarray:
        eps = self.cfg.reg_eps
        uprime = np.gradient(u, r)
        D = self.cfg.D0 + self.cfg.closure_strength * (uprime / (u + eps))
        # clip to keep the toy solver stable; this is numerical regularization, not fitting
        return np.clip(D, self.cfg.D0 - 4.0, self.cfg.D0 + 4.0)

    def _compute_J(self, u: np.ndarray, core_channel_count: int) -> np.ndarray:
        c_core = float(core_channel_count)
        return self.cfg.kappa * c_core * (u ** 2)

    def _effective_potential(self, r: np.ndarray, D: np.ndarray, J: np.ndarray, sigma: float) -> np.ndarray:
        reg = self.cfg.lambda_reg / (r ** 2 + self.cfg.reg_eps)
        UD = self.cfg.lambda_D * (D - self.cfg.D0) ** 2
        UJ = self.cfg.lambda_J * (J ** 2)
        sigma_term = self.cfg.potential_mix * sigma ** 2 * np.abs(D - self.cfg.D0)
        return reg + UD + UJ + sigma_term

    def _ground_state(self, U: np.ndarray, dr: float) -> Tuple[float, np.ndarray]:
        n = U.size
        main = (2.0 / dr ** 2) + U
        off = -np.ones(n - 1) / dr ** 2

        if SCIPY_AVAILABLE:
            H = diags([off, main, off], [-1, 0, 1], shape=(n, n), format="csr")
            vals, vecs = eigsh(H, k=1, which="SA")
            vec = np.asarray(vecs[:, 0]).reshape(-1)
            return float(vals[0]), vec

        H = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
        vals, vecs = np.linalg.eigh(H)
        return float(vals[0]), np.asarray(vecs[:, 0]).reshape(-1)


class DSIShellProjector:
    def __init__(self, cfg: SimulatorConfig):
        self.cfg = cfg

    def shell_radii(self) -> List[float]:
        return [self.cfg.r0_shell * (PHI ** (2 * n)) for n in range(self.cfg.n_shells)]

    def project(self, sol: RadialSolution) -> Dict[int, Dict[str, float]]:
        r = sol.r
        radii = self.shell_radii()
        payload: Dict[int, Dict[str, float]] = {}

        prev_edge = 0.0
        for n, edge in enumerate(radii):
            if n == len(radii) - 1:
                mask = r >= prev_edge
            else:
                mask = (r >= prev_edge) & (r < edge)
            if not np.any(mask):
                payload[n] = {"rho": 0.0, "J_sc": 0.0, "energy": 0.0, "D_var": 0.0}
                prev_edge = edge
                continue

            rr = r[mask]
            rho = float(np.trapz(sol.rho[mask], rr))
            J_sc = float(np.trapz(sol.J_sc[mask], rr))
            energy = float(np.trapz(sol.energy_density[mask], rr))
            D_var = float(np.var(sol.D[mask]))

            payload[n] = {"rho": rho, "J_sc": J_sc, "energy": energy, "D_var": D_var}
            prev_edge = edge
        return payload


class QSTMatrixBuilder:
    def __init__(self, cfg: SimulatorConfig, channels: List[ChannelInfo], graph: E8SelectionGraph):
        self.cfg = cfg
        self.channels = channels
        self.graph = graph

    def build_dense(self, sigma: float, shell_weights: np.ndarray) -> np.ndarray:
        """
        Build a minimal dense block-tridiagonal Hamiltonian.
        Suitable for small toy sizes in the demo.
        """
        N = self.cfg.n_shells
        d = self.cfg.n_channels
        M = N * d
        H = np.zeros((M, M), dtype=float)

        C = self.graph.adj.copy()
        for n in range(N):
            E_n = self.cfg.kappa * sigma ** 2 * (PHI ** (-2 * n))
            idx = slice(n * d, (n + 1) * d)
            H[idx, idx] = np.eye(d) * E_n

            if n < N - 1:
                J_n = self.cfg.g_s * sigma * (PHI ** (-(2 * n + 1)))
                idx2 = slice((n + 1) * d, (n + 2) * d)
                H[idx, idx2] = J_n * C
                H[idx2, idx] = (J_n * C).T
        return H

    def lowest_state(self, sigma: float, shell_weights: np.ndarray) -> Tuple[float, np.ndarray]:
        H = self.build_dense(sigma, shell_weights)
        vals, vecs = np.linalg.eigh(H)
        return float(vals[0]), np.asarray(vecs[:, 0]).reshape(-1)

    def sigma_self_consistent(self, sigma_in: float, shell_weights: np.ndarray) -> SigmaAudit:
        sigma = float(sigma_in)
        c_n_last: List[float] = []
        converged = False

        for it in range(1, self.cfg.sigma_max_iter + 1):
            _, psi = self.lowest_state(sigma, shell_weights)
            psi_by_shell = self._slice_by_shell(psi)

            c_n = []
            for n in range(len(psi_by_shell) - 1):
                a = psi_by_shell[n]
                b = psi_by_shell[n + 1]
                denom = np.linalg.norm(a) * np.linalg.norm(b) + 1e-12
                cn = float(abs(np.vdot(a, b)) / denom)
                c_n.append(cn)

            c_lock = float(np.mean(c_n)) if c_n else 0.0
            sigma_new = math.tanh(PHI * c_lock)

            c_n_last = c_n
            if abs(sigma_new - sigma) < self.cfg.sigma_tol:
                sigma = sigma_new
                converged = True
                return SigmaAudit(
                    sigma_in=sigma_in,
                    sigma_out=sigma,
                    c_lock=c_lock,
                    c_n=c_n_last,
                    converged=converged,
                    iterations=it,
                )

            sigma = 0.6 * sigma + 0.4 * sigma_new

        c_lock = float(np.mean(c_n_last)) if c_n_last else 0.0
        return SigmaAudit(
            sigma_in=sigma_in,
            sigma_out=sigma,
            c_lock=c_lock,
            c_n=c_n_last,
            converged=converged,
            iterations=self.cfg.sigma_max_iter,
        )

    def _slice_by_shell(self, psi: np.ndarray) -> List[np.ndarray]:
        d = self.cfg.n_channels
        N = self.cfg.n_shells
        return [psi[n * d:(n + 1) * d] for n in range(N)]


class DominantShellAuditor:
    def __init__(self, cfg: SimulatorConfig, channels: List[ChannelInfo]):
        self.cfg = cfg
        self.channels = channels

    def branch_audit(self, sol: RadialSolution) -> BranchAuditResult:
        u = sol.u
        r = sol.r

        passed_regularity = bool(abs(u[0] - np.max(u[:5])) < 0.2)
        passed_finite_core = bool(np.isfinite(sol.D[0]) and np.isfinite(sol.J_sc[0]))
        tail = u[int(0.8 * len(u)):]
        passed_monotonic_decay = bool(np.mean(np.diff(tail)) <= 1e-3)

        score = float(passed_regularity + passed_finite_core + passed_monotonic_decay) / 3.0
        notes = []
        if not passed_regularity:
            notes.append("core regularity weak")
        if not passed_finite_core:
            notes.append("non-finite core closure")
        if not passed_monotonic_decay:
            notes.append("tail decay not monotone enough")

        return BranchAuditResult(
            branch_name="electron_filled_core_v0",
            score=score,
            passed_regularity=passed_regularity,
            passed_finite_core=passed_finite_core,
            passed_monotonic_decay=passed_monotonic_decay,
            notes=notes,
        )

    def dominant_shell_payload(
        self,
        shell_data: Dict[int, Dict[str, float]],
        sigma_audit: SigmaAudit,
        graph_channels: List[ChannelInfo],
    ) -> ShellPayload:
        metric = self.cfg.shell_payload_metric
        shell_scores = [shell_data[n][metric] for n in range(self.cfg.n_shells)]
        n_star = int(np.argmax(shell_scores))
        shell_radius = self.cfg.r0_shell * (PHI ** (2 * n_star))

        # channel coefficients: discrete allocation by class
        class_groups: Dict[str, List[int]] = {"core-like": [], "edge-like": [], "flow-arm/chiral": []}
        for ch in graph_channels:
            class_groups[ch.vortex_class].append(ch.index)

        alpha_coeffs = {
            cls: float(len(idxs)) / max(1, len(graph_channels))
            for cls, idxs in class_groups.items()
        }

        base = shell_data[n_star]
        delta_psi = {
            cls: alpha_coeffs[cls] * base["rho"] for cls in class_groups
        }
        delta_J = {
            cls: alpha_coeffs[cls] * base["J_sc"] for cls in class_groups
        }
        sigma_psi = {
            cls: sigma_audit.sigma_out * alpha_coeffs[cls] for cls in class_groups
        }
        sigma_J = {
            cls: sigma_audit.sigma_out * alpha_coeffs[cls] for cls in class_groups
        }

        return ShellPayload(
            n_star=n_star,
            shell_radius=shell_radius,
            sigma_out=sigma_audit.sigma_out,
            alpha_coeffs=alpha_coeffs,
            delta_psi=delta_psi,
            delta_J=delta_J,
            sigma_psi=sigma_psi,
            sigma_J=sigma_J,
            shell_weights=shell_scores,
        )


class EffectiveCurrentEmitter:
    def __init__(self, cfg: SimulatorConfig):
        self.cfg = cfg

    def emit_scalar_observable(self, sol: RadialSolution, sigma_out: float) -> Dict[str, float]:
        """
        Reduced scalar version of:
            J_eff,e^mu = Pi_M Pi_H [ (1 + sigma^2 D_e) \bar Psi_e gamma^mu Psi_e ]
        Since the simulator is lowest-order radial, use rho = |u|^2 as the current-density carrier.
        """
        j_eff = (1.0 + sigma_out ** 2 * sol.D) * sol.rho

        r = sol.r
        total = float(np.trapz(j_eff, r))
        peak = float(np.max(j_eff))
        mean = float(np.mean(j_eff))

        # qst energy-flow surrogate
        gradD = np.gradient(sol.D, r)
        s_qst = self.cfg.kappa * sigma_out ** 2 * np.abs(gradD * sol.J_sc)
        total_s_qst = float(np.trapz(s_qst, r))

        return {
            "J_eff_total_scalar": total,
            "J_eff_peak_scalar": peak,
            "J_eff_mean_scalar": mean,
            "S_QST_total_surrogate": total_s_qst,
            "residue_normalization_open": 1.0,  # placeholder only
        }


class QSTElectronIntegratedAuditSimulator:
    def __init__(self, cfg: Optional[SimulatorConfig] = None):
        self.cfg = cfg or SimulatorConfig()

    def run(self) -> FinalReport:
        graph = E8SelectionGraph(self.cfg.n_channels, self.cfg.e8_mode)
        channels = graph.classify_channels()

        core_channel_count = sum(1 for c in channels if c.vortex_class == "core-like")

        radial_solver = RadialBVPSolver(self.cfg)
        sol = radial_solver.solve(sigma=self.cfg.sigma_init, core_channel_count=core_channel_count)

        projector = DSIShellProjector(self.cfg)
        shell_data = projector.project(sol)
        shell_weights = np.array([shell_data[n]["energy"] for n in range(self.cfg.n_shells)], dtype=float)

        matrix = QSTMatrixBuilder(self.cfg, channels, graph)
        sigma_audit = matrix.sigma_self_consistent(self.cfg.sigma_init, shell_weights)

        # rerun the radial side once with emitted sigma
        sol = radial_solver.solve(sigma=sigma_audit.sigma_out, core_channel_count=core_channel_count)
        shell_data = projector.project(sol)

        auditor = DominantShellAuditor(self.cfg, channels)
        branch_audit = auditor.branch_audit(sol)
        payload = auditor.dominant_shell_payload(shell_data, sigma_audit, channels)

        emitter = EffectiveCurrentEmitter(self.cfg)
        observables = emitter.emit_scalar_observable(sol, sigma_audit.sigma_out)

        closure_status = {
            "formal_microscopic_source_chain": "PASS (upstream assumption)",
            "branch_to_field_lowest_order_working_closure": "PASS",
            "effective_current_emission": "PASS",
            "actual_selected_branch_emitted_dominant_shell_payload": "EVALUATED IN TOY/V0 FORM ONLY",
            "half_spin_inevitability": "OPEN",
            "unique_practical_E8_basis_chi_a": "OPEN",
            "full_geometry_closure_for_selected_branch": "OPEN",
            "one_loop_normalization_N(s_*)": "OPEN",
            "oscillatory_residue_amplitude_A": "OPEN",
        }

        report = FinalReport(
            config=asdict(self.cfg),
            sigma_audit=asdict(sigma_audit),
            dominant_shell_payload=asdict(payload),
            branch_audit=asdict(branch_audit),
            closure_status=closure_status,
            observables=observables,
        )
        return report

    def save_report(self, report: FinalReport, path: str | Path) -> None:
        p = Path(path)
        p.write_text(json.dumps(asdict(report), indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    sim = QSTElectronIntegratedAuditSimulator()
    report = sim.run()
    out = Path(__file__).with_name("demo_report.json")
    sim.save_report(report, out)
    print(f"Wrote report to: {out}")
    print(json.dumps(report.closure_status, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
