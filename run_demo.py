from pathlib import Path
import json

from qst_electron_simulator import QSTElectronIntegratedAuditSimulator, SimulatorConfig

cfg = SimulatorConfig(
    sigma_init=0.40,
    n_shells=10,
    n_channels=10,
    n_r=384,
    r_max=16.0,
)

sim = QSTElectronIntegratedAuditSimulator(cfg)
report = sim.run()

out = Path(__file__).with_name("demo_report.json")
sim.save_report(report, out)

print("Demo completed.")
print(f"Report saved to: {out}")
print(json.dumps(report.dominant_shell_payload, indent=2, ensure_ascii=False))
