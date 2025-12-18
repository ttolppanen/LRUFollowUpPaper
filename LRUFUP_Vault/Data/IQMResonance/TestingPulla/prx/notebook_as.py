# %%
import os
from dotenv import load_dotenv
load_dotenv()
from iqm.pulla.pulla import Pulla
from iqm.qiskit_iqm import IQMProvider

iqm_server_url = "https://resonance.meetiqm.com"
quantum_computer = "garnet"
os.environ["IQM_TOKEN"] = os.environ.get("RESONANCE_API_TOKEN")

p = Pulla(iqm_server_url, quantum_computer=quantum_computer)
provider = IQMProvider(iqm_server_url, quantum_computer=quantum_computer)
backend = provider.get_backend(use_metrics=False)
compiler = p.get_standard_compiler()

# %%
qubits = compiler.builder.chip_topology.qubits_sorted
q1 = qubits[0]
print(f"Using qubit {q1}")

from iqm.pulse import Circuit
from iqm.pulse import CircuitOperation as Op
import numpy as np


instructions = []
instructions.append(Op("prx", (q1, ), args={"angle": np.pi, "phase": 0.0})) # |0> -> |1>
instructions.append(Op("prx_12", (q1, ), args={"angle": np.pi, "phase": 0.0})) # |1> -> |2>
instructions.append(Op("measure", (q1, ), args={"key": "M"})) # measure

circuit = Circuit("2 preparation", instructions)

compiler.amend_calibration_for_gate_implementation("measure", "constant", (q1, ), {"acquisition_type": "complex"}) # Enable complex readout to observe |2> state.

playlist, context = compiler.compile([circuit])
settings, context = compiler.build_settings(context, shots=100)

# settings.options.averaging_bins = 1

from iqm.pulse.playlist.visualisation.base import inspect_playlist
from IPython.core.display import HTML

HTML(inspect_playlist(playlist, [0]))

# %%
job = p.submit_playlist(playlist, settings, context=context)
job.wait_for_completion()

# %%
job.result()[0]["M"]

# %%
# p.get_channel_properties()[0]['TC-16-20__flux.awg']
# for k, v in compiler.get_calibration().items():
#     if k.startswith(""):
#         print(f"{k}: {v}")
# print({k: v for k, v in compiler.get_calibration().items() if k.startswith('gates.prx_12')})

# %%
compiler.print_implementations_trees(compiler.builder.op_table["measure"])


