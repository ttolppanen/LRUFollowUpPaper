from qiskit import QuantumCircuit, transpile
from iqm.qiskit_iqm.fake_backends.fake_garnet import IQMFakeGarnet
import matplotlib.pyplot as plt

num_qb = 11
qc = QuantumCircuit(num_qb, num_qb)

qc.x(1)
qc.barrier()
qc.measure(range(2), range(2))

backend = IQMFakeGarnet()

qc_transpiled = transpile(qc, backend)
job = backend.run(qc_transpiled, shots=1000)


res=job.result()
counts=res.get_counts()

#Plot the results
from qiskit.visualization import plot_histogram
import numpy as np

fig, ax = plt.subplots(figsize=(20, 6))
plot_histogram(counts, ax=ax)

desired_ticks = 2**6   ### decreasing the total number of tick labels to have a clean figure
x_ticks = np.linspace(0, len(counts)-1, desired_ticks, dtype=int)
ax.set_xticks(x_ticks)

plt.show()

print(qc.draw(output='text'))