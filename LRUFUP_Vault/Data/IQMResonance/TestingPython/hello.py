from qiskit import QuantumCircuit
from qiskit import transpile
from iqm.qiskit_iqm import IQMFakeGarnet
import matplotlib.pyplot as plt

num_qb = 11
qc = QuantumCircuit(num_qb)

qc.h(0)
qc.cx(0, 1)
qc.cx(0, 2)
qc.cx(0, 3)
qc.cx(0, 4)
qc.cx(1, 5)
qc.cx(1, 6)
qc.cx(1, 7)
qc.cx(2, 8)
qc.cx(2, 9)
qc.cx(2, 10)
qc.barrier()


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

qc.measure_all()
print(qc.draw(output='text'))