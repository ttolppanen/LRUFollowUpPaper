### Abstract
We look at leakage mobility in tunable couplers. We find that in the $g_{eff} = 0$ position the coupler halts qubit excitation swapping, but does not halt leakage swapping. We show this in detail, and we show the conditions in which this happens. We also look at leakage tunneling, and find that it can be relevant. We demonstrate numerically how tunable couplers could provide a passive leakage removal method. Lastly, we apply similar ideas to the junction readout scheme. We find that junction readout offers natural dissipation, which can be made to target only the $|f\rangle$ level of the system.
# Introduction
- Why tunable couplers are relevant.
- Why is leakage relevant?
- Why are we interested in the leakage dynamics/passive propagation of leakage here?
- Summary of the chapters?
# Model
- Introduce the Hamiltonian for the system.
- Mention that we use the coupler below setup.
- Mention that we neglect the $\hat{a}\hat{a}$ etc. term. We do not see much difference between having them or not. Not having allows numerics to be done in total boson number subspaces.
- We focus on the leakage dynamics in the sense of the $|f\rangle$ state, not with leakage in the coupler.
- We also look at the tunneling. Here the middle transmon is higher in resonance. 
- We also look with adding dissipation to the last site.
# Leakage Dynamics with Tunable Couplers
- We can have leakage mobility while canceling the qubit interaction, if the two transmons are on resonance and the tunable coupler is tuned in the $g_{eff} = 0$ position.
- The resonance of the two transmons is crucial.
- We assume that the sites have the same anharmonicity. In practice this is not true, but it does not seem to matter too much.
- If the coupler is tuned in the $g_{zz} = 0$ position, the resonance condition is more strict.
- We get faster dynamics if the target site is populated, $|201\rangle$.
- Here also the requirements for U and resonance are a lot more relaxed.
- But in the $g_{zz} = 0$ position the resonance is more strict, and since here the qubits are operated with large detuning, in practice you would not expect leakage propagation.
- I think it would be nice if we could plot the coupler frequencies where $g_{eff}=0$, $g_{zz}=0$ and $g_{leakage}=0$. 
# Leakage Tunneling
- Normally the qubits are operated in the $g_{zz} = 0$ position where they are very detuned. Due to this there can be structure, so that the next nearest neighbors can be close to resonance.
- We find that in specific cases this resonance is enough for there to be leakage tunneling. The conditions for this are:
	- In the three transmon setup the middle transmon should be higher in frequency than the left and right one. Here then the couplers are closer in frequency to the left and right one, compared to the case where the middle transmon is lower in frequency than left and right.
	- The target site, the right site, should be populated. Population in the middle site can make the dynamics faster, but if the right site is empty, we do not see leakage hopping.
	- The anharmonicity does not matter too much, we saw that it can make the dynamics faster or slower, but the point is that this leakage tunneling still happens.
- Larger detuning from the middle transmon causes faster leakage mobility. Possibly because closer in resonance the right transmon can give its single excitation to the middle one, which blocks the leakage from propagating.
- The $|211\rangle$ state is less strict about the resonance condition than $|201\rangle$, but they are both still quite sensitive. 
- The conclusion is that leakage tunneling can be relevant in tunable coupler setups.
- It causes passive leakage propagation, so the only source of leakage mobility/propagation is not gates, but also tunneling.
- The tunneling can be difficult to avoid in practice, since each qubit has eight nearest neighbors.
# Passive Leakage Removal
- Since the tunable coupler is capable of passing leakage while blocking qubit excitations, it is a natural candidate for a passive LRU.
- We turn the tunable coupler setup to a passive leakage removal unit by adding dissipation to the other transmon. The dissipation is added so that we get the fastest leakage removal. Too strong dissipation also hinders leakage removal, due to zeno-effect.
- We benchmark the performance of it, by plotting the leakage removal and qubit population together with what we would get if we would just capacitively couple two transmons (no coupler in between, $|20\rangle \rightarrow |11\rangle$ in resonance).
- There is also the case where we have $|200\rangle \rightarrow |101\rangle$ in resonance with the coupler, but somehow I feel like maybe this does not need to be in the paper?
- Leakage removal using the $|200\rangle \rightarrow |002\rangle$ seems rather slow.
- Leakage removal using the $|201\rangle \rightarrow |102\rangle$ by pumping the second transmon seems fast, almost as fast as the direct $|20\rangle \rightarrow |11\rangle$ case.
- In the qubit lifetimes the $|20\rangle \rightarrow |11\rangle$ performs badly, and all the other cases quite well. The pumped case has very minor qubit loss.
- In the pumped case, in short times there is additional population going in to the coupler. This should not be a problem, since the system is just going in the logical $|11\rangle$ state. But can it explain why there is only little loss, since it repopulates the coupler little bit?
- Anyways, the pumped removal seems to work well, it removes leakage quite fast, and does not remove the qubit. But since we assume resonance, and they are neighbors, in practice cross-talk might be a problem.
# Passive Leakage Removal via Junction Readout
- Leakage removal with tunable coupled transmons implies that you have to either sacrifice existing transmons in a QPU for leakage removal, or you would have to engineer additional transmon in the QPU. Here we look at a LRU that takes advantage of already existing readout. This then would not cause additional connectivity or complexity in the QPU.
- We adopt a junction readout scheme, where a readout resonator and a transmon qubit are coupled capacitively and through a Josephson junction. The capacitance (?) is chosen specifically so that $g_{eff} = 0$ between the transmon and the resonator. Then we have a similar setup as discussed previously in the paper.
- We can tune the transmon away from the readout frequency, so that the $|f\rangle$ levels of the transmon and resonator are on resonance.
- Then we see that there is leakage propagation, while the qubit propagation is minimal, due to the choice of the capacitive coupling.
- The readout resonator has intrinsic dissipation, giving the source of leakage removal for a passive LRU.
- The leakage removal time is on par with the previously shown passive leakage removal methods.
- The two-mode squeezing can be a problem here, it causes error in the first excited state, and even in the ground state. For the ground state it causes evolution $|00\rangle \rightarrow |11\rangle$. Dissipation breaks this coherent evolution and causes population build up in the ground state. For the first excited state $|10\rangle$, the coupling is not fully blocked, and there can also be evolution to $|21\rangle$. Due to these the dissipation causes qubit population loss. But when compared with realistic $T_1$ on the qubit, it does not cause too much extra dissipation on it.
- Therefore, in principle we could use a junction readout setup as passive LRU, by tuning the qubit so that the second levels are on resonance between it and the resonator. Since every qubit is connected to a readout resonator, every qubit would have available LRU.